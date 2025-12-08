#----------- 1. Create Cross-Validation Datasets ----------------
# Leave-One-Out CV without buffer exclusion
loo_cv <- function(spatial_data) {
  
  # Initialize a list to store train-test splits
  train_test_splits <- vector("list", length = nrow(spatial_data))
  
  # Loop through each observation
  for (i in 1:nrow(spatial_data)) {
    
    # Define test set (leave one point out)
    test_set <- spatial_data[i, , drop = FALSE]
    
    # Define train set (all other points)
    train_set <- spatial_data[-i, ]
    
    # Store train-test split
    train_test_splits[[i]] <- list(train = train_set, test = test_set)
  }
  
  return(train_test_splits)
}








# Leave-One-Out CV with buffer exclusion
bloo_cv <- function(spatial_data, buffer_size_km) {
  
  # Convert buffer size to meters
  buffer_size_m <- buffer_size_km * 1000
  
  # Initialize a list to store train-test splits
  train_test_splits <- list()
  
  # Extract coordinates
  coords <- st_coordinates(spatial_data)  # Extract spatial coordinates
  
  # Loop through each observation
  for (i in 1:nrow(spatial_data)) {
    
    # Define the test point (one point per iteration)
    test_set <- spatial_data[i, , drop = FALSE]  
    
    # Compute distances from the test point to all other points
    distances <- distVincentySphere(
      p1 = coords, 
      p2 = coords[i, , drop = FALSE]  # Reference test point
    )
    
    # Identify training points (exclude points within the buffer)
    train_indices <- which(distances > buffer_size_m)  # Keep only points outside the buffer
    
    # Define training set (excluding the test point and buffer)
    train_set <- spatial_data[train_indices, ]
    
    # Store train-test split
    train_test_splits[[i]] <- list(train = train_set, test = test_set)
  }
  
  return(train_test_splits)
}




























#----------- 2. Fit Models ----------------
fit_models <- function(dataset, response_var, predictors, cv_configs = NULL,
                       models = c("GLM", "GAM", "RF", "GLM.nb", "spaMM", "XGB", "GLSRF"),
                       distribution = NULL, mtry = NULL) {
  
  library(RandomForestsGLS)
  
  cat("\n -------------------","\n Starting Model Evaluation for:", response_var, "\n")
  
  fitted_models <- list()
  
  use_distribution <- any(models %in% c("GLM", "GAM", "spaMM", "GLM.nb"))
  
  if (use_distribution && is.null(distribution)) {
    stop("The 'distribution' parameter is required for GLM, GAM, GLM.nb, and spaMM models!")
  }
  
  if (!is.null(distribution)){
    d <- distribution$Distribution[distribution$response_var == response_var]  
    
    if (length(d) != 1) {
      stop(paste("Error: Found", length(d), "distribution values for", response_var, "instead of one."))
    }
    select_family <- function(d, gam_call = FALSE, spamm_call = FALSE) {
      if (d == "poisson") {
        return(poisson(link = "log"))
      } else if (d == "nb") {
        if (spamm_call) return(negbin2())
        if (gam_call) return(nb())
        return(MASS::negative.binomial(theta = 1))
      } else if (d == "quasipoisson") {
        return(quasipoisson(link = "log"))
      } else if (d == "gaussian") {
        return(gaussian())
      } else {
        stop(paste("Unknown distribution type:", d))
      }
    }
    
  }
  
  gam_call <- "GAM" %in% models
  spamm_call <- "spaMM" %in% models
  family <- if (use_distribution) select_family(d, gam_call, spamm_call) else NULL
  
  for (model_name in models) {
    fitted_models[[model_name]] <- list()
  }
  
  for (cv_name in names(cv_configs)) {
    cat("\n🔹 Running", cv_name, "Cross-Validation\n")
    
    train_test_splits <- cv_configs[[cv_name]]
    
    for (model_name in models) {
      fitted_models[[model_name]][[cv_name]] <- list()
    }
    
    for (i in seq_along(train_test_splits)) {
      cat("   ➡️ Fold:", i, "/", length(train_test_splits), "\n")
      
      train_set <- train_test_splits[[i]]$train
      test_set <- train_test_splits[[i]]$test
      
      train_df <- as.data.frame(st_drop_geometry(train_set))
      test_df  <- as.data.frame(st_drop_geometry(test_set))
      
      if (!(response_var %in% colnames(train_df))) {
        stop(paste("Response variable", response_var, "is missing!"))
      }
      
      formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
      
      # --- GLM ---
      if ("GLM" %in% models) {
        cat("\nFitting GLM for fold", i, "...")
        glm_model <- tryCatch(
          glm(formula = formula, family = family, data = train_df),
          error = function(e) {
            cat("\n GLM failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(glm_model)) {
          predictions <- predict(glm_model, newdata = test_df, type = "response")
          aic_value <- AIC(glm_model)
          fitted_models$GLM[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = glm_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- GAM ---
      if ("GAM" %in% models) {
        cat("\nFitting GAM for fold", i, "...")
        gam_model <- tryCatch(
          mgcv::gam(formula = formula, family = family, data = train_df, method = 'ML'),
          error = function(e) {
            cat("\n GAM failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(gam_model)) {
          predictions <- predict(gam_model, newdata = test_df, type = "response")
          aic_value <- AIC(gam_model)
          fitted_models$GAM[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = gam_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- GLM.nb ---
      if ("GLM.nb" %in% models) {
        cat("\nFitting GLM.nb for fold", i, "...")
        glm_nb_model <- tryCatch(
          MASS::glm.nb(formula = formula, data = train_df, init.theta = 1, control = glm.control(maxit = 50)),
          error = function(e) {
            cat("\n GLM.nb failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(glm_nb_model)) {
          predictions <- predict(glm_nb_model, newdata = test_df, type = "response")
          aic_value <- AIC(glm_nb_model)
          fitted_models$GLM.nb[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = glm_nb_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- RF ---
      if ("RF" %in% models) {
        cat("\nFitting RF for fold", i, "...")
        rf_model <- tryCatch(
          randomForest(formula = formula, data = train_df, ntree = 500, importance = TRUE, mtry = 2, nodesize = 1),
          error = function(e) {
            cat("\n RF failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(rf_model)) {
          predictions <- predict(rf_model, newdata = test_df)
          fitted_models$RF[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = rf_model,
            predictions = predictions,
            AIC = NA,
            test_set = test_set
          )
        }
      }
      
      
      # --- GLS RF ---
      if ("GLSRF" %in% models) {
        cat("\nFitting GLSRF for fold", i, "...")
        
        # Coordinates: as matrix (X, Y)
        coords_train <- sf::st_coordinates(train_set)
        coords_test  <- sf::st_coordinates(test_set)
        
        # Covariate matrix (no geometry, only predictors)
        X_train <- as.matrix(train_df[, predictors, drop = FALSE])
        X_test  <- as.matrix(test_df[, predictors, drop = FALSE])
        
        # Response vector
        y_train <- train_df[[response_var]]
        
        # Fit model 
        glsrf_model <- RandomForestsGLS::RFGLS_estimate_spatial(
          coords = coords_train,
            y = y_train,
            X = X_train,
            #ntree = 500, # Default is 50
            mtry = 2
        )
        
        if (!is.null(glsrf_model)) {
          # Mean-function prediction (ignores spatial dependence)
          glsrf_pred_mean <- RandomForestsGLS::RFGLS_predict(
            RFGLS_out = glsrf_model,
            Xtest     = X_test
          )
          
          # Spatial prediction (uses coords + covariates)
          glsrf_pred_spatial <- RandomForestsGLS::RFGLS_predict_spatial(
            RFGLS_out = glsrf_model,
            coords.0  = coords_test,
            Xtest     = X_test
          )
          
          # By default, store the spatial response as 'predictions'
          # but keep the mean prediction as well if you want to compare later
          fitted_models$GLSRF[[cv_name]][[paste0("Fold_", i)]] <- list(
            model              = glsrf_model,
            predictions        = glsrf_pred_spatial$prediction,  # spatial response
            predictions_mean   = glsrf_pred_mean$predicted,      # mean function
            AIC                = NA,
            test_set           = test_set
          )
        }
      }
      
      
      
      # --- XGBoost ---
      if ("XGB" %in% models) {
        cat("\nFitting XGBoost for fold", i, "...")
        
        # Convert character columns to factors
        train_df <- train_df %>% mutate(across(where(is.character), as.factor))
        test_df  <- test_df %>% mutate(across(where(is.character), as.factor))
        
        # Drop unused levels
        train_x_raw <- droplevels(train_df[, predictors])
        test_x_raw  <- droplevels(test_df[, predictors])
        
        # Create dummy-encoded numeric matrices
        train_x <- model.matrix(~ . - 1, data = train_x_raw)
        test_x  <- model.matrix(~ . - 1, data = test_x_raw)
        
        # Build DMatrix
        train_matrix <- tryCatch({
          xgb.DMatrix(data = train_x, label = train_df[[response_var]])
        }, error = function(e) {
          cat("\n XGB train matrix creation failed:", e$message, "\n")
          return(NULL)
        })
        
        test_matrix <- tryCatch({
          xgb.DMatrix(data = test_x, label = test_df[[response_var]])
        }, error = function(e) {
          cat("\n XGB test matrix creation failed:", e$message, "\n")
          return(NULL)
        })
        
        # Train and predict
        if (!is.null(train_matrix) && !is.null(test_matrix)) {
          xgb_model <- tryCatch(
            xgb.train(
              params = list(
                objective = "reg:squarederror",
                eval_metric = "rmse",
                eta = 0.1,
                max_depth = 6
              ),
              data = train_matrix,
              watchlist = list(train = train_matrix, test = test_matrix), # !!!!!!!! ERROR DATA LEACKAGE !!!!!!!!!
              nrounds = 100,
              early_stopping_rounds = 10,
              verbose = 0
            ),
            error = function(e) {
              cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
              return(NULL)
            }
          )
          
          if (!is.null(xgb_model)) {
            predictions <- predict(xgb_model, newdata = test_matrix)
            fitted_models$XGB[[cv_name]][[paste0("Fold_", i)]] <- list(
              model = xgb_model,
              predictions = predictions,
              AIC = NA,
              test_set = test_set
            )
          }
        }
      }
    }
  }
  
  return(fitted_models)
}


fit_models_eval <- function(dataset, response_var, predictors, cv_configs = NULL,
                            models = c("GLM", "GAM", "RF", "GLM.nb", "spaMM", "XGB"),
                            distribution = NULL, mtry =if (!is.null(y) && !is.factor(y))
                              max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))) {
  
  cat("\n -------------------","\n▶ Starting Model Evaluation for:", response_var, "\n")
  
  fitted_models <- list()
  
  use_distribution <- any(models %in% c("GLM", "GAM", "spaMM", "GLM.nb"))
  
  if (use_distribution && is.null(distribution)) {
    stop("❌ The 'distribution' parameter is required for GLM, GAM, GLM.nb, and spaMM models!")
  }
  
  if (!is.null(distribution)){
    d <- distribution$Distribution[distribution$response_var == response_var]  
    
    if (length(d) != 1) {
      stop(paste("❌ Error: Found", length(d), "distribution values for", response_var, "instead of one."))
    }
    select_family <- function(d, gam_call = FALSE, spamm_call = FALSE) {
      if (d == "poisson") {
        return(poisson(link = "log"))
      } else if (d == "nb") {
        if (spamm_call) return(negbin2())
        if (gam_call) return(nb())
        return(MASS::negative.binomial(theta = 1))
      } else if (d == "quasipoisson") {
        return(quasipoisson(link = "log"))
      } else if (d == "gaussian") {
        return(gaussian())
      } else {
        stop(paste("❌ Unknown distribution type:", d))
      }
    }
    
  }
  
  gam_call <- "GAM" %in% models
  spamm_call <- "spaMM" %in% models
  family <- if (use_distribution) select_family(d, gam_call, spamm_call) else NULL
  
  for (model_name in models) {
    fitted_models[[model_name]] <- list()
  }
  
  for (cv_name in names(cv_configs)) {
    cat("\n🔹 Running", cv_name, "Cross-Validation\n")
    
    train_test_splits <- cv_configs[[cv_name]]
    
    for (model_name in models) {
      fitted_models[[model_name]][[cv_name]] <- list()
    }
    
    for (i in seq_along(train_test_splits)) {
      cat("   ➡️ Fold:", i, "/", length(train_test_splits), "\n")
      
      train_set <- train_test_splits[[i]]$train
      test_set <- train_test_splits[[i]]$test
      eval_set <- train_test_splites[[i]]$eval
      
      train_df <- as.data.frame(st_drop_geometry(train_set))
      test_df  <- as.data.frame(st_drop_geometry(test_set))
      eval_df <- as.data.frame(st_drop_geometry(eval_set))
      
      if (!(response_var %in% colnames(train_df))) {
        stop(paste("❌ Response variable", response_var, "is missing!"))
      }
      
      formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
      
      # --- GLM ---
      if ("GLM" %in% models) {
        cat("\nFitting GLM for fold", i, "...")
        glm_model <- tryCatch(
          glm(formula = formula, family = family, data = train_df),
          error = function(e) {
            cat("\n❌ GLM failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(glm_model)) {
          predictions <- predict(glm_model, newdata = test_df, type = "response")
          aic_value <- AIC(glm_model)
          fitted_models$GLM[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = glm_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- GAM ---
      if ("GAM" %in% models) {
        cat("\nFitting GAM for fold", i, "...")
        gam_model <- tryCatch(
          mgcv::gam(formula = formula, family = family, data = train_df, method = 'ML'),
          error = function(e) {
            cat("\n❌ GAM failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(gam_model)) {
          predictions <- predict(gam_model, newdata = test_df, type = "response")
          aic_value <- AIC(gam_model)
          fitted_models$GAM[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = gam_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- GLM.nb ---
      if ("GLM.nb" %in% models) {
        cat("\nFitting GLM.nb for fold", i, "...")
        glm_nb_model <- tryCatch(
          MASS::glm.nb(formula = formula, data = train_df, init.theta = 1, control = glm.control(maxit = 50)),
          error = function(e) {
            cat("\n❌ GLM.nb failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(glm_nb_model)) {
          predictions <- predict(glm_nb_model, newdata = test_df, type = "response")
          aic_value <- AIC(glm_nb_model)
          fitted_models$GLM.nb[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = glm_nb_model,
            predictions = predictions,
            AIC = aic_value,
            test_set = test_set
          )
        }
      }
      
      # --- RF ---
      if ("RF" %in% models) {
        cat("\nFitting RF for fold", i, "...")
        rf_model <- tryCatch(
          randomForest(formula = formula, data = train_df, ntree = 500, mtry = mtry, importance = TRUE),
          error = function(e) {
            cat("\n RF failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(rf_model)) {
          predictions <- predict(rf_model, newdata = test_df)
          fitted_models$RF[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = rf_model,
            predictions = predictions,
            AIC = NA,
            test_set = test_set
          )
        }
      }
      
      # --- XGBoost ---
      if ("XGB" %in% models) {
        cat("\nFitting XGBoost for fold", i, "...")
        
        # Convert character columns to factors
        train_df <- train_df %>% mutate(across(where(is.character), as.factor))
        test_df  <- test_df %>% mutate(across(where(is.character), as.factor))
        eval_df <- eval_df %>% mutate(across(where(is.character), as.factor))
        
        # Drop unused levels
        train_x_raw <- droplevels(train_df[, predictors])
        test_x_raw  <- droplevels(test_df[, predictors])
        eval_x_raw  <- droplevels(eval_df[, predictors])
        
        # Create dummy-encoded numeric matrices
        train_x <- model.matrix(~ . - 1, data = train_x_raw)
        test_x  <- model.matrix(~ . - 1, data = test_x_raw)
        
        # Build DMatrix
        train_matrix <- 
          xgb.DMatrix(data = train_x, label = train_df[[response_var]])
        
        
        test_matrix <- 
          xgb.DMatrix(data = test_x, label = test_df[[response_var]])
        
        eval_matrix <- 
          xgb.DMatrix(data = test_x, label = eval_df[[response_var]])
        
        
        
        # Train and predict
        if (!is.null(train_matrix) && !is.null(test_matrix)) {
          xgb_model <- tryCatch(
            xgb.train(
              params = list(
                objective = "reg:squarederror",
                eval_metric = "rmse",
                eta = 0.1,
                max_depth = 6
              ),
              data = train_matrix,
              watchlist = list(train = train_matrix, test = eval_matrix),
              nrounds = 100,
              early_stopping_rounds = 10,
              verbose = 0
            ),
            error = function(e) {
              cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
              return(NULL)
            }
          )
          
          if (!is.null(xgb_model)) {
            predictions <- predict(xgb_model, newdata = test_matrix)
            fitted_models$XGB[[cv_name]][[paste0("Fold_", i)]] <- list(
              model = xgb_model,
              predictions = predictions,
              AIC = NA,
              test_set = test_set
            )
          }
        }
      }
    }
  }
  
  return(fitted_models)
}






tune_xgb <- function(dataset, response_var, 
                     predictors, 
                     cv_configs = NULL,
                     models = "XGB", 
                     params = list(
                       objective = "reg:squarederror",
                       eval_metric = "rmse",
                       eta = 0.1,
                       max_depth = 6
                     ), 
                     nrounds = 100) {
  
  cat("\n -------------------","\n▶ Starting Model Evaluation for:", response_var, "\n")
  
  fitted_models <- list()
  
  
  for (model_name in models) {
    fitted_models[[model_name]] <- list()
  }
  
  for (cv_name in names(cv_configs)) {
    cat("\n🔹 Running", cv_name, "Cross-Validation\n")
    
    train_test_splits <- cv_configs[[cv_name]]
    
    for (model_name in models) {
      fitted_models[[model_name]][[cv_name]] <- list()
    }
    
    for (i in seq_along(train_test_splits)) {
      cat("   ➡️ Fold:", i, "/", length(train_test_splits), "\n")
      
      train_set <- train_test_splits[[i]]$train
      test_set <- train_test_splits[[i]]$test
      
      train_df <- as.data.frame(st_drop_geometry(train_set))
      test_df  <- as.data.frame(st_drop_geometry(test_set))
      
      if (!(response_var %in% colnames(train_df))) {
        stop(paste("❌ Response variable", response_var, "is missing!"))
      }
      
      formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
      
      
      # --- XGBoost ---
      cat("\nFitting XGBoost for fold", i, "...")
      
      # Convert character columns to factors
      train_df <- train_df %>% mutate(across(where(is.character), as.factor))
      test_df  <- test_df %>% mutate(across(where(is.character), as.factor))
      
      # Drop unused levels
      train_x_raw <- droplevels(train_df[, predictors])
      test_x_raw  <- droplevels(test_df[, predictors])
      
      # Create dummy-encoded numeric matrices
      train_x <- model.matrix(~ . - 1, data = train_x_raw)
      test_x  <- model.matrix(~ . - 1, data = test_x_raw)
      
      # Build DMatrix
      train_matrix <- tryCatch({
        xgb.DMatrix(data = train_x, label = train_df[[response_var]])
      }, error = function(e) {
        cat("\n❌ XGB train matrix creation failed:", e$message, "\n")
        return(NULL)
      })
      
      test_matrix <- tryCatch({
        xgb.DMatrix(data = test_x, label = test_df[[response_var]])
      }, error = function(e) {
        cat("\n❌ XGB test matrix creation failed:", e$message, "\n")
        return(NULL)
      })
      
      # Train and predict
      if (!is.null(train_matrix) && !is.null(test_matrix)) {
        xgb_model <- tryCatch(
          xgb.train(
            params = params,
            data = train_matrix,
            watchlist = list(train = train_matrix, test = test_matrix),
            nrounds = nrounds,
            early_stopping_rounds = 10,
            verbose = 0
          ),
          error = function(e) {
            cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
            return(NULL)
          }
        )
        
        if (!is.null(xgb_model)) {
          predictions <- predict(xgb_model, newdata = test_matrix)
          fitted_models$XGB[[cv_name]][[paste0("Fold_", i)]] <- list(
            model = xgb_model,
            predictions = predictions,
            AIC = NA,
            test_set = test_set
          )
        }
      }
      
    }
  }
  
  return(fitted_models)
}



#----------- 3. Evaluate Models ----------------
evaluate_models <- function(models_list,
                            response_var = NULL,
                            models       = NULL,
                            cv_names     = NULL,   # names of CV methods to use (e.g. "env_k6", "bloo50")
                            cv_splits    = NULL,   # optional: list of CV splits to compute train sizes
                            output       = c("summary", "fold", "plot", "all")) {
  
  output <- match.arg(output)
  cat("\n* Starting Model Evaluation\n")
  
  #-----------------------------------------
  # Helper to compute metrics
  #-----------------------------------------
  compute_metrics <- function(observed, predicted, aic_values = NULL) {
    data.frame(
      R2            = cor(observed, predicted, use = "complete.obs")^2,
      Pearson_Corr  = cor(observed, predicted, method = "pearson",  use = "complete.obs"),
      Spearman_Corr = cor(observed, predicted, method = "spearman", use = "complete.obs"),
      MAE           = mean(abs(observed - predicted), na.rm = TRUE),
      RMSE          = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
      AIC           = if (!is.null(aic_values)) mean(aic_values, na.rm = TRUE) else NA_real_
    )
  }
  
  # Storage
  performance_summary  <- data.frame()
  performance_per_fold <- data.frame()
  plot_list            <- list()
  
  # Filter response variables
  response_vars <- if (!is.null(response_var)) {
    intersect(response_var, names(models_list))
  } else {
    names(models_list)
  }
  
  for (resp_var in response_vars) {
    
    available_models <- names(models_list[[resp_var]])
    
    model_types <- if (!is.null(models)) {
      intersect(models, available_models)
    } else {
      available_models
    }
    
    for (model_name in model_types) {
      
      cat("\n* Evaluating model:", model_name, "for response variable:", resp_var, "\n")
      
      available_cvs <- names(models_list[[resp_var]][[model_name]])
      
      cv_methods <- if (!is.null(cv_names)) {
        intersect(cv_names, available_cvs)
      } else {
        available_cvs
      }
      
      for (cv_name in cv_methods) {
        
        cat("\n* Processing CV Method:", cv_name, "\n")
        
        model_folds <- models_list[[resp_var]][[model_name]][[cv_name]]
        # Optional CV splits to get train sizes
        this_cv_split <- if (!is.null(cv_splits)) cv_splits[[cv_name]] else NULL
        
        all_observed    <- c()
        all_predictions <- c()
        aic_values      <- c()
        train_sizes     <- c()
        
        for (i in seq_along(model_folds)) {
          
          cat("   - Fold:", i, "/", length(model_folds), "\n")
          
          fold_name <- paste0("Fold_", i)
          fold_data <- model_folds[[fold_name]]
          
          if (is.null(fold_data)) {
            cat("      Model missing for", fold_name, "- Skipping.\n")
            next
          }
          
          test_set <- fold_data$test_set
          
          # Allow multi-row test sets; skip only if empty
          if (is.null(test_set) || nrow(test_set) == 0) {
            cat("      Test set missing or empty for", fold_name, "- Skipping.\n")
            next
          }
          
          test_df <- as.data.frame(st_drop_geometry(test_set))
          
          if (!(resp_var %in% colnames(test_df))) {
            stop(paste("Error:", fold_name, "- Response variable missing in test set!"))
          }
          
          fold_observed <- test_df[[resp_var]]
          
          #-----------------------------
          # Predictions
          #-----------------------------
          fold_pred <- tryCatch({
            if (model_name == "XGB") {
              test_df <- test_df %>% dplyr::mutate(across(where(is.character), as.factor))
              
              xgb_features   <- fold_data$model$feature_names
              available_cols <- intersect(colnames(test_df), xgb_features)
              
              test_x_raw <- test_df[, available_cols, drop = FALSE]
              
              valid_cols <- sapply(test_x_raw, function(col) {
                if (is.factor(col)) {
                  nlevels(col) > 1
                } else if (is.numeric(col)) {
                  length(unique(col)) > 1
                } else {
                  TRUE
                }
              })
              
              test_x_clean <- test_x_raw[, valid_cols, drop = FALSE]
              
              if (ncol(test_x_clean) == 0) {
                test_matrix <- matrix(0, nrow = nrow(test_df), ncol = length(xgb_features))
                colnames(test_matrix) <- xgb_features
              } else {
                test_matrix <- model.matrix(~ . - 1, data = test_x_clean)
                
                missing_cols <- setdiff(xgb_features, colnames(test_matrix))
                for (col in missing_cols) {
                  test_matrix <- cbind(test_matrix, setNames(data.frame(0), col))
                }
                
                test_matrix <- test_matrix[, xgb_features, drop = FALSE]
              }
              
              test_matrix <- xgb.DMatrix(data = test_matrix)
              predict(fold_data$model, newdata = test_matrix)
              
            } else if (model_name == "RF") {
              predict(fold_data$model, newdata = test_df)
              
            } else {
              predict(fold_data$model, newdata = test_df, type = "response")
            }
          }, error = function(e) {
            cat("      Prediction failed for", model_name, "on", fold_name, ":", e$message, "\n")
            return(rep(NA_real_, length(fold_observed)))
          })
          
          # Append to global vectors
          all_observed    <- c(all_observed,    fold_observed)
          all_predictions <- c(all_predictions, fold_pred)
          
          # AIC per fold if relevant
          if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM")) {
            aic_values <- c(aic_values, AIC(fold_data$model))
          }
          
          # Train size (if we have cv_splits)
          this_train_size <- NA_real_
          if (!is.null(this_cv_split) &&
              length(this_cv_split) >= i &&
              !is.null(this_cv_split[[i]]$train)) {
            this_train_size <- nrow(this_cv_split[[i]]$train)
            train_sizes     <- c(train_sizes, this_train_size)
          }
          
          #-----------------------------
          # Per-fold metrics
          #-----------------------------
          if (!all(is.na(fold_pred))) {
            fold_metrics <- compute_metrics(
              fold_observed, fold_pred,
              if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM"))
                AIC(fold_data$model) else NULL
            )
            
            fold_metrics$Response_Var <- resp_var
            fold_metrics$Model        <- model_name
            fold_metrics$CV           <- cv_name
            fold_metrics$Fold         <- i
            fold_metrics$Train_Size   <- this_train_size
            
            performance_per_fold <- rbind(performance_per_fold, fold_metrics)
            
            cat(sprintf(
              "      Fold %d metrics: R2 = %.3f | RMSE = %.3f | MAE = %.3f | Pearson = %.3f | Spearman = %.3f\n",
              i,
              fold_metrics$R2,
              fold_metrics$RMSE,
              fold_metrics$MAE,
              fold_metrics$Pearson_Corr,
              fold_metrics$Spearman_Corr
            ))
          } else {
            cat("      All predictions NA for", fold_name, "- no fold metrics.\n")
          }
        } # end folds
        
        #-----------------------------
        # Global (across folds) metrics
        #-----------------------------
        if (length(all_observed) == 0 || length(all_predictions) == 0) {
          cat("\nNo valid predictions for model", model_name, "and CV", cv_name, "- skipping global metrics.\n")
        } else {
          global_metrics <- compute_metrics(all_observed, all_predictions, aic_values)
          global_metrics$Response_Var <- resp_var
          global_metrics$Model        <- model_name
          global_metrics$CV           <- cv_name
          global_metrics$Train_Size   <- if (length(train_sizes) > 0) mean(train_sizes) else NA_real_
          
          performance_summary <- rbind(performance_summary, global_metrics)
          
          #-----------------------------
          # Plot observed vs predicted
          #-----------------------------
          if (output %in% c("plot", "all")) {
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              warning("ggplot2 not available; plots will not be created.")
            } else {
              df_plot <- data.frame(
                Observed  = all_observed,
                Predicted = all_predictions
              )
              
              # build annotation label from global metrics
              metrics_label <- sprintf(
                "R2 = %.2f\nPearson = %.2f\nSpearman = %.2f",
                global_metrics$R2,
                global_metrics$Pearson_Corr,
                global_metrics$Spearman_Corr
              )
              
              p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Observed, y = Predicted)) +
                ggplot2::geom_point(alpha = 0.6) +
                ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                ggplot2::annotate(
                  "text",
                  x     = Inf,
                  y     = -Inf,
                  hjust = 1.1,
                  vjust = -0.5,
                  label = metrics_label,
                  size  = 3
                ) +
                ggplot2::labs(
                  title = paste("Observed vs Predicted -", resp_var, "-", model_name, "-", cv_name),
                  x = "Observed",
                  y = "Predicted"
                ) +
                ggplot2::theme_minimal()
              
              plot_name <- paste(resp_var, model_name, cv_name, sep = "_")
              plot_list[[plot_name]] <- p
            }
          }
        } # end global metrics
        
      } # end CV methods
    } # end model types
  } # end response vars
  
  #-----------------------------
  # Return according to 'output'
  #-----------------------------
  if (output == "summary") {
    return(performance_summary)
  } else if (output == "fold") {
    return(performance_per_fold)
  } else if (output == "plot") {
    return(plot_list)
  } else {  # "all"
    return(list(
      summary  = performance_summary,
      per_fold = performance_per_fold,
      plots    = plot_list
    ))
  }
}


evaluate_models_170425 <- function(models_list, response_var = NULL, models = NULL, cv_configs = NULL) {
  
  cat("\n🔹 Starting Model Evaluation\n")
  
  performance_metrics <- data.frame()
  
  # Determine which response variables to evaluate
  response_vars <- if (!is.null(response_var)) {
    intersect(response_var, names(models_list))
  } else {
    names(models_list)
  }
  
  for (resp_var in response_vars) {
    
    available_models <- names(models_list[[resp_var]])
    
    model_types <- if (!is.null(models)) {
      intersect(models, available_models)
    } else {
      available_models
    }
    
    for (model_name in model_types) {
      
      cat("\n🔹 Evaluating model:", model_name, "for response variable:", resp_var, "\n")
      
      available_cvs <- names(models_list[[resp_var]][[model_name]])
      
      cv_methods <- if (!is.null(cv_configs)) {
        intersect(cv_configs, available_cvs)
      } else {
        available_cvs
      }
      
      for (cv_name in cv_methods) {
        
        cat("\n🔹 Processing CV Method:", cv_name, "\n")
        
        model_folds <- models_list[[resp_var]][[model_name]][[cv_name]]
        cv_config <- cv_configs[[cv_name]]
        
        all_observed <- c()
        predictions <- c()
        aic_values <- c()
        train_size_list <- c()
        
        for (i in seq_along(model_folds)) {
          cat("   ➡️ Fold:", i, "/", length(model_folds), "\n")
          
          fold_data <- model_folds[[paste0("Fold_", i)]]
          
          if (is.null(fold_data) || is.null(fold_data$model)) {
            cat("   ⚠️ Fold", i, ": Missing model — skipping.\n")
            next
          }
          
          model <- fold_data$model
          test_set <- fold_data$test_set
          
          if (is.null(test_set) || nrow(test_set) != 1) {
            cat("   ⚠️ Fold", i, ": Invalid test set — skipping.\n")
            next
          }
          
          test_df <- as.data.frame(sf::st_drop_geometry(test_set))
          
          if (!(resp_var %in% names(test_df))) {
            cat("   ❌ Fold", i, ": Response variable missing in test set!\n")
            next
          }
          
          all_observed <- c(all_observed, test_df[[resp_var]])
          
          # Prediction block
          pred <- tryCatch({
            if (model_name == "XGB") {
              test_df <- test_df |> dplyr::mutate(across(where(is.character), as.factor))
              xgb_features <- model$feature_names
              available_cols <- intersect(names(test_df), xgb_features)
              test_x_raw <- test_df[, available_cols, drop = FALSE]
              
              valid_cols <- sapply(test_x_raw, function(col) {
                if (is.factor(col)) nlevels(col) > 1
                else if (is.numeric(col)) length(unique(col)) > 1
                else TRUE
              })
              
              test_x_clean <- test_x_raw[, valid_cols, drop = FALSE]
              
              if (ncol(test_x_clean) == 0) {
                test_matrix <- matrix(0, nrow = 1, ncol = length(xgb_features))
                colnames(test_matrix) <- xgb_features
              } else {
                test_matrix <- model.matrix(~ . - 1, data = test_x_clean)
                missing_cols <- setdiff(xgb_features, colnames(test_matrix))
                for (col in missing_cols) {
                  test_matrix <- cbind(test_matrix, setNames(data.frame(0), col))
                }
                test_matrix <- test_matrix[, xgb_features, drop = FALSE]
              }
              
              test_matrix <- xgboost::xgb.DMatrix(data = test_matrix)
              predict(model, newdata = test_matrix)
              
            } else if (model_name == "RF") {
              predict(model, newdata = test_df)
              
            } else {
              predict(model, newdata = test_df, type = "response")
            }
          }, error = function(e) {
            cat("   ❌ Prediction failed for fold", i, ":", e$message, "\n")
            return(NA)
          })
          
          predictions <- c(predictions, pred)
          
          if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM")) {
            aic_values <- c(aic_values, AIC(model))
          }
          
          train_size_list <- c(train_size_list, nrow(cv_config[[i]]$train))
        }
        
        # Compute metrics safely
        compute_metrics <- function(observed, predicted, aic_values) {
          if (length(observed) == 0 || length(predicted) == 0 || all(is.na(predicted))) {
            cat("   ⚠️ No valid predictions — skipping metrics.\n")
            return(data.frame(
              R2 = NA,
              Pearson_Corr = NA,
              Spearman_Corr = NA,
              MAE = NA,
              RMSE = NA,
              AIC = if (!is.null(aic_values)) mean(aic_values, na.rm = TRUE) else NA
            ))
          }
          
          data.frame(
            R2 = cor(observed, predicted, use = "complete.obs")^2,
            Pearson_Corr = cor(observed, predicted, method = "pearson", use = "complete.obs"),
            Spearman_Corr = cor(observed, predicted, method = "spearman", use = "complete.obs"),
            MAE = mean(abs(observed - predicted), na.rm = TRUE),
            RMSE = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
            AIC = if (!is.null(aic_values)) mean(aic_values, na.rm = TRUE) else NA
          )
        }
        
        cat("   📊 Evaluating metrics for", length(predictions), "predictions.\n")
        
        model_metrics <- compute_metrics(all_observed, predictions, aic_values)
        model_metrics$Response_Var <- resp_var
        model_metrics$Model <- model_name
        model_metrics$CV <- cv_name
        model_metrics$Train_Size <- if (length(train_size_list) > 0) mean(train_size_list) else NA
        
        performance_metrics <- bind_rows(performance_metrics, model_metrics)
      }
    }
  }
  
  return(performance_metrics)
}





plot_performance_metrics <- function(performance_metrics) {
  
  # Ensure CV and Model are factors for proper ordering
  performance_metrics <- performance_metrics %>%
    mutate(CV = factor(CV, levels = unique(CV)),
           Model = factor(Model, levels = unique(Model)))  # Ensure correct model order
  
  # Define a color scheme for models
  model_colors <- c("GLM" = "#1F78B4", "GLM.nb" = "#A6CEE3", "GAM" = "#33A02C", 
                    "RF" = "#E31A1C", "spaMM" = "#FF7F00")
  
  # Define offsets for Train_Size labels
  text_offset_rmse <- min(performance_metrics$RMSE, na.rm = TRUE) * 0.1
  text_offset_r2 <- min(performance_metrics$R2, na.rm = TRUE) * 0.05
  text_offset_mae <- min(performance_metrics$MAE, na.rm = TRUE) * 0.1
  
  # --- RMSE Plot ---
  p1 <- ggplot(performance_metrics, aes(x = CV, y = RMSE, color = Model)) +
    geom_segment(aes(xend = CV, y = 0, yend = RMSE), linewidth = 1) +  
    geom_point(size = 5, alpha = 0.9) +  
    scale_color_manual(values = model_colors) +
    labs(title = "RMSE Across Cross-Validation Methods", y = "RMSE", x = "Cross-Validation Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +  
    facet_wrap(~ Model) +  
    geom_text(aes(label = round(Train_Size, 0), y = min(RMSE, na.rm = TRUE) - text_offset_rmse), 
              size = 4, color = "black")  
  
  # --- R² Plot ---
  p2 <- ggplot(performance_metrics, aes(x = CV, y = R2, color = Model)) +
    geom_segment(aes(xend = CV, y = 0, yend = R2), linewidth = 1) +  
    geom_point(size = 5, alpha = 0.9) +  
    scale_color_manual(values = model_colors) +
    labs(title = "R² Across Cross-Validation Methods", y = "R²", x = "Cross-Validation Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +  
    facet_wrap(~ Model) +  
    geom_text(aes(label = round(Train_Size, 0), y = min(R2, na.rm = TRUE) - text_offset_r2), 
              size = 4, color = "black")  
  
  # --- MAE Plot ---
  p3 <- ggplot(performance_metrics, aes(x = CV, y = MAE, color = Model)) +
    geom_segment(aes(xend = CV, y = 0, yend = MAE), linewidth = 1) +  
    geom_point(size = 5, alpha = 0.9) +  
    scale_color_manual(values = model_colors) +
    labs(title = "MAE Across Cross-Validation Methods", y = "MAE", x = "Cross-Validation Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +  
    facet_wrap(~ Model) +  
    geom_text(aes(label = round(Train_Size, 0), y = min(MAE, na.rm = TRUE) - text_offset_mae), 
              size = 4, color = "black")  
  
  # Print plots
  print(p1)
  print(p2)
  print(p3)
}


#----------- spBlock_cv : env_xgb_modified ----------------
# source : https://github.com/Disha0903/spBlock_cv/tree/main

env_xgb <- function(hyperparams_xgb, cluster_data, cluster_count, in_time_data, out_of_time_data, model_name) {
  
  results <- data.frame()
  results_past <- data.frame()
  results_past_lf <- data.frame()
  
  feature_names <- colnames(in_time_data)[-ncol(in_time_data)]
  
  for (j in 1:nrow(hyperparams_xgb)) {
    print(paste0("Iteration: ", j))
    fold_ROC_AUC <- c()
    
    for (test_cluster in 1:cluster_count) {
      train_clusters <- setdiff(1:cluster_count, test_cluster)
      
      for (fold in train_clusters) {
        train <- cluster_data[[fold]]
        test  <- cluster_data[[test_cluster]]
        
        X_train <- as.matrix(train[, -ncol(train)])
        colnames(X_train) <- feature_names
        
        xgb_model <- xgboost(
          data = X_train,
          label = train$occurrenceStatus,
          nrounds = hyperparams_xgb$nrounds[j],
          max_depth = hyperparams_xgb$max_depth[j],
          eta = hyperparams_xgb$eta[j],
          subsample = hyperparams_xgb$subsample[j],
          min_child_weight = hyperparams_xgb$min_child_weight[j],
          gamma = hyperparams_xgb$gamma[j],
          colsample_bylevel = hyperparams_xgb$colsample_bylevel[j],
          objective = "binary:logistic",
          eval_metric = "auc",
          verbose = 0
        )
        
        X_test <- as.matrix(test[, -ncol(test)])
        colnames(X_test) <- feature_names
        
        prob_predictions <- predict(xgb_model, newdata = X_test)
        ROC_AUC <- pROC::auc(test$occurrenceStatus, prob_predictions)
        fold_ROC_AUC <- c(fold_ROC_AUC, ROC_AUC)
        
        if (fold == tail(train_clusters, 1)) {
          last_fold_model <- xgb_model
        }
      }
    }
    
    # Save cross-validation results
    results <- rbind(results, data.frame(
      iteration = j,
      nrounds = hyperparams_xgb$nrounds[j],
      max_depth = hyperparams_xgb$max_depth[j],
      eta = hyperparams_xgb$eta[j],
      subsample = hyperparams_xgb$subsample[j],
      min_child_weight = hyperparams_xgb$min_child_weight[j],
      gamma = hyperparams_xgb$gamma[j],
      colsample_bylevel = hyperparams_xgb$colsample_bylevel[j],
      mean_ROC_AUC = mean(fold_ROC_AUC),
      fold_ROC_AUC = toString(round(fold_ROC_AUC, 4))
    ))
    
    # Retrain on all in_time_data
    X_full_train <- as.matrix(in_time_data[, -ncol(in_time_data)])
    colnames(X_full_train) <- feature_names
    
    xgb_model_full <- xgboost(
      data = X_full_train,
      label = in_time_data$occurrenceStatus,
      nrounds = hyperparams_xgb$nrounds[j],
      max_depth = hyperparams_xgb$max_depth[j],
      eta = hyperparams_xgb$eta[j],
      subsample = hyperparams_xgb$subsample[j],
      min_child_weight = hyperparams_xgb$min_child_weight[j],
      gamma = hyperparams_xgb$gamma[j],
      colsample_bylevel = hyperparams_xgb$colsample_bylevel[j],
      objective = "binary:logistic",
      eval_metric = "auc",
      verbose = 0
    )
    
    X_valid <- as.matrix(out_of_time_data[, -ncol(out_of_time_data)])
    colnames(X_valid) <- feature_names
    
    valid_predictions <- predict(xgb_model_full, newdata = X_valid)
    ROC_AUC_valid <- pROC::auc(out_of_time_data$occurrenceStatus, valid_predictions)
    
    results_past <- rbind(results_past, data.frame(
      iteration = j,
      nrounds = hyperparams_xgb$nrounds[j],
      max_depth = hyperparams_xgb$max_depth[j],
      eta = hyperparams_xgb$eta[j],
      subsample = hyperparams_xgb$subsample[j],
      min_child_weight = hyperparams_xgb$min_child_weight[j],
      gamma = hyperparams_xgb$gamma[j],
      colsample_bylevel = hyperparams_xgb$colsample_bylevel[j],
      ROC_AUC_valid = ROC_AUC_valid
    ))
    
    valid_predictions_lf <- predict(last_fold_model, newdata = X_valid)
    ROC_AUC_valid_lf <- pROC::auc(out_of_time_data$occurrenceStatus, valid_predictions_lf)
    
    results_past_lf <- rbind(results_past_lf, data.frame(
      iteration = j,
      nrounds = hyperparams_xgb$nrounds[j],
      max_depth = hyperparams_xgb$max_depth[j],
      eta = hyperparams_xgb$eta[j],
      subsample = hyperparams_xgb$subsample[j],
      min_child_weight = hyperparams_xgb$min_child_weight[j],
      gamma = hyperparams_xgb$gamma[j],
      colsample_bylevel = hyperparams_xgb$colsample_bylevel[j],
      ROC_AUC_valid_lf = ROC_AUC_valid_lf  
    ))
  }
  
  write.csv(results,        paste0('C:/Users/User/Downloads/Downloads/phd_project/results/hpm/', model_name, '_xgb_env_results.csv'), row.names = FALSE)
  write.csv(results_past,   paste0('C:/Users/User/Downloads/Downloads/phd_project/results/hpm/', model_name, '_xgb_env_results_past.csv'), row.names = FALSE)
  write.csv(results_past_lf,paste0('C:/Users/User/Downloads/Downloads/phd_project/results/hpm/', model_name, '_xgb_env_results_past_lf.csv'), row.names = FALSE)
}
