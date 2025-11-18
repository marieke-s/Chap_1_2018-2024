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
      
      train_df <- as.data.frame(st_drop_geometry(train_set))
      test_df  <- as.data.frame(st_drop_geometry(test_set))
      
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
            cat("\n❌ RF failed for fold", i, "of", cv_name, ": ", e$message, "\n")
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
              params = list(
                objective = "reg:squarederror",
                eval_metric = "rmse",
                eta = 0.1,
                max_depth = 6
              ),
              data = train_matrix,
              watchlist = list(train = train_matrix, test = test_matrix),
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
evaluate_models <- function(models_list, response_var = NULL, models = NULL, cv_configs = NULL) {
  
  cat("\n🔹 Starting Model Evaluation\n")
  
  # Initialize performance results storage
  performance_metrics <- data.frame()
  
  # Filter by response variable (if provided)
  response_vars <- if (!is.null(response_var)) {
    intersect(response_var, names(models_list))
  } else {
    names(models_list)  # Use all response variables
  }
  
  for (resp_var in response_vars) {
    
    available_models <- names(models_list[[resp_var]])
    
    model_types <- if (!is.null(models)) {
      intersect(models, available_models)
    } else {
      available_models  # Use all models
    }
    
    for (model_name in model_types) {
      
      cat("\n🔹 Evaluating model:", model_name, "for response variable:", resp_var, "\n")
      
      available_cvs <- names(models_list[[resp_var]][[model_name]])
      
      cv_methods <- if (!is.null(cv_configs)) {
        intersect(cv_configs, available_cvs)
      } else {
        available_cvs  # Use all CV configurations
      }
      
      for (cv_name in cv_methods) {
        
        cat("\n🔹 Processing CV Method:", cv_name, "\n")
        
        model_folds <- models_list[[resp_var]][[model_name]][[cv_name]]
        cv_config <- cv_configs[[cv_name]]  # Get CV splits
        
        all_observed <- c()
        predictions <- c()
        aic_values <- c()
        train_size_list <- c()
        
        for (i in seq_along(model_folds)) {
          
          cat("   ➡️ Fold:", i, "/", length(model_folds), "\n")
          
          fold_data <- model_folds[[paste0("Fold_", i)]]
          
          if (is.null(fold_data)) {
            cat("\n⚠️ Warning: Model missing for fold", i, "- Skipping.\n")
            next
          }
          
          model <- fold_data$model
          test_set <- fold_data$test_set
          
          if (is.null(test_set) || nrow(test_set) != 1) {
            cat("\n⚠️ Warning: Test set missing or incorrect for fold", i, "- Skipping.\n")
            next
          }
          
          test_df <- as.data.frame(st_drop_geometry(test_set))
          
          if (!(resp_var %in% colnames(test_df))) {
            stop(paste("❌ Fold", i, "- Response variable missing in test set!"))
          }
          
          all_observed <- c(all_observed, test_df[[resp_var]])
          
          # --- Prediction handling
          pred <- tryCatch({
            if (model_name == "XGB") {
              test_df <- test_df %>% mutate(across(where(is.character), as.factor))
              
              # Retrieve training feature names
              xgb_features <- model$feature_names
              available_cols <- intersect(colnames(test_df), xgb_features)
              
              test_x_raw <- test_df[, available_cols, drop = FALSE]
              
              # Drop constant columns
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
              
              # If empty, create a zero matrix
              if (ncol(test_x_clean) == 0) {
                test_matrix <- matrix(0, nrow = 1, ncol = length(xgb_features))
                colnames(test_matrix) <- xgb_features
              } else {
                test_matrix <- model.matrix(~ . - 1, data = test_x_clean)
                
                # Add missing columns as 0
                missing_cols <- setdiff(xgb_features, colnames(test_matrix))
                for (col in missing_cols) {
                  test_matrix <- cbind(test_matrix, setNames(data.frame(0), col))
                }
                
                # Reorder
                test_matrix <- test_matrix[, xgb_features, drop = FALSE]
              }
              
              test_matrix <- xgb.DMatrix(data = test_matrix)
              predict(model, newdata = test_matrix)
              
            } else if (model_name == "RF") {
              predict(model, newdata = test_df)
              
            } else {
              predict(model, newdata = test_df, type = "response")
            }
          }, error = function(e) {
            cat("\n❌ Prediction failed for", model_name, "on fold", i, ": ", e$message, "\n")
            return(NA)
          })
          
          predictions <- c(predictions, pred)
          
          if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM")) {
            aic_values <- c(aic_values, AIC(model))
          }
          
          train_size_list <- c(train_size_list, nrow(cv_config[[i]]$train))
        }
        
        compute_metrics <- function(observed, predicted, aic_values) {
          data.frame(
            R2 = cor(observed, predicted, use = "complete.obs")^2,
            Pearson_Corr = cor(observed, predicted, method = "pearson", use = "complete.obs"),
            Spearman_Corr = cor(observed, predicted, method = "spearman", use = "complete.obs"),
            MAE = mean(abs(observed - predicted), na.rm = TRUE),
            RMSE = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
            AIC = if (!is.null(aic_values)) mean(aic_values, na.rm = TRUE) else NA
          )
        }
        
        model_metrics <- compute_metrics(all_observed, predictions, aic_values)
        model_metrics$Response_Var <- resp_var
        model_metrics$Model <- model_name
        model_metrics$CV <- cv_name
        model_metrics$Train_Size <- if (length(train_size_list) > 0) mean(train_size_list) else NA
        
        performance_metrics <- rbind(performance_metrics, model_metrics)
      }
    }
  }
  
  return(performance_metrics)
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

