#------------- Description ---------------------
# Purpose: 

# Author: Marieke Schultz

# Date script created: 18/11/2025

#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(AER)
library(ape)
library(automap)
library(biomod2)
library(blockCV)
library(dplyr)
library(fastDummies)
library(fishtree)
library(geiger)
library(geosphere)  # For distance calculations
library(ggplot2)
library(ggpubr)
library(glmnet)
library(grid)
library(gstat)
library(gridExtra)
library(igraph)
library(MASS)
library(moments)    # For D’Agostino test
library(nortest)   # For Anderson-Darling test
library(pMEM)
library(picante)
library(Rarity)
library(readr)
library(reshape2)
library(rlang)
library(rnaturalearth)  # For land map background
library(sampbias)
library(sf)
library(spdep)
library(superml)
library(svglite)
library(tidyr)
library(units)
library(virtualspecies)
library(xgboost)
library(geosphere)
library(tibble)
library(purrr)




# Load functions
source("./utils/Fct_Modeling.R")

#------------- Load and prep data ------------------
# predictors_sel_v.1.3 ----
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.3.gpkg")

# remove space and majuscule in habitat names
pred$grouped_main_habitat <- gsub(" ", "_", pred$grouped_main_habitat)
pred$grouped_main_habitat <- tolower(pred$grouped_main_habitat)

# set character as factor
pred$grouped_main_habitat <- as.factor(pred$grouped_main_habitat)

# check levels
levels(pred$grouped_main_habitat)
table(pred$grouped_main_habitat)

unique(pred$grouped_main_habitat)

# mtdt_7_sel_v1.1 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

# div_indices_sel_v1.1.gpkg ----
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# tot ----
tot <- pred %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(mtdt), by = "replicates") %>%
  left_join(div, by = "replicates")

#-------------------------- T0.0 : XGBOOST :  ---------------------
## Description ----
#--- Model : XGBOOST

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R, Crypto and Elasmo

#---- Predictors : 
# from predictors_sel_v.1.3
# [1] "northness"                  "eastness"                   "tpi_mean_log"               "port_dist_m_weight"        
# [5] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [9] "temp_mean_1m"               "sal_mean_1m"                "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log" 

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2              Pearson_Corr Spearman_Corr      MAE      RMSE AIC Response_Var Model     CV Train_Size
# 1 0.59889763    0.7738848     0.7738419 8.598846 11.212882  NA            R   XGB bloo50         NA
# 2 0.56089502    0.7489293     0.8188332 2.619474  3.964392  NA       Crypto   XGB bloo50         NA
# 3 0.07769429    0.7489293      0.3864764 1.014452  1.470286  NA       Elasmo   XGB bloo50         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select(-c("replicates"))

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
buffer_sizes <- seq(from = 10, to = 200, by = 10)

# Initialize list to store BLOO-CV results
bloo_results_list <- list()

# Loop through each buffer size and compute BLOO-CV
for (buffer_size in buffer_sizes) {
  
  # Compute BLOO-CV
  bloo_results <- bloo_cv(st_as_sf(tot_sf_4326), buffer_size)
  
  # Store results in list with a named index
  bloo_results_list[[paste0("bloo", buffer_size)]] <- bloo_results
}


# ALL CV
# Combine all CV configurations
cv_configs <- c(list(loo = loo), bloo_results_list)



rm(tot_sf_4326, loo, bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Check and map CV config ----
# Extract BLOO with 50 km buffer
bloo50 <- cv_configs[["bloo50"]]

# Fold sizes: number of points in test (should be 1) and train
fold_sizes_bloo50 <- purrr::map_dfr(
  seq_along(bloo50),
  ~ tibble(
    Fold    = .x,
    n_test  = nrow(bloo50[[.x]]$test),
    n_train = nrow(bloo50[[.x]]$train)
  )
)

print(head(fold_sizes_bloo50))

# Fold n_test n_train
# <int>  <int>   <int>
#   1     1      1     578
# 2     2      1     550
# 3     3      1     550
# 4     4      1     550
# 5     5      1     578
# 6     6      1     555


summary(fold_sizes_bloo50)

# Fold         n_test     n_train     
# Min.   :  1   Min.   :1   Min.   :497.0  
# 1st Qu.:160   1st Qu.:1   1st Qu.:539.0  
# Median :319   Median :1   Median :561.0  
# Mean   :319   Mean   :1   Mean   :558.6  
# 3rd Qu.:478   3rd Qu.:1   3rd Qu.:578.0  
# Max.   :637   Max.   :1   Max.   :599.0  



# --- Map folds ---
extract_fold_df <- function(cv_object, fold_id, keep_geometry = FALSE) {
  
  if (fold_id < 1 || fold_id > length(cv_object)) {
    stop("fold_id is out of range for this CV object.")
  }
  
  train_sf <- cv_object[[fold_id]]$train
  test_sf  <- cv_object[[fold_id]]$test
  
  if (!keep_geometry) {
    train_df <- train_sf |> sf::st_drop_geometry()
    test_df  <- test_sf  |> sf::st_drop_geometry()
  } else {
    train_df <- train_sf
    test_df  <- test_sf
  }
  
  train_df$set <- "train"
  test_df$set  <- "test"
  
  # merge into one data.frame
  df <- dplyr::bind_rows(test_df, train_df)
  
  return(df)
}

for (X in 1:10) {
  f <- extract_fold_df(bloo50, fold_id = X) %>% 
    st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)
  p <- ggplot() +
    geom_sf(data = f, aes(color = set), size = 1) +
    scale_color_manual(values = c("train" = "blue", "test" = "red")) +
    labs(title = "BLOO CV Fold 3: Train (blue) vs Test (red)") +
    theme_minimal()
  print(p)
}



# --- Check distances ---
set.seed(123)
x <- 50  # number of folds to inspect
fold_ids_to_check <- sample(seq_along(bloo50), size = x)

dist_stats <- purrr::map_dfr(fold_ids_to_check, function(i) {
  fold  <- bloo50[[i]]
  test  <- st_coordinates(fold$test)   # 1 x 2 (lon, lat)
  train <- st_coordinates(fold$train)  # n_train x 2
  
  # Distances in metres (Vincenty, same as in bloo_cv)
  d_m  <- geosphere::distVincentySphere(train, test)
  d_km <- d_m / 1000
  
  tibble(
    Fold         = i,
    min_dist_km  = min(d_km, na.rm = TRUE),
    mean_dist_km = mean(d_km, na.rm = TRUE),
    n_train      = nrow(fold$train)
  )
})

head(dist_stats)
summary(dist_stats$min_dist_km)
summary(dist_stats$mean_dist_km)

# Hard check: any training points inside 50 km?
any(dist_stats$min_dist_km < 50)






## Fit model ----------
# Initialize list to store models
models_list <- list()

# Iterate over each response variable
for (response_var in colnames(ind)) {
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["bloo50"],
    models = c("XGB"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}
beepr::beep()

rm(model)

## Save model ----
saveRDS(models_list, "./output/models/xgboost/T0.0.rds") # Response var = R, Crypto and Elasmo

##  Evaluate performance ----
perf <- evaluate_models(models_list)


# Compare with baseline
md_baseline <- readRDS("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1/output/models/MEM_tests/R_bloo30_XGB_T8_pred_4.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline # 0.6613284
perf




#-------------------------- Plot variable importance -------------------
md <- models_list
# perf <- evaluate_models(md)


# 1.  Extract and Aggregate Variable Importance
# Get folds
folds <- md$Crypto$XGB$bloo50

# Extract importance from each fold's model
importance_list <- lapply(folds, function(fold) {
  xgb.importance(model = fold$model)  # <-- change 'model' to actual element name
})

# Combine all importances into one table
all_importance <- bind_rows(importance_list, .id = "Fold")

# Average importance per variable
mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(mean_gain = mean(Gain, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_gain))
















# Add a flag column for highlight -----


# Variables to highlight 
X <- c(
  "TPI_min", "TPI_max", "TPI_mean",
  "TRI_min", "TRI_max", "TRI_mean",
  "slope_min", "slope_max", "slope_mean",
  "aspect_min", "aspect_max", "aspect_mean",
  "roughness_min", "roughness_max", "roughness_mean"
)

X <-c(
  "Infralitoral_Algae_mean", "Other_Seagrass_mean", "Posidonia_mean",
  "Sand_mean", "Dead_Matte_mean", "Rocks_mean", "Coralligeneous_mean"
)

X <- "pMEM"

X <- "bathy"

X <- "gravity"

X <- 'dist_shore'

starts_with_any <- function(feature_name, prefix_vec) {
  any(startsWith(feature_name, prefix_vec))
}

mean_importance <- mean_importance |>
  mutate(
    highlight = ifelse(
      sapply(Feature, starts_with_any, prefix_vec = X),
      "X Variables", "Other"
    )
  )


# Make the plot -----
ggplot(mean_importance, aes(x = reorder(Feature, mean_gain), y = mean_gain)) +
  geom_col(aes(fill = highlight)) +
  geom_text(
    aes(label = round(mean_gain, 3), color = highlight),
    hjust = -0.1,
    size = 3
  ) +
  scale_fill_manual(values = c("X Variables" = "tomato", "Other" = "grey70")) +
  scale_color_manual(values = c("X Variables" = "tomato", "Other" = "grey30")) +
  coord_flip() +
  labs(
    title = "Mean Variable Importance Across Folds",
    x = "Feature",
    y = "Mean Gain",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ylim(0, max(mean_importance$mean_gain) * 1.15)  # Add space for text























#-------------------------- T0.1 : XGBOOST :  ENV-k6-CV---------------------
## Description ----
#--- Model : XGBOOST

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R                                                                   <== CHANGE

#---- Predictors : 
# from predictors_sel_v.1.3
# [1] "northness"                  "eastness"                   "tpi_mean_log"               "port_dist_m_weight"        
# [5] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [9] "temp_mean_1m"               "sal_mean_1m"                "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log" 

#---- CV config :
# ENV CV with k = 6 (see 8_Cross_validation.R for k-means clustering) <== CHANGE

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.1699394     0.412237     0.3962696 12.88808 15.89576  NA            R   XGB env_k6   530.8333



## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)
if (st_crs(tot_sf_4326)$epsg != 4326) {            
  stop("Error: tot_sf_4326 is not in EPSG:4326 coordinate reference system.")
}
stopifnot(nrow(predictors) == nrow(tot_sf_4326))




# Load CV cluster assignments csv
cv <- readr::read_csv("./data/processed_data/Cross-val/cluster_assignments_k6_envCV_T0.0.csv")
stopifnot(nrow(cv) == nrow(tot_sf_4326))
if (!all(cv$replicates %in% tot_sf_4326$replicates)) {
  stop("Some 'replicates' in the CV file are not found in 'tot_sf_4326'.")
}




# Attach cluster info to tot_sf_4326 by 'replicates'
tot_sf_4326 <- tot_sf_4326 %>%
  dplyr::left_join(cv, by = "replicates")
if (any(is.na(tot_sf_4326$cluster))) {
  stop("NA values in 'cluster' after join – check replicates matching.")
}


# Cluster assignments (k = 6 for env CV)
cluster_assignments <- tot_sf_4326$cluster
unique_clusters     <- sort(unique(cluster_assignments))

# Build environmental CV folds in the same format as loo/bloo
env_cv <- vector("list", length = length(unique_clusters))

for (i in seq_along(unique_clusters)) {
  k            <- unique_clusters[i]
  test_idx     <- which(cluster_assignments == k)
  train_idx    <- setdiff(seq_len(nrow(tot_sf_4326)), test_idx)
  
  env_cv[[i]] <- list(
    train = tot_sf_4326[train_idx, ],
    test  = tot_sf_4326[test_idx, ]
  )
}

names(env_cv) <- paste0("fold_", seq_along(env_cv))

# Optional: print fold sizes
total_points <- nrow(tot_sf_4326)
for (i in seq_along(env_cv)) {
  fold_points      <- nrow(env_cv[[i]]$test)
  remaining_points <- total_points - fold_points
  cat("Fold ", i, ": test : ", fold_points, ", train : ", remaining_points, "\n", sep = "")
}

# Create cv_configs object using only environmental CV (k=6)
cv_configs <- list(env_k6 = env_cv)

# Clean up
rm(tot_sf_4326, cv, cluster_assignments, unique_clusters, env_cv)

# Fold 1: test : 79, train : 558
# Fold 2: test : 137, train : 500
# Fold 3: test : 39, train : 598
# Fold 4: test : 74, train : 563
# Fold 5: test : 123, train : 514
# Fold 6: test : 185, train : 452









## Fit model ----------
# Initialize list to store models
models_list <- list()

# Iterate over each response variable
for (response_var in colnames(ind)) {   # Just R
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["env_k6"],
    models = c("XGB"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}
beepr::beep()



## Save model ----
saveRDS(models_list, "./output/models/xgboost/T0.1.rds") # Response var = R


##  Evaluate performance ----
perf <- evaluate_models(models_list)


# Compare with baseline
md_baseline <- readRDS("./output/models/xgboost/T0.0.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline # 0.7738848
perf # 0.412237 

perf









# models_list as you already have
# cv_configs is your list with env_k6, bloo50, etc.

# 1) Summary performance (average over folds) – same spirit as before
perf_summary <- evaluate_models(
  models_list,
  cv_splits = cv_configs,      # to compute mean train size; optional
  output    = "summary"
)

# 2) Per-fold performance table
perf_folds <- evaluate_models(
  models_list,
  cv_splits = cv_configs,
  output    = "fold"
)

# 3) Plots (observed vs predicted)
plots <- evaluate_models(
  models_list,
  cv_splits = cv_configs,
  output    = "plot"
)

# e.g. print a given plot in your R session:
print(plots[["R_XGB_env_k6"]])   # adjust name depending on your resp/model/cv


























#-------------------------- T0.2 : XGBOOST : RANDOM-CV ---------------------
## Description ----
#--- Model : XGBOOST

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R                                                                   


#---- Predictors : 
# from predictors_sel_v.1.3
# [1] "northness"                  "eastness"                   "tpi_mean_log"               "port_dist_m_weight"        
# [5] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [9] "temp_mean_1m"               "sal_mean_1m"                "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log" 


#---- CV config :
# Random CV with 20% test / 80% train, repeated n_folds times         <== CHANGE

#---- Perf :
# R2 Pearson_Corr Spearman_Corr     MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4142866     0.643651     0.6389859 10.0907 13.32401  NA            R   XGB rand20         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>%
  dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variable
ind <- div %>% dplyr::select("R")


## CV config -----
# Convert tot to sf object 
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE)
if (st_crs(tot_sf_4326)$epsg != 4326) {            
  stop("Error: tot_sf_4326 is not in EPSG:4326 coordinate reference system.")
}
stopifnot(nrow(predictors) == nrow(tot_sf_4326))





# --- Random CV: 20% test / 80% train ---
set.seed(123)  # for reproducibility

n        <- nrow(tot_sf_4326)
test_prop <- 0.2
n_test    <- floor(test_prop * n)
n_folds   <- 5  # number of random splits / repeats (change if you wish)

rand_cv <- vector("list", length = n_folds)

for (i in seq_len(n_folds)) {
  test_idx  <- sample(seq_len(n), size = n_test, replace = FALSE)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  rand_cv[[i]] <- list(
    train = tot_sf_4326[train_idx, ],
    test  = tot_sf_4326[test_idx, ]
  )
}

names(rand_cv) <- paste0("fold_", seq_len(n_folds))

# Optional: check fold sizes
total_points <- nrow(tot_sf_4326)
for (i in seq_along(rand_cv)) {
  fold_points      <- nrow(rand_cv[[i]]$test)
  remaining_points <- total_points - fold_points
  cat("Random CV - Fold ", i, ": test : ", fold_points,
      ", train : ", remaining_points, "\n", sep = "")
}

# Random CV - Fold 1: test : 127, train : 510
# Random CV - Fold 2: test : 127, train : 510
# Random CV - Fold 3: test : 127, train : 510
# Random CV - Fold 4: test : 127, train : 510
# Random CV - Fold 5: test : 127, train : 510


# CV configs object for this run (only random CV)
cv_configs <- list(rand20 = rand_cv)



## Check and map cv config ----
# Map of folds (points colored by fold_id)
ggplot(tot_sf_4326) +
  geom_sf(aes(color = fold_id), size = 2, alpha = 0.8) +
  labs(
    title = "Random 5-fold CV (20% test each fold)",
    color = "Fold"
  ) +
  theme_minimal()


# Clean up
rm(tot_sf_4326, n, n_test, test_prop, n_folds, rand_cv, fold_points, remaining_points, i)













## Fit model ----------
# Initialize list to store models
models_list <- list()

# Iterate over each response variable
for (response_var in colnames(ind)) { # Just R
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["rand20"],
    models = c("XGB"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

beepr::beep()


## Save model ----
saveRDS(models_list, "./output/models/xgboost/T0.2.rds") # Response var = R

##  Evaluate performance ----
perf <- evaluate_models(models_list)
perf

# Per fold performance
perf_folds <- evaluate_models(
  models_list,
  cv_splits = cv_configs,
  output    = "fold"
)






















