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
library(lubridate)
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
library(RandomForestsGLS)
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
# dev.off()
#===============================================================================
#-------------------------- T.1 : Load and prep data ------------------
# predictors_sel_v.1.4 ----
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.4.gpkg")

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


#===============================================================================
#-------------------------- T.0 : Load and prep data ------------------
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

#===============================================================================

#-------------------------- T0.0 : XGBOOST : bloo50  ---------------------
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
# 3 0.07769429    0.2787370     0.3864764 1.014452  1.470286  NA       Elasmo   XGB bloo50         NA



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
  
  # Build the fold df + convert back to sf
  f <- extract_fold_df(bloo50, fold_id = X) %>% 
    st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)
  
  # Build plot
  p <- ggplot() +
    geom_sf(data = f, aes(color = set), size = 1) +
    scale_color_manual(values = c("train" = "blue", "test" = "red")) +
    labs(title = paste("BLOO CV", X, ": Train (blue) vs Test (red)")) +
    theme_minimal()
  
  # Print to RStudio plotting window
  print(p)
  
  # Save plot
  filename <- sprintf("./figures/Cross-val/Map_Bloo50_fold%d_T0.0.png", X)
  ggsave(filename, plot = p, width = 7, height = 6, dpi = 300)
}



# Panel plot of 12 randomly selected folds
set.seed(123)  # for reproducibility

# 1) Randomly select 12 distinct folds from BLOO CV
n_folds_total <- length(bloo50)
n_sample      <- 12

fold_ids_sample <- sample(seq_len(n_folds_total), size = n_sample)

# 2) Build combined sf object for these folds
panel_data <- map_dfr(fold_ids_sample, function(fid) {
  df <- extract_fold_df(bloo50, fold_id = fid)
  df$panel_fold <- paste0("Fold ", fid)  # label for facet
  df
})

panel_data_sf <- panel_data %>%
  st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)

# 3) Build 4 x 3 panel plot
p_panel <- ggplot(panel_data_sf) +
  geom_sf(aes(color = set), size = 0.8) +
  scale_color_manual(values = c(train = "blue", test = "red")) +
  facet_wrap(~ panel_fold, ncol = 4) +
  labs(
    title    = "BLOO CV (50 km) – Sample of 12 Folds",
    subtitle = "Train (blue) vs Test (red)",
    color    = "Set"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

# 4) Print to device
print(p_panel)

# 5) Save single figure
outfile <- "./figures/Cross-val/Map_Bloo50_12folds_T0.0.png"
ggsave(outfile, plot = p_panel, width = 10, height = 8, dpi = 300)


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
models_list <- readRDS("./output/models/xgboost/T0.0.rds")
perf <- evaluate_models(models_list)
perf

# Save perf
saveRDS(perf, "./output/models/xgboost/T0.0_performance.rds")

# Compare with baseline
md_baseline <- readRDS("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1/output/models/MEM_tests/R_bloo30_XGB_T8_pred_4.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline # 0.6613284
perf




#-------------------------- T0.0.1 : XGBOOST : bloo50  [on marbec-gpu] ---------------------
## Description ----
#--- Model : XGBOOST

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   "AngelShark"  <== CHANGE

#---- Predictors : 
# from predictors_sel_v.1.3
# [1] "northness"                  "eastness"                   "tpi_mean_log"               "port_dist_m_weight"        
# [5] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [9] "temp_mean_1m"               "sal_mean_1m"                "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log" 

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr        MAE       RMSE AIC Response_Var Model     CV Train_Size
# 1 0.1234221    0.3513149     0.6578672 0.17644896 0.35649337  NA        DeBRa   XGB bloo50         NA
# 2 0.4868173    0.6977230     0.7671870 0.88346354 1.31262646  NA      RedList   XGB bloo50         NA
# 3 0.5956218    0.7717654     0.7126899 6.19833542 8.06845462  NA         LRFI   XGB bloo50         NA
# 4 0.5864132    0.7657762     0.7642129 1.66474540 2.36872499  NA      TopPred   XGB bloo50         NA
# 5 0.6279663    0.7924432     0.7570144 4.27136071 5.54144664  NA   Commercial   XGB bloo50         NA
# 6 0.3171457    0.5631569     0.5327619 0.30929279 0.37687384  NA      Grouper   XGB bloo50         NA
# 7 0.2952172    0.5433390     0.1751339 0.01780742 0.09190353  NA   AngelShark   XGB bloo50         NA

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























# [[MODEL FITTED ON MARBEC-GPU]] ----
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T0.0.1.rds")
perf <- evaluate_models(models_list)
perf

# Save perf 
saveRDS(perf, "./output/models/xgboost/T0.0.1_performance.rds")




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


























#-------------------------- T0.2 : XGBOOST : RANDOM-CV 80-20 ---------------------
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

# --- Random K-fold CV: 5 folds (~20% each) ---
set.seed(123)  # for reproducibility

n      <- nrow(tot_sf_4326)
K      <- 5   # number of folds

# Random fold assignment: each point gets exactly one fold_id
fold_id <- sample(rep(1:K, length.out = n))

# Attach fold_id to sf (for mapping etc.)
tot_sf_4326$fold_id <- factor(fold_id)

# Build CV object: for fold k, test = fold k, train = all other folds
rand_cv <- vector("list", length = K)

for (k in seq_len(K)) {
  test_idx  <- which(fold_id == k)
  train_idx <- which(fold_id != k)
  
  rand_cv[[k]] <- list(
    train = tot_sf_4326[train_idx, ],
    test  = tot_sf_4326[test_idx, ]
  )
}

names(rand_cv) <- paste0("fold_", seq_len(K))

# Optional: check fold sizes
total_points <- nrow(tot_sf_4326)
for (k in seq_along(rand_cv)) {
  fold_points      <- nrow(rand_cv[[k]]$test)
  remaining_points <- total_points - fold_points
  cat("Random K-fold - Fold ", k, ": test : ", fold_points,
      ", train : ", remaining_points, "\n", sep = "")
}


# Random K-fold - Fold 1: test : 128, train : 509
# Random K-fold - Fold 2: test : 128, train : 509
# Random K-fold - Fold 3: test : 127, train : 510
# Random K-fold - Fold 4: test : 127, train : 510
# Random K-fold - Fold 5: test : 127, train : 510



# CV configs object for this run (only random CV)
cv_configs <- list(rand20 = rand_cv)

## Check and map cv config ----
# ensure fold_id is a factor
tot_sf_4326$fold_id <- factor(tot_sf_4326$fold_id)

fold_levels <- levels(tot_sf_4326$fold_id)

ggplot(panel_data) +
  geom_sf(aes(color = set), size = 1) +
  facet_wrap(~ panel_fold, ncol = 3) +
  scale_color_manual(
    values = c(train = "grey80", test = "red")
  ) +
  labs(
    title = "Random 5-fold CV",
    subtitle = "Each panel: test fold in red, others as training (grey)",
    color = "Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )










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
perf # 0.6446923

# Per fold performance
perf_folds <- evaluate_models(
  models_list,
  cv_splits = cv_configs,
  output    = "fold"
)
































#-------------------------- T0.3 : XGBOOST : SPATIAL-CV-50, 100, 200km ---------------------
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


#--- CV config :                                                      <== CHANGE
# spatial CV using blockCV::cv_spatial with spatial block sizes = 50 km, 100 km, 200 km


#--- Perf :
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model    CV Train_Size
# 1 0.2856066    0.5344218     0.5291810 11.65594 14.61688  NA            R   XGB  sp50         NA
# 2 0.2068828    0.4548438     0.4501214 12.37006 15.33875  NA            R   XGB sp100         NA
# 3 0.1839316    0.4288725     0.4002176 12.57053 15.69570  NA            R   XGB sp200         NA





## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select(("R"))




## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)
if (st_crs(tot_sf_4326)$epsg != 4326) {            
  stop("Error: tot_sf_4326 is not in EPSG:4326 coordinate reference system.")
}
stopifnot(nrow(predictors) == nrow(tot_sf_4326))


# Reproject to a metric CRS for correct km units (e.g. EPSG:3035)
tot_sf_3035 <- st_transform(tot_sf_4326, 3035)

library(blockCV)
## Helper: build cv_config from cv_spatial() ---
build_spatial_cv <- function(sf_3035,
                             sf_4326,
                             block_km,
                             k_folds = 5,
                             n_iter  = 100,
                             seed    = 123) {
  
  message("\n--- Building spatial CV with block size = ", block_km, " km ---")
  
  block_size_m <- block_km * 1000
  
  spat_cv_obj <- cv_spatial(
    x         = sf_3035,
    column    = NULL,           # no response column used to balance classes
    size      = block_size_m,   # metres, since CRS is metric
    k         = k_folds,
    hexagon   = TRUE,
    selection = "random",
    iteration = n_iter,
    biomod2   = FALSE,
    plot      = FALSE,          # <- avoid cv_plot / species argument issue
    progress  = TRUE,
    report    = TRUE,
    seed      = seed
  )
  
  # Summary table (same as what is printed)
  print(spat_cv_obj$records)
  
  # Convert folds_list (indices) into your standard list(list(train=sf, test=sf))
  folds_list <- spat_cv_obj$folds_list   # <--- THIS is the correct element
  
  this_cv <- lapply(seq_along(folds_list), function(i) {
    idx_train <- folds_list[[i]][[1]]
    idx_test  <- folds_list[[i]][[2]]
    
    list(
      train = sf_4326[idx_train, ],
      test  = sf_4326[idx_test, ]
    )
  })
  
  names(this_cv) <- paste0("fold_", seq_along(this_cv))
  
  return(this_cv)
}

## Build the three spatial CV configs ---
sp50  <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 50)
sp100 <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 100)
sp200 <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 200)

# Fold sizes for a quick sanity check
fold_sizes_sp50 <- purrr::map_dfr(
  seq_along(sp50),
  ~ tibble(
    Fold    = .x,
    n_test  = nrow(sp50[[.x]]$test),
    n_train = nrow(sp50[[.x]]$train)
  )
)
print(fold_sizes_sp50)

#     Fold n_test n_train
# 1     1    117     520
# 2     2    125     512
# 3     3    125     512
# 4     4    117     520
# 5     5    153     484





fold_sizes_sp100 <- purrr::map_dfr(
  seq_along(sp100),
  ~ tibble(
    Fold    = .x,
    n_test  = nrow(sp100[[.x]]$test),
    n_train = nrow(sp100[[.x]]$train)
  )
)

print(fold_sizes_sp100)


# Fold n_test n_train
# 1     1    133     504
# 2     2    121     516
# 3     3    146     491
# 4     4    108     529
# 5     5    129     508





fold_sizes_sp200 <- purrr::map_dfr(
  seq_along(sp200),
  ~ tibble(
    Fold    = .x,
    n_test  = nrow(sp200[[.x]]$test),
    n_train = nrow(sp200[[.x]]$train)
  )
)
print(fold_sizes_sp200)

# Fold n_test n_train
# 1     1     65     572
# 2     2    214     423
# 3     3    187     450
# 4     4     33     604
# 5     5    138     499





## Combine into cv_configs for T0.3 ---
cv_configs <- list(
  sp50  = sp50,
  sp100 = sp100,
  sp200 = sp200
)

## Map  CV (e.g. 50 km, fold 1) ----
# Combine test/train into one sf for a given fold
extract_spatial_fold_df <- function(cv_object, fold_id) {
  train_sf <- cv_object[[fold_id]]$train
  test_sf  <- cv_object[[fold_id]]$test
  
  train_sf$set <- "train"
  test_sf$set  <- "test"
  
  dplyr::bind_rows(train_sf, test_sf)
}

fold_id_example <- 2
fold_sf <- extract_spatial_fold_df(sp50, fold_id_example)

ggplot(fold_sf) +
  geom_sf(aes(color = set), size = 1) +
  scale_color_manual(values = c(train = "blue", test = "red")) +
  labs(
    title = paste0("Spatial CV (50 km) – Fold ", fold_id_example),
    subtitle = "Train (blue) vs Test (red)"
  ) +
  theme_minimal()


# Panel plot
make_spatial_panel <- function(cv_obj, outfile, title_prefix) {
  panel_sf <- purrr::map_dfr(seq_along(cv_obj), function(fid) {
    sf_fold <- extract_spatial_fold_df(cv_obj, fold_id = fid)
    sf_fold$fold_lab <- paste0("Fold ", fid)
    sf_fold
  })
  
  p <- ggplot(panel_sf) +
    geom_sf(aes(color = set), size = 0.8) +
    scale_color_manual(values = c(train = "blue", test = "red")) +
    facet_wrap(~ fold_lab, ncol = 3) +
    labs(
      title    = paste0(title_prefix, " – 5-fold"),
      subtitle = "Train (blue) vs Test (red)",
      color    = "Set"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  ggsave(outfile, plot = p, width = 10, height = 8, dpi = 300)
}

make_spatial_panel(sp50,  "./figures/Cross-val/Map_spatialCV_50km_5folds_T0.3.png", "Spatial CV (50 km)")
make_spatial_panel(sp100, "./figures/Cross-val/Map_spatialCV_100km_5folds_T0.3.png", "Spatial CV (100 km)")
make_spatial_panel(sp200, "./figures/Cross-val/Map_spatialCV_200km_5folds_T0.3.png", "Spatial CV (200 km)")




















## Plots folds -----
build_spatial_cv <- function(sf_3035,
                             sf_4326,
                             block_km,
                             k_folds = 5,
                             n_iter  = 100,
                             seed    = 123) {
  
  message("\n--- Building spatial CV with block size = ", block_km, " km ---")
  
  block_size_m <- block_km * 1000
  
  spat_cv_obj <- cv_spatial(
    x         = sf_3035,
    column    = NULL,
    size      = block_size_m,
    k         = k_folds,
    hexagon   = TRUE,
    selection = "random",
    iteration = n_iter,
    biomod2   = FALSE,
    plot      = FALSE,
    progress  = TRUE,
    report    = TRUE,
    seed      = seed
  )
  
  folds_list <- spat_cv_obj$folds_list
  
  this_cv <- lapply(seq_along(folds_list), function(i) {
    idx_train <- folds_list[[i]][[1]]
    idx_test  <- folds_list[[i]][[2]]
    list(
      train = sf_4326[idx_train, ],
      test  = sf_4326[idx_test, ]
    )
  })
  
  names(this_cv) <- paste0("fold_", seq_along(this_cv))
  
  # return both items
  return(list(
    folds = this_cv,
    obj   = spat_cv_obj
  ))
}



sp50_full  <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 50)
sp100_full <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 100)
sp200_full <- build_spatial_cv(tot_sf_3035, tot_sf_4326, block_km = 200)

sp50  <- sp50_full$folds
sp100 <- sp100_full$folds
sp200 <- sp200_full$folds

sp50_obj  <- sp50_full$obj
sp100_obj <- sp100_full$obj
sp200_obj <- sp200_full$obj


blocks50_sf <- st_as_sf(sp50_obj$blocks)

ggplot() +
  geom_sf(data = blocks50_sf, fill = NA, color = "black") +
  geom_sf(data = tot_sf_4326, aes(color = "Points"), size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Points" = "blue")) +
  labs(
    title = "Spatial Blocks (50 km)",
    subtitle = "Underlying hexagonal blocks generated by cv_spatial()"
  ) +
  theme_minimal()


# add fold ID to each point
fold_id_vec <- rep(NA, nrow(tot_sf_4326))

for (k in seq_along(sp50)) {
  idx_test <- rownames(sp50[[k]]$test)
  idx_test <- as.numeric(idx_test)
  fold_id_vec[idx_test] <- k
}

tot_sf_4326$fold_id <- factor(fold_id_vec)

ggplot() +
  geom_sf(data = blocks50_sf, fill = NA, color = "grey50") +
  geom_sf(data = tot_sf_4326, aes(color = fold_id), size = 1) +
  scale_color_viridis_d(na.value = "lightgrey") +
  labs(
    title = "Spatial CV (50 km) – Blocks + Fold Assignment",
    color = "Fold"
  ) +
  theme_minimal()


## Fit model ----------
# Initialize list to store models
models_list_T0.3 <- list()

for (response_var in colnames(ind)) {
  model <- fit_models(
    dataset    = tot,
    response_var = response_var,
    predictors = colnames(predictors),
    cv_configs = cv_configs,       # sp50, sp100, sp200
    models     = c("XGB"),
    distribution = NULL
  )
  models_list_T0.3[[response_var]] <- model
}




## Save model -----
saveRDS(models_list_T0.3, "./output/models/xgboost/T0.3.rds") # Response var = R
models_list <- models_list_T0.3
rm(models_list_T0.3)



##  Evaluate performance ----
perf <- evaluate_models(models_list)
perf 

# 1 0.2856066    0.5344218     0.5291810 11.65594 14.61688  NA            R   XGB  sp50         NA
# 2 0.2068828    0.4548438     0.4501214 12.37006 15.33875  NA            R   XGB sp100         NA
# 3 0.1839316    0.4288725     0.4002176 12.57053 15.69570  NA            R   XGB sp200         NA

# Per fold performance
perf_folds <- evaluate_models(
  models_list,
  cv_splits = cv_configs,
  output    = "fold"
)


































#-------------------------- T0.4 : XGBOOST : RANDOM-CV 90-10---------------------
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
# Random CV with 10% test / 90% train, repeated n_folds times         <== CHANGE

#---- Perf :
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4153888    0.6445066      0.632044 10.36372 13.15706  NA            R   XGB rand20         NA



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

# --- Random K-fold CV: 10 folds (~10% each) ---
set.seed(123)  # reproducible

n <- nrow(tot_sf_4326)
K <- 10        # 10 folds → 10% test each

# Random fold assignment
fold_id <- sample(rep(1:K, length.out = n))

# Attach to sf for mapping
tot_sf_4326$fold_id <- factor(fold_id)

# Build CV object
rand_cv <- vector("list", length = K)

for (k in seq_len(K)) {
  test_idx  <- which(fold_id == k)
  train_idx <- which(fold_id != k)
  
  rand_cv[[k]] <- list(
    train = tot_sf_4326[train_idx, ],
    test  = tot_sf_4326[test_idx, ]
  )
}

names(rand_cv) <- paste0("fold_", seq_len(K))


# Optional: check fold sizes
total_points <- nrow(tot_sf_4326)
for (k in seq_along(rand_cv)) {
  fold_points      <- nrow(rand_cv[[k]]$test)
  remaining_points <- total_points - fold_points
  cat("Random K-fold - Fold ", k, ": test : ", fold_points,
      ", train : ", remaining_points, "\n", sep = "")
}


# Random K-fold - Fold 1: test : 64, train : 573
# Random K-fold - Fold 2: test : 64, train : 573
# Random K-fold - Fold 3: test : 64, train : 573
# Random K-fold - Fold 4: test : 64, train : 573
# Random K-fold - Fold 5: test : 64, train : 573
# Random K-fold - Fold 6: test : 64, train : 573
# Random K-fold - Fold 7: test : 64, train : 573
# Random K-fold - Fold 8: test : 63, train : 574
# Random K-fold - Fold 9: test : 63, train : 574
# Random K-fold - Fold 10: test : 63, train : 574


# CV configs object for this run (only random CV)
cv_configs <- list(rand20 = rand_cv)

## Check and map cv config ----
panel_data <- map_dfr(seq_along(rand_cv), function(k) {
  fold <- rand_cv[[k]]
  
  train_sf <- fold$train
  test_sf  <- fold$test
  
  train_sf$set <- "train"
  test_sf$set  <- "test"
  
  bind_rows(train_sf, test_sf) %>%
    mutate(panel_fold = paste0("Fold ", k))
})

# make sure it’s sf (it already is, but this is safe)
panel_data <- st_as_sf(panel_data)

# ensure fold_id is a factor
tot_sf_4326$fold_id <- factor(tot_sf_4326$fold_id)

fold_levels <- levels(tot_sf_4326$fold_id)

ggplot(panel_data) +
  geom_sf(aes(color = set), size = 1) +
  facet_wrap(~ panel_fold, ncol = 3) +
  scale_color_manual(
    values = c(train = "grey80", test = "red")
  ) +
  labs(
    title = "Random K-fold CV",
    subtitle = "Each panel: test fold in red, others as training (grey)",
    color = "Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave("./figures/Cross-val/Map_RandomCV_10-90_T0.2.png",
       width = 10, height = 8, dpi = 300)







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
saveRDS(models_list, "./output/models/xgboost/T0.4.rds") # Response var = R

##  Evaluate performance ----
perf <- evaluate_models(models_list)
perf # 0.6445066


































#-------------------------- T1.0 : XGBOOST : bloo50  [on marbec-gpu]---------------------
## Description ----
#--- Model : XGBOOST

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R, Crypto,  Elasmo "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   "AngelShark"

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2              Pearson_Corr Spearman_Corr        MAE        RMSE AIC Response_Var Model     CV Train_Size
# 1  0.60698014    0.7790893     0.7695691 8.46826331 11.10390361  NA            R   XGB bloo50         NA
# 2  0.58232803    0.7631042     0.8305654 2.57250228  3.84527522  NA       Crypto   XGB bloo50         NA
# 3  0.10702283    0.3271434     0.4081011 0.99986696  1.45600642  NA       Elasmo   XGB bloo50         NA
# 4  0.04856717    0.2203796     0.5050083 0.20353596  0.36573106  NA        DeBRa   XGB bloo50         NA
# 5  0.51729240    0.7192304     0.7875041 0.88075750  1.29827096  NA      RedList   XGB bloo50         NA
# 6  0.58726090    0.7663295     0.7025712 6.33019028  8.16209304  NA         LRFI   XGB bloo50         NA
# 7  0.57787778    0.7601827     0.7598383 1.68161264  2.38770971  NA      TopPred   XGB bloo50         NA
# 8  0.62866074    0.7928813     0.7665601 4.26526788  5.58965117  NA   Commercial   XGB bloo50         NA
# 9  0.27898353    0.5281889     0.4995802 0.31641806  0.38295882  NA      Grouper   XGB bloo50         NA
# 10 0.29125019    0.5396760     0.1751339 0.01821245  0.09242783  NA   AngelShark   XGB bloo50         NA





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


## Check CV config ----
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








# [[MODEL FITTED ON MARBEC-GPU]] ----
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.0.rds")
perf <- evaluate_models(models_list)
class(perf)

# Save perf 
saveRDS(perf, "./output/models/xgboost/T1.0_performance.rds")
perf <- readRDS("./output/models/xgboost/T1.0_performance.rds")

# Compare with baseline
md_baseline <- readRDS("./output/models/xgboost/T0.0.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline 
perf






#-------------------------- T1.1 : XGBOOST : no early stopping no watchlist [on marbec-gpu]---------------------
## Description ----
#--- Model : XGBOOST
# parameter "early_stopping_rounds" disabled (default = NULL)         <== CHANGE

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R, Crypto,  Elasmo "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   "AngelShark"

#---- Predictors : 
# from predictors_sel_v.1.3 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 




##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.1.rds")
perf <- evaluate_models(models_list)
beepr::beep()

# Save perf 
saveRDS(perf, "./output/models/xgboost/T1.1_performance.rds")

# Compare with baseline
md_baseline <- readRDS("./output/models/xgboost/T1.0.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline 
perf


md_baseline$R$XGB$bloo50$Fold_1$model



²# > perf_baseline 
# R2 Pearson_Corr Spearman_Corr        MAE        RMSE AIC Response_Var Model     CV Train_Size
# 1  0.60698014    0.7790893     0.7695691 8.46826331 11.10390361  NA            R   XGB bloo50         NA
# 2  0.58232803    0.7631042     0.8305654 2.57250228  3.84527522  NA       Crypto   XGB bloo50         NA
# 3  0.10702283    0.3271434     0.4081011 0.99986696  1.45600642  NA       Elasmo   XGB bloo50         NA
# 4  0.04856717    0.2203796     0.5050083 0.20353596  0.36573106  NA        DeBRa   XGB bloo50         NA
# 5  0.51729240    0.7192304     0.7875041 0.88075750  1.29827096  NA      RedList   XGB bloo50         NA
# 6  0.58726090    0.7663295     0.7025712 6.33019028  8.16209304  NA         LRFI   XGB bloo50         NA
# 7  0.57787778    0.7601827     0.7598383 1.68161264  2.38770971  NA      TopPred   XGB bloo50         NA
# 8  0.62866074    0.7928813     0.7665601 4.26526788  5.58965117  NA   Commercial   XGB bloo50         NA
# 9  0.27898353    0.5281889     0.4995802 0.31641806  0.38295882  NA      Grouper   XGB bloo50         NA
# 10 0.29125019    0.5396760     0.1751339 0.01821245  0.09242783  NA   AngelShark   XGB bloo50         NA
# > perf
# R2 Pearson_Corr Spearman_Corr         MAE       RMSE AIC Response_Var Model     CV Train_Size
# 1  3.747176e-03  0.061214185   0.055522392 15.56891857 18.7699076  NA            R   XGB bloo50         NA
# 2  2.402349e-05 -0.004901376   0.042826156  4.75582933  5.5277563  NA       Crypto   XGB bloo50         NA
# 3  1.778358e-03  0.042170585   0.039250498  1.10332687  1.5376973  NA       Elasmo   XGB bloo50         NA
# 4  1.147989e-02 -0.107144239  -0.049152797  0.21548838  0.3773059  NA        DeBRa   XGB bloo50         NA
# 5  2.942464e-03 -0.054244487  -0.033827228  1.32372439  1.6531809  NA      RedList   XGB bloo50         NA
# 6  7.312943e-03 -0.085515748  -0.074975227 11.28926087 13.9305861  NA         LRFI   XGB bloo50         NA
# 7  3.918910e-04  0.019796237   0.049026022  2.96086254  3.5461398  NA      TopPred   XGB bloo50         NA
# 8  1.180931e-02  0.108670630   0.091861660  7.66376832  9.3817172  NA   Commercial   XGB bloo50         NA
# 9  2.953504e-04 -0.017185761   0.006184936  0.31614660  0.4716645  NA      Grouper   XGB bloo50         NA
# 10 1.311917e-02  0.114538960   0.123745494  0.01099401  0.1048279  NA   AngelShark   XGB bloo50         NA



#-------------------------- T1.2 : XGBOOST : no early stopping (watchlist = list(train = train_matrix)) [on marbec-gpu]---------------------
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.2.rds")
perf <- evaluate_models(models_list)
beepr::beep()
models_list$R$XGB$bloo50$Fold_1$model

# Save perf 
saveRDS(perf, "./output/models/xgboost/T1.2_performance.rds")

# Compare with baseline
md_baseline <- readRDS("./output/models/xgboost/T1.0.rds")
perf_baseline <- evaluate_models(md_baseline)
perf_baseline 
perf


md_baseline$R$XGB$bloo50$Fold_1$model














#-------------------------- T1.3 : ------------------------------------------
# xgb_model <- tryCatch(
#   xgb.train(
#     params = list(
#       colsample_bytree = 0.8,
#       objective = "reg:squarederror",
#       eval_metric = "rmse",
#       eta = 0.1,
#       max_depth = 5
#     ),
#     data = train_matrix,
#     watchlist = list(train = train_matrix),
#     nrounds = 100,
#     early_stopping_rounds = 10,
#     verbose = 0
#   )

##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.3.rds")
perf <- evaluate_models(models_list)
models_list$R$XGB$bloo50$Fold_1$model
# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.009963534    0.0998175     0.1042987 15.38325 18.5497  NA            R   XGB bloo50         NA





#-------------------------- T1.4 ----------------------------------
# xgb_model <- tryCatch(
#   xgb.train(
#     params = list(
#       subsample = 0.9
#       colsample_bytree = 0.8,
#       objective = "reg:squarederror",
#       eval_metric = "rmse",
#       eta = 0.1,
#       max_depth = 6
#     ),
#     data = train_matrix,
#     watchlist = list(train = train_matrix),
#     nrounds = 200,
#     early_stopping_rounds = 10,
#     verbose = 0



##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.4.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model




#-------------------------- T1.5 ----------------------------------
# xgb_model <- tryCatch(
#   xgb.cv(
#     nfold = 10,
#     params = list(
#       objective = "reg:squarederror",
#       eval_metric = "rmse",
#       eta = 0.1,
#       max_depth = 6
#     ),
#     data = train_matrix,
#     nrounds = 100,
#     early_stopping_rounds = 10,
#     verbose = 0
#   ),
#   error = function(e) {
#     cat("\n XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
#     return(NULL)
#   }
# )


##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.5.rds")
models_list$R$XGB$bloo50$Fold_1$model




















#-------------------------- T1.6 ----------------------------------
# # Train and predict
# if (!is.null(train_matrix) && !is.null(test_matrix)) {
#   xgb_cv <- tryCatch(
#     xgb.cv(
#       nfold = 10,
#       params = list(
#         objective = "reg:squarederror",
#         eval_metric = "rmse",
#         eta = 0.1,
#         max_depth = 6
#       ),
#       data = train_matrix,
#       nrounds = 100,
#       early_stopping_rounds = 10,
#       verbose = 0
#     ),
#     error = function(e) {
#       cat("\n XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
#       return(NULL)
#     }
#   )
#   
#   if (!is.null(xgb_cv)) {
#     
#     niter = xgb_cv$best_iteration
#     
#     xgb_model <- tryCatch(
#       xgb.train(
#         params = list(
#           objective = "reg:squarederror",
#           eval_metric = "rmse",
#           eta = 0.1,
#           max_depth = 6
#         ),
#         data = train_matrix,
#         watchlist = list(train = train_matrix),
#         nrounds = niter,
#         early_stopping_rounds = NULL,
#         verbose = 0
#       )
#     )
#     
#     
#     
#     predictions <- predict(xgb_model, newdata = test_matrix)
#     fitted_models$XGB[[cv_name]][[paste0("Fold_", i)]] <- list(
#       model = xgb_model,
#       predictions = predictions,
#       AIC = NA,
#       test_set = test_set


##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.6.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model




















#-------------------------- T1.7 : with evalution test  ----------------------------------
# 
# 
# 
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.7.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model









#-------------------------- T1.8 : : back to watchlist on test----------------------------------
# 
# 
# 
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.8.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model











#-------------------------- T1.9 : ------------------------------------------
# # Create dummy-encoded numeric matrices
# train_x <- model.matrix(~ . - 1, data = train_x_raw)
# test_x  <- model.matrix(~ . - 1, data = test_x_raw)
# 
# # --- Create a 10% validation split from train_x ---
# n <- nrow(train_x)
# val_idx <- sample(seq_len(n), size = floor(0.10 * n))
# 
# val_x <- train_x[val_idx, ]
# val_y <- train_df[[response_var]][val_idx]
# 
# train_x_sub <- train_x[-val_idx, ]
# train_y_sub <- train_df[[response_var]][-val_idx]
# 
# # Build DMatrix
# train_matrix <- 
#   xgb.DMatrix(data = train_x_sub, label = train_y_sub)
# 
# 
# val_matrix <- 
#   xgb.DMatrix(data = val_x, label = val_y)
# 
# 
# test_matrix <- 
#   xgb.DMatrix(data = test_x, label = test_df[[response_var]])
# 
# 
# 
# # Train and predict
# if (!is.null(train_matrix) && !is.null(test_matrix)) {
#   xgb_model <- tryCatch(
#     xgb.train(
#       params = list(
#         objective = "reg:squarederror",
#         eval_metric = "rmse",
#         eta = 0.1,
#         max_depth = 6
#       ),
#       data = train_matrix,
#       watchlist = list(test = val_matrix),
#       nrounds = 100,
#       early_stopping_rounds = 10,
#       verbose = 0
#     ),
#     error = function(e) {
#       cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
#       return(NULL)
#     }
#   )
##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.9.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model



#-------------------------- T1.10 : ------------------------------------------
# # Train and predict
# if (!is.null(train_matrix) && !is.null(test_matrix)) {
#   xgb_model <- tryCatch(
#     xgb.train(
#       params = list(
#         objective = "reg:squarederror",
#         eval_metric = "rmse",
#         eta = 0.1,
#         max_depth = 6
#       ),
#       data = train_matrix,
#       watchlist = list(train = train_matrix, test = test_matrix),
#       nrounds = 100,
#       # early_stopping_rounds = 10,
#       verbose = 0
#     ),
#     error = function(e) {
#       cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
#       return(NULL)
#     }
#   )


models_list <- readRDS("./output/models/xgboost/T1.10.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model













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



##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.10.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model






















#================================= Bloo with eval test ============================================
## bloo_cv_eval function ----
# Leave-One-Out CV with buffer exclusion
bloo_cv_eval <- function(spatial_data, buffer_size_km) {
  
  # Convert buffer size to meters
  buffer_size_m <- buffer_size_km * 1000
  
  # Initialize a list to store train-test_eval splits
  train_test_eval_splits <- list()
  
  # Extract coordinates
  coords <- st_coordinates(spatial_data)  # Extract spatial coordinates
  
  # Loop through each observation
  for (i in 1:nrow(spatial_data)) {
    
    # Define the test point (one point per iteration)
    test_set <- spatial_data[i, , drop = FALSE]  
    
    # Compute distances from the test point to all other points
    distances <- distVincentySphere(
      p1 = coords, # coords of all points
      p2 = coords[i, , drop = FALSE]  # coords of test point
    )
    
    # Identify remaining points (exclude points within the buffer)
    rmn_indices <- which(distances > buffer_size_m)  # Keep only points outside the buffer
    
    # Define remaining set 
    rmn_set <- spatial_data[rmn_indices, ]
    
    # Define the eval point in the remaining set
    eval_set <- rmn_set[which.max(distances[rmn_indices]), , drop = FALSE] # point furthest away from test point
    # eval_set <- rmn_set[sample(1:nrow(rmn_set), 1), , drop = FALSE] # random point
    
    # Compute distance from eval point to all other points in remaining set
    distances_eval <- distVincentySphere(
      p1 = st_coordinates(rmn_set), # coords of all points in remaining set
      p2 = st_coordinates(eval_set)  # coords of eval point
    )
    
    # Identify training points (exclude eval point and points within the buffer of eval point)
    train_indices <- which(distances_eval > buffer_size_m & rmn_indices != which.max(distances[rmn_indices]))
    
    # Define training set
    train_set <- rmn_set[train_indices, ]
    
    # Store train-test split
    train_test_eval_splits[[i]] <- list(train = train_set, test = test_set, eval = eval_set)
  }
  
  return(train_test_eval_splits)
}





## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# BLOO CV 50 with eval 

buffer_sizes <- c(50)  # only bloo50

# Initialize list to store BLOO-CV results
bloo_results_list <- list()

# Loop through each buffer size and compute BLOO-CV
for (buffer_size in buffer_sizes) {
  
  # Compute BLOO-CV
  bloo_results <- bloo_cv_eval(st_as_sf(tot_sf_4326), buffer_size)
  
  # Store results in list with a named index
  bloo_results_list[[paste0("bloo", buffer_size)]] <- bloo_results
}


# ALL CV
# Combine all CV configurations
cv_configs <- c( bloo_results_list)



rm(tot_sf_4326, bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Check and map CV config ----
# Extract BLOO with 50 km buffer
bloo50 <- cv_configs[["bloo50"]]

# Fold sizes: number of points in test (should be 1) and train
fold_sizes_bloo50 <- purrr::map_dfr(
  seq_along(bloo50),
  ~ tibble(
    Fold    = .x,
    n_test  = nrow(bloo50[[.x]]$test),
    n_train = nrow(bloo50[[.x]]$train),
    n_eval = nrow(bloo50[[.x]]$eval)
  )
)

print(head(fold_sizes_bloo50))

# Fold n_test n_train n_eval
# <int>  <int>   <int>  <int>
#   1     1      1     515      1
# 2     2      1     480      1
# 3     3      1     485      1
# 4     4      1     485      1
# 5     5      1     515      1
# 6     6      1     485      1


# Fold n_test n_train n_eval
# <int>  <int>   <int>  <int>
#   1     1      1     464      1
# 2     2      1     485      1
# 3     3      1     473      1
# 4     4      1     475      1
# 5     5      1     488      1
# 6     6      1     470      1

summary(fold_sizes_bloo50)

# Fold         n_test     n_train          n_eval 
# Min.   :  1   Min.   :1   Min.   :427.0   Min.   :1  
# 1st Qu.:160   1st Qu.:1   1st Qu.:485.0   1st Qu.:1  
# Median :319   Median :1   Median :506.0   Median :1  
# Mean   :319   Mean   :1   Mean   :503.9   Mean   :1  
# 3rd Qu.:478   3rd Qu.:1   3rd Qu.:531.0   3rd Qu.:1  
# Max.   :637   Max.   :1   Max.   :555.0   Max.   :1  



# Fold         n_test     n_train          n_eval 
# Min.   :  1   Min.   :1   Min.   :394.0   Min.   :1  
# 1st Qu.:160   1st Qu.:1   1st Qu.:463.0   1st Qu.:1  
# Median :319   Median :1   Median :484.0   Median :1  
# Mean   :319   Mean   :1   Mean   :482.4   Mean   :1  
# 3rd Qu.:478   3rd Qu.:1   3rd Qu.:504.0   3rd Qu.:1  
# Max.   :637   Max.   :1   Max.   :559.0   Max.   :1  



# --- Map folds ---
extract_fold_df <- function(cv_object, fold_id, keep_geometry = FALSE) {
  
  if (fold_id < 1 || fold_id > length(cv_object)) {
    stop("fold_id is out of range for this CV object.")
  }
  
  train_sf <- cv_object[[fold_id]]$train
  test_sf  <- cv_object[[fold_id]]$test
  eval_sf <- cv_object[[fold_id]]$eval
  
  if (!keep_geometry) {
    train_df <- train_sf |> sf::st_drop_geometry()
    test_df  <- test_sf  |> sf::st_drop_geometry()
    eval_df <- eval_sf |> sf::st_drop_geometry()
  } else {
    train_df <- train_sf
    test_df  <- test_sf
    eval_df <- eval_sf
  }
  
  train_df$set <- "train"
  test_df$set  <- "test"
  eval_df$set <- "eval"
  
  # merge into one data.frame
  df <- dplyr::bind_rows(test_df, train_df, eval_df)
  
  return(df)
}

for (X in 1:10) {

  X = 1 
  # Build the fold df + convert back to sf
  f <- extract_fold_df(bloo50, fold_id = X) %>% 
    st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)
  
  # Build plot
  p <- ggplot() +
    geom_sf(data = f, aes(color = set), size = 1) +
    scale_color_manual(values = c("train" = "blue", "test" = "red", "eval" = "yellow")) +
    labs(title = paste("BLOO CV", X, ": Train (blue), Test (red), Eval (yellow)")) +
    theme_minimal()
  
  # Print to RStudio plotting window
  print(p)
  
  # Save plot
   #filename <- sprintf("./figures/Cross-val/Map_Bloo50_EVAL_fold%d_T0.0.png", X)
  ggsave(filename, plot = p, width = 7, height = 6, dpi = 300)
}



# Panel plot of 12 randomly selected folds
set.seed(1)  # for reproducibility

# 1) Randomly select 12 distinct folds from BLOO CV
n_folds_total <- length(bloo50)
n_sample      <- 12

fold_ids_sample <- sample(seq_len(n_folds_total), size = n_sample)

# 2) Build combined sf object for these folds
panel_data <- map_dfr(fold_ids_sample, function(fid) {
  df <- extract_fold_df(bloo50, fold_id = fid)
  df$panel_fold <- paste0("Fold ", fid)  # label for facet
  df
})

panel_data_sf <- panel_data %>%
  st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)

# 3) Build 4 x 3 panel plot
p_panel <- ggplot(panel_data_sf) +
  geom_sf(aes(color = set), size = 0.8) +
  scale_color_manual(values = c("train" = "blue", "test" = "red", "eval" = "green")) +
  facet_wrap(~ panel_fold, ncol = 4) +
  labs(
    title    = "BLOO CV with EVAL - furthest distance (50 km) – Sample of 12 Folds",
    subtitle = "Train (blue), Test (red), Eval (green)",
    color    = "Set"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

# 4) Print to device
print(p_panel)

# 5) Save single figure
outfile <- "./figures/Cross-val/Map_Bloo50_12folds_T1._EVAL.png"
ggsave(outfile, plot = p_panel, width = 10, height = 8, dpi = 300)


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







#-------------------------- T1.11 : Eval test in bloo_cv_eval (random) ------------------------------------------

# Using bloo_cv_eval to compute bloo50 with eval test & fit_model_eval to use eval test : 
# xgb_model <- tryCatch(
#   xgb.train(
#     params = list(
#       objective = "reg:squarederror",
#       eval_metric = "rmse",
#       eta = 0.1,
#       max_depth = 6
#     ),
#     data = train_matrix,
#     watchlist = list(train = train_matrix, test = eval_matrix),
#     nrounds = 100,
#     early_stopping_rounds = 10,
#     verbose = 0
#   )

##  Evaluate performance ----
models_list <- readRDS("./output/models/xgboost/T1.11.rds")
perf <- evaluate_models(models_list)
perf
models_list$R$XGB$bloo50$Fold_1$model


















#-------------------------- T1.12 : With eval set from Gaétan  ------------------------------------------
# # Train and predict
# if (!is.null(train_matrix) && !is.null(test_matrix)) {
#   xgb_model <- tryCatch(
#     xgb.train(
#       params = list(
#         objective = "reg:squarederror",
#         eval_metric = "rmse",
#         eta = 0.1,
#         max_depth = 6
#       ),
#       data = train_matrix,
#       watchlist = list(train = train_matrix, test = test_matrix),
#       nrounds = 100,
#       # early_stopping_rounds = 10,
#       verbose = 0
#     ),
#     error = function(e) {
#       cat("\n❌ XGB failed for fold", i, "of", cv_name, ": ", e$message, "\n")
#       return(NULL)
#     }
#   )


models_list <- readRDS("./output/models/xgboost/T1.12.rds")
perf <- evaluate_models(models_list)
perf
saveRDS("./output/models/xgboost/T1.12-performance.rds")
models_list$R$XGB$bloo50$Fold_1$model


# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var
# 1 0.003095748  -0.05563945   -0.05822617 16.69027 20.37551  NA            R













#================================ GLM ============================================
# Resp var distribution ---------------------
#----------- CHECKS : normality + poisson --------------------------
# Explanation : We check if response variables are normally distributed in order to determine which distribution to use in glm.
# Documentation for normality tests : https://www.geeksforgeeks.org/how-to-test-for-normality-in-r/

# Initialize results dataframe
normality_results <- data.frame(Column = character(), 
                                Shapiro_p = numeric(), 
                                KS_p = numeric(), 
                                Normality = character(), 
                                stringsAsFactors = FALSE)
ind <- div[-1]
for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Compute normality tests
  shapiro_res <- shapiro.test(ind[[col]])
  
  # Handle ties issue in KS test
  ks_res <- tryCatch({
    ks.test(ind[[col]], "pnorm", mean(ind[[col]]), sd(ind[[col]]))
  }, warning = function(w) {
    cat("⚠️ Warning: Ties detected in KS test. Using jittered values.\n")
    ks.test(jitter(ind[[col]]), "pnorm", mean(ind[[col]]), sd(ind[[col]]))
  })
  
  ad_res <- ad.test(ind[[col]])   # Anderson-Darling test
  dag_res <- agostino.test(ind[[col]])  # D’Agostino test
  
  # Print results
  print(shapiro_res)
  print(ks_res)
  print(ad_res)
  print(dag_res)
  
  # Determine normality based on multiple tests
  normality_status <- ifelse(shapiro_res$p.value > 0.05 & ks_res$p.value > 0.05 & 
                               ad_res$p.value > 0.05 & dag_res$p.value > 0.05, 
                             "Normal", "Not-Normal")
  
  # Store results in dataframe
  normality_results <- rbind(normality_results, 
                             data.frame(Column = col, 
                                        Shapiro_p = shapiro_res$p.value, 
                                        KS_p = ks_res$p.value, 
                                        AD_p = ad_res$p.value, 
                                        DAg_p = dag_res$p.value, 
                                        Normality = normality_status))
}

# Print results
print(normality_results)

# Save results
write.csv(normality_results, "./output/div_indices_sel_v1.1_normality_test.csv", row.names = FALSE)





#----------- FIGURE : Histograms + Normality test ---------------------
# Function to create histogram with Shapiro-Wilk normality test
create_hist_plot <- function(var_name, values) {
  
  # Remove NAs to avoid errors
  values <- na.omit(values)
  
  # Perform Shapiro-Wilk Normality Test
  shapiro_test <- shapiro.test(values)
  p_value <- shapiro_test$p.value
  interpretation <- ifelse(p_value >= 0.05, "Normal", "Not Normal")
  
  # Create histogram
  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_histogram(bins = 15, fill = "gray", color = "black", alpha = 0.6) +
    labs(title = var_name, x = "", y = "Frequency") +
    theme_minimal(base_size = 12) +
    annotate("text", x = max(values, na.rm = TRUE) * 0.7, 
             y = max(hist(values, plot = FALSE)$counts) * 0.9,
             label = paste0("SW p = ", round(p_value, 3), "\n(", interpretation, ")"),
             size = 4, hjust = 0, vjust = 1, color = "black", fontface = "bold",
             label.size = 0.5)
  
  return(p)
}


# Generate histograms for each variable in 'ind'
plots <- lapply(names(ind), function(var) create_hist_plot(var, ind[[var]]))

# Arrange all plots in a grid
grid.arrange(grobs = plots, ncol = 3, top = "Histograms with Shapiro-Wilk Normality Test Results")




#----------- FIGURE : QQ Plot + Poisson test --------------------------
# See doc here on how we check if data is Poisson distributed : https://www.geeksforgeeks.org/how-to-know-if-a-data-follows-a-poisson-distribution-in-r/

# Function to create QQ plot for Poisson distribution with chi-squared test
create_qq_plot <- function(var_name, values) {
  
  # Remove NA values
  values <- na.omit(values)
  
  # Convert non-integer values to nearest integers
  values <- round(values)
  
  # Compute mean & variance
  lambda <- mean(values)
  variance <- var(values)
  
  # Observed frequencies
  obs_freq <- table(values)
  
  # Expected frequencies based on Poisson distribution
  exp_freq <- dpois(as.numeric(names(obs_freq)), lambda) * length(values)
  
  # Perform Chi-Squared Test with Monte Carlo approximation
  chisq_test <- tryCatch(
    chisq.test(obs_freq, p = exp_freq, rescale.p = TRUE, simulate.p.value = TRUE, B = 10000),
    error = function(e) return(NULL)
  )
  
  # Interpretation
  if (!is.null(chisq_test)) {
    p_value <- chisq_test$p.value
    interpretation <- ifelse(p_value >= 0.05, "Poisson Distributed", "Not Poisson Distributed")
  } else {
    p_value <- NA
    interpretation <- "Chi-Squared Test Failed"
  }
  
  # Generate QQ Plot
  qqplot_data <- data.frame(
    x = qpois(ppoints(length(values)), lambda),
    y = sort(values)
  )
  
  # Fixed annotation position using `Inf` (top-right of plot)
  p <- ggplot(qqplot_data, aes(x = x, y = y)) +
    geom_point(color = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = var_name, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal(base_size = 12) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Mean: ", round(lambda, 2), 
                            "\nVar: ", round(variance, 2),
                            "\nChi-sq p = ", ifelse(is.na(p_value), "NA", round(p_value, 3)),
                            "\n(", interpretation, ")"),
             hjust = 1, vjust = 1, size = 4, fontface = "bold", color = "black")
  
  return(p)
}






# Generate QQ plots for each variable in 'ind'
qq_plots <- lapply(names(ind), function(var) create_qq_plot(var, ind[[var]]))

# Arrange all QQ plots in a grid using patchwork
final_plot <- patchwork::wrap_plots(qq_plots, ncol = 3) +
  patchwork::plot_annotation(title = "QQ Plots for Poisson Distribution with Chi-Squared Test")

# Display plot
print(final_plot)



#----------- CHECKS : 0 inflation --------------------------
# Load required package
if (!require("pscl")) install.packages("pscl", dependencies = TRUE)
library(pscl)

# Initialize results dataframe
zero_inflation_results <- data.frame(Column = character(), 
                                     Zero_Proportion = numeric(), 
                                     Vuong_p = numeric(), 
                                     Zero_Inflated = character(), 
                                     stringsAsFactors = FALSE)

# Loop over each column
for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Count zeros
  zero_count <- sum(ind[[col]] == 0, na.rm = TRUE)
  zero_proportion <- zero_count / length(ind[[col]])
  
  # Perform Vuong test for zero inflation in Poisson regression
  tryCatch({
    poisson_model <- glm(ind[[col]] ~ 1, family = poisson)
    zip_model <- zeroinfl(ind[[col]] ~ 1, dist = "poisson")
    vuong_test <- vuongs(poisson_model, zip_model)  # Vuong test
    
    # Extract p-value
    vuong_p <- as.numeric(vuong_test$p.value["Raw"])
  }, error = function(e) {
    cat("⚠️ Error: Vuong test failed for", col, "\n")
    vuong_p <<- NA
  })
  
  # Print results
  cat("   - Zero Proportion:", round(zero_proportion, 3), "\n")
  cat("   - Vuong Test p-value:", ifelse(is.na(vuong_p), "Failed", round(vuong_p, 3)), "\n")
  
  # Determine zero inflation based on Vuong test (p < 0.05 suggests zero inflation)
  zero_inflation_status <- ifelse(!is.na(vuong_p) & vuong_p < 0.05, "Zero-Inflated", "Not Zero-Inflated")
  
  # Store results in dataframe
  zero_inflation_results <- rbind(zero_inflation_results, 
                                  data.frame(Column = col, 
                                             Zero_Proportion = zero_proportion, 
                                             Vuong_p = vuong_p, 
                                             Zero_Inflated = zero_inflation_status))
}

# Print final results
print(zero_inflation_results)

#----------- CHECKS : Overdispersion --------------------------
check_overdispersion <- function(dataset, response_vars, predictors) {
  
  # Store all plots
  plots <- list()
  
  # Loop through each response variable
  for (response_var in response_vars) {
    
    cat("\n🔹 Checking Overdispersion for:", response_var, "\n")
    
    # Define formula
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    
    # Fit Poisson Model
    poisson_model <- glm(formula, data = dataset, family = poisson(link = "log"))
    
    # Run Overdispersion Test
    aer_test <- AER::dispersiontest(poisson_model)
    dispersion_ratio <- aer_test$estimate["dispersion"]
    p_value <- aer_test$p.value
    overdisp_status <- ifelse(p_value < 0.05, "Significant Overdispersion", "No Strong Evidence of Overdispersion")
    
    # Generate Residual Plots
    residual_data <- data.frame(
      Fitted = poisson_model$fitted.values,
      Residuals = residuals(poisson_model, type = "pearson")
    )
    
    p <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "blue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Overdispersion Check -", response_var), 
           x = "Fitted Values", y = "Pearson Residuals") +
      theme_minimal(base_size = 12) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold",
               label = paste0(
                 "Dispersion: ", round(dispersion_ratio, 2), 
                 "\nP-value: ", round(p_value, 3),
                 "\n", overdisp_status
               ))
    
    plots[[response_var]] <- p
  }
  
  # Arrange all plots in a grid
  grid.arrange(grobs = plots, ncol = 2, top = "Overdispersion Analysis Across Variables")
}
check_overdispersion <- function(dataset, response_vars, predictors) {
  
  # Store all plots
  plots <- list()
  
  # Initialize a data frame to store results
  results_df <- data.frame(
    Response_Variable = character(),
    Dispersion_Index = numeric(),
    P_Value = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each response variable
  for (response_var in response_vars) {
    
    cat("\n🔹 Checking Overdispersion for:", response_var, "\n")
    
    # Define formula
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    
    # Fit Poisson Model
    poisson_model <- glm(formula, data = dataset, family = poisson(link = "log"))
    
    # Run Overdispersion Test
    aer_test <- AER::dispersiontest(poisson_model)
    dispersion_ratio <- aer_test$estimate["dispersion"]
    p_value <- aer_test$p.value
    overdisp_status <- ifelse(p_value < 0.05, "Significant Overdispersion", "No Strong Evidence of Overdispersion")
    
    # Save results in the data frame
    results_df <- rbind(results_df, data.frame(
      Response_Variable = response_var,
      Dispersion_Index = round(dispersion_ratio, 2),
      P_Value = round(p_value, 3),
      Interpretation = overdisp_status,
      stringsAsFactors = FALSE
    ))
    
    # Generate Residual Plots
    residual_data <- data.frame(
      Fitted = poisson_model$fitted.values,
      Residuals = residuals(poisson_model, type = "pearson")
    )
    
    p <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "blue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Overdispersion Check -", response_var), 
           x = "Fitted Values", y = "Pearson Residuals") +
      theme_minimal(base_size = 12) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold",
               label = paste0(
                 "Dispersion: ", round(dispersion_ratio, 2), 
                 "\nP-value: ", round(p_value, 3),
                 "\n", overdisp_status
               ))
    
    plots[[response_var]] <- p
  }
  
  # Print results table
  print(results_df)
  
  # Arrange all plots in a grid
  grid.arrange(grobs = plots, ncol = 2, top = "Overdispersion Analysis Across Variables")
  
  # Return the results dataframe
  return(results_df)
}

# Usage
check_overdispersion(dataset = tot, response_vars = colnames(ind), predictors = predictors)


#----------- CHECKS : mean and variance --------------------------
# Initialize results dataframe
mean_variance_results <- data.frame(Column = character(), 
                                    Mean = numeric(), 
                                    Variance = numeric(), 
                                    Overdispersed = character(), 
                                    stringsAsFactors = FALSE)

# Loop over each column
for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Compute mean & variance
  mean_val <- mean(ind[[col]], na.rm = TRUE)
  var_val <- var(ind[[col]], na.rm = TRUE)
  
  # Print results
  cat("   - Mean:", round(mean_val, 3), "\n")
  cat("   - Variance:", round(var_val, 3), "\n")
  
  # Determine overdispersion (Variance > Mean)
  overdispersion_status <- ifelse(var_val > mean_val, "Yes", "No")
  
  # Store results in dataframe
  mean_variance_results <- rbind(mean_variance_results, 
                                 data.frame(Column = col, 
                                            Mean = mean_val, 
                                            Variance = var_val, 
                                            Overdispersed = overdispersion_status))
}

# Print final results
print(mean_variance_results)



#----------- DOC : quasipoisson vs nb distribution ---------------------
# how do I now if my overdispersion is "mild" or "severe" ?
#   ChatGPT said:
#   
#   ​Overdispersion in Poisson regression occurs when the observed variance in the count data exceeds the mean, violating the Poisson model's assumption that the mean equals the variance. Assessing the severity of overdispersion is crucial for selecting an appropriate statistical model. Here's how you can evaluate it:​
# en.wikipedia.org
# 
# 1. Calculate the Dispersion Statistic:
#   
#   The dispersion statistic (also known as the dispersion parameter or index of dispersion) is calculated as:​
# 
# Dispersion Statistic=∑(Yi−Y^i)2/Y^in−pDispersion Statistic=n−p∑(Yi​−Y^i​)2/Y^i​​​
# 
# where YiYi​ are the observed counts, Y^iY^i​ are the predicted counts from the model, nn is the number of observations, and pp is the number of parameters estimated in the model.​
# 
# 2. Interpret the Dispersion Statistic:
#   
#   Close to 1: Indicates that the variance is approximately equal to the mean, suggesting no overdispersion.​
# 
# Greater than 1: Indicates overdispersion.​
# 
# 3. Assess the Severity of Overdispersion:
#   
#   Mild Overdispersion: Dispersion statistic is slightly above 1 (e.g., between 1 and 2). In such cases, the standard errors may be slightly underestimated, but the Poisson model might still be adequate with robust standard errors.​
# 
# Moderate to Severe Overdispersion: Dispersion statistic significantly exceeds 2. This suggests that the Poisson model is likely inappropriate, and alternative modeling approaches should be considered.​
# 
# 4. Alternative Modeling Approaches:
#   
#   Quasi-Poisson Model: Adjusts the standard errors to account for overdispersion without altering the mean structure. Suitable for mild to moderate overdispersion.​
# 
# Negative Binomial Model: Introduces an additional parameter to model the overdispersion explicitly. Preferred for moderate to severe overdispersion.​
# 
# 5. Formal Testing:
#   
#   While calculating the dispersion statistic provides a straightforward assessment, formal statistical tests, such as Dean's test for overdispersion, can also be employed to evaluate the presence and extent of overdispersion. ​



#----------- DOC : over and underdispersion issues --------------------------
# is underdispersion an issue for normally distributed variables ? 
#   ChatGPT said:
#   🔍 Is Underdispersion an Issue for Normally Distributed Variables?
#   
#   No, underdispersion is not an issue for normally distributed variables. The concept of underdispersion primarily applies to count data models (Poisson, Negative Binomial, etc.), where variance is expected to follow a specific structure (e.g., equal to the mean in Poisson models).

# 🔹 Why Doesn't Underdispersion Matter for Normal Distributions?
# 
#     Normal Distribution Assumes Constant Variance
#         In a normal distribution:
#         Y∼N(μ,σ2)
#         Y∼N(μ,σ2) variance (σ2σ2) is a free parameter that is estimated from the data.
#         There's no requirement that variance follows the mean, unlike in Poisson models.
# Whether variance is "too small" or "too large" relative to expectations does not break any assumptions.
# 
# Overdispersion and Underdispersion are Concerns for Poisson-Based Models
# In Poisson models, variance must equal the mean (λλ).
# Underdispersion (Var<MeanVar<Mean) causes issues because standard Poisson models cannot accommodate it.
# Overdispersion (Var>MeanVar>Mean) causes similar problems.
# This is why alternative models (Quasi-Poisson, Negative Binomial, Generalized Poisson) exist.
# 
# In Normal Regression (GLM, GAM, etc.), Variance is Simply Estimated
# In a Gaussian GLM/GAM (family = gaussian()), we assume:
#   Y=Xβ+ε,ε∼N(0,σ2)
# Y=Xβ+ε,ε∼N(0,σ2)
# The model does not expect variance to follow a specific function of the mean, so underdispersion is not problematic.



# 🔍 Is Overdispersion an Issue for Normally Distributed Variables?
#   
#   No, overdispersion is not an issue for normally distributed variables because:
#   
#   The Normal (Gaussian) distribution does not assume a strict relationship between the mean and variance.
# Unlike Poisson or Negative Binomial models, where variance is linked to the mean, a Gaussian model explicitly estimates variance (σ2σ2) as a separate parameter.




# 
# is underdispersion or overdispersion an issue for negative binomial distribution ? 
#   ChatGPT said:
#   🔍 Is Underdispersion or Overdispersion an Issue for Negative Binomial (NB) Distribution?
#   🔹 Overdispersion: NO, NB is designed for it
# 
# Overdispersion is NOT an issue for the Negative Binomial distribution because it is explicitly designed to handle count data where variance exceeds the mean.
# 🔹 Underdispersion: YES, NB is NOT suited for it
# 
# Underdispersion is a problem because the Negative Binomial assumes variance is always greater than the mean:
#   Var(Y)=μ+κμ2
# Var(Y)=μ+κμ2
# 
# where κ>0κ>0 controls the dispersion. This means the variance must be greater than the mean. If the observed variance is less than expected, the NB model will likely perform poorly.
# 🔹 When is Overdispersion or Underdispersion a Problem?
#   Dispersion Type	Impact on NB Model	Solution
# Overdispersion (Var>MeanVar>Mean)	✅ NB models this well	Use glm.nb() or mgcv::nb()
# Underdispersion (Var<MeanVar<Mean)	❌ NB does not fit well	Use Generalized Poisson or COM-Poisson





#----------- RUN : Assign distribution to each resp var --------------------------

# Create a df to store the distribution of each response variable
distrib <- data.frame(response_var = colnames(ind))

# Assign distribution for each response variable : 
# poisson for count variables,
# gaussian for continuous variables,
# nb for overdispersed count variables,
# (using nb for overdispersion because quasipoisson is better suited for "mild" overdispersion (ie between 1 and 2) and nb for "severe" overdispersion (ie > 2))



distrib$Distribution <- c("nb", "nb", "nb",  "gaussian", "nb",  "gaussian", "nb", "nb", "poisson", "poisson")
distrib

names(ind)

# "R" -->  nb (overdispersed count var)
# "Crypto"  -->  nb (overdispersed count var)
# "Elasmo"  -->  nb (overdispersed count var)
# "DeBRa" --> gaussian ? 
# "RedList"  -->  nb (overdispersed count var)
# "LRFI"   -->  gaussian (continuous count var)
# "TopPred"  -->  nb (overdispersed count var)
# "Commercial -->  nb (overdispersed count var)
# "Grouper"   --> poisson
# "AngelShark" --> poisson









#----------- check model assumptions -------------------------


check_model_assumptions <- function(models_list) {

require(ggplot2)
require(patchwork)
require(lmtest)
require(car)
require(ggpubr)
require(grid)
require(mgcv)
require(lme4)
require(spaMM)  # ✅ Ensure spaMM is loaded

diagnostic_results <- data.frame(
  Response = character(),
  Model_Type = character(),
  Shapiro_p = numeric(),
  Poisson_p = numeric(),
  Spearman_p = numeric(),
  Spearman_rho = numeric(),
  BP_p = numeric(),
  Levene_p = numeric(),
  Dispersion_Index = numeric(),
  Dispersion_p = numeric(),
  Dispersion_Interpretation = character(),
  Normality = character(),
  Poisson_Distribution = character(),
  Homoscedasticity = character(),
  Linearity = character(),
  stringsAsFactors = FALSE
)

all_plots <- list()

for (response_var in names(models_list)) {
  model_obj <- models_list[[response_var]]
  
  if (!is.list(model_obj)) next # Skip if model object is not a list
  
  # Detect model type
  if ("GLM" %in% names(model_obj)) {
    
    model <- model_obj$GLM
    model_type <- "GLM"
  } else if ("GAM" %in% names(model_obj)) {
    
    model <- model_obj$GAM
    model_type <- "GAM"
  } else if ("spaMM" %in% names(model_obj)) {
    
    model <- model_obj$spaMM  # ✅ Extract correct spaMM model
    model_type <- "spaMM"
  } else {
    next
  }
  
  # ✅ Fix: Extract residuals correctly for spaMM models
  if (inherits(model, "HLfit")) {  # ✅ Check if it's a spaMM model
    residuals <- resid(model, type = "std_dev_res")  # ✅ Use resid() WITHOUT type argument
    fitted_values <- predict(model, type = "response")  # ✅ Extract fitted values properly
  } else {
    residuals <- residuals(model, type = "pearson")
    fitted_values <- fitted(model)
  }
  
  # ✅ Fix: Handle standard residuals correctly
  if (inherits(model, "gam")) {
    std_residuals <- residuals(model, type = "deviance")  
  } else if (inherits(model, "HLfit")) {  # ✅ Handle spaMM models
    std_residuals <- resid(model)  # ✅ Use resid() without type argument
  } else {
    std_residuals <- rstandard(model)  
  }
  
  shapiro_test <- shapiro.test(residuals)
  normality_status <- ifelse(shapiro_test$p.value > 0.05, "Normal", "Not Normal")
  
  spearman_test <- cor.test(fitted_values, residuals, method = "spearman")
  linearity_status <- ifelse(abs(spearman_test$estimate) < 0.2, "Linear", "Non-Linear")
  
  if (model_type == "GLM" || model_type == "spaMM") {
    bp_test <- bptest(model)
    homoscedasticity_status <- ifelse(bp_test$p.value > 0.05, "Homoscedastic", "Heteroscedastic")
    levene_p_value <- NA
  } else {
    fitted_groups <- cut(fitted_values, breaks = 5, include.lowest = TRUE)
    levene_test <- tryCatch({
      leveneTest(residuals, as.factor(fitted_groups))
    }, error = function(e) return(NULL))
    
    if (!is.null(levene_test) && is.numeric(levene_test$`Pr(>F)`[1])) {
      levene_p_value <- round(levene_test$`Pr(>F)`[1], 4)
      homoscedasticity_status <- ifelse(levene_p_value > 0.05, "Homoscedastic", "Heteroscedastic")
    } else {
      levene_p_value <- NA
      homoscedasticity_status <- "Test Not Valid"
    }
    
    bp_test <- list(p.value = NA)
  }
  
  
  
  # Compute df_residual manually if missing
  if (is.null(model$df.residual)) {
    df_residual <- nobs(model) - length(fixef(model))  # Estimate df.residual manually
  } else {
    df_residual <- model$df.residual
  }
  
  # Ensure df_residual is numeric and valid
  if (is.na(df_residual) || df_residual <= 0) {
    dispersion_index <- NA
    dispersion_p <- NA
    dispersion_interpretation <- "Invalid df_residual"
  } else {
    # Compute dispersion index
    dispersion_index <- sum(residuals^2) / df_residual
    
    # Compute dispersion p-value
    dispersion_p <- pchisq(dispersion_index * df_residual, df_residual, lower.tail = FALSE)
    
    # Interpretation of dispersion
    dispersion_interpretation <- ifelse(
      dispersion_index < 1, "Underdispersion",
      ifelse(dispersion_index > 2, "Overdispersed", "No Overdispersion")
    )
  }
  
  diagnostic_results <- rbind(diagnostic_results, 
                              data.frame(
                                Response = response_var,
                                Model_Type = model_type,
                                Shapiro_p = round(shapiro_test$p.value, 4),
                                Poisson_p = NA,
                                Spearman_p = round(spearman_test$p.value, 4),
                                Spearman_rho = round(spearman_test$estimate, 4),
                                BP_p = round(bp_test$p.value, 4),
                                Levene_p = levene_p_value,
                                Dispersion_Index = round(dispersion_index, 4),
                                Dispersion_p = round(dispersion_p, 4),
                                Dispersion_Interpretation = dispersion_interpretation,
                                Normality = normality_status,
                                Poisson_Distribution = NA,
                                Homoscedasticity = homoscedasticity_status,
                                Linearity = linearity_status
                              ))
  
  # ✅ Generate diagnostic plots
  p1 <- ggplot(data.frame(Residuals = residuals), aes(x = Residuals)) +
    geom_histogram(bins = 15, fill = "deepskyblue", color = "black", alpha = 0.6, aes(y = after_stat(density))) +
    geom_density(color = "red", linewidth = 1) +
    labs(title = paste(response_var, "- Residuals Histogram")) +
    theme_minimal()
  
  p2 <- ggplot(data.frame(Fitted = fitted_values, Residuals = residuals), aes(x = Fitted, y = Residuals)) +
    geom_point(color = "deepskyblue", alpha = 0.6) +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    labs(title = paste(response_var, "- Residuals vs Fitted")) +
    theme_minimal()
  
  all_plots[[response_var]] <- (p1 | p2)
  
  print((p1 | p2))  # ✅ Explicitly print each plot
}

print(diagnostic_results)  # ✅ Explicitly print the diagnostics table

return(list(Diagnostics = diagnostic_results, Plots = all_plots))
}

#-------------------------------- T1.0 : Fit GLM models : bloo50 ------------------------------------------
## Description ----
#--- Model : GLM

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE      AIC Response_Var Model     CV Train_Size
# 1 0.103464    0.3216582     0.3449345 13.79105 17.14771 5121.264            R   GLM bloo50         NA




## Prep data for model ----
# Assign distribution to each response variable
distrib <- data.frame(response_var = colnames(ind))
distrib$Distribution <- c("nb")



# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select(("R"))

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
#loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
# buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
cv_configs <- c( bloo_results_list)



rm(tot_sf_4326, bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("GLM"), 
    distribution = distrib
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}


## Save model ----

saveRDS(models_list, "./output/models/glm/T1.0.rds")


##  Evaluate performance ----
models_list <- readRDS("./output/models/glm/T1.0.rds")
perf <- evaluate_models(models_list)
models_list$R$XGB$bloo50$Fold_1$model

# Save perf 
saveRDS(perf, "./output/models/glm/T1.0s_performance.rds")

















## Check model assumption ----
residual_results_glm <- check_model_assumptions(models_list)

class(models_list)
str(models_list, max.level = 3)


#================================ GAM ============================================
#-------------------------------- T1.0 : Fit GAM models : bloo50 ------------------------------------------
## Description ----
#--- Model : GAM

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE     AIC Response_Var Model     CV Train_Size
# 1 0.1012086     0.318133     0.3413883 13.78959 17.13608 4669.13            R   GAM bloo50         NA



## Prep data for model ----
# Assign distribution to each response variable
distrib <- data.frame(response_var = colnames(ind))
distrib$Distribution <- c("nb")



# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select(("R"))

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
#loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
# buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
cv_configs <- c( bloo_results_list)



rm(tot_sf_4326, bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("GAM"), 
    distribution = distrib
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}


## Save model ----
saveRDS(models_list, "./output/models/gam/T1.0.rds")


##  Evaluate performance ----
models_list <- readRDS("./output/models/gam/T1.0.rds")
perf <- evaluate_models(models_list)

# Save perf 
saveRDS(perf, "./output/models/gam/T1.0s_performance.rds")




















#================================ GLS RF ============================================

# Functions for parallelisation + evaluation on GLS -----

# Function : fit_models_parallel_glsrf ------
fit_models_parallel_glsrf <- function(dataset, response_var, predictors, cv_configs = NULL,
                                      models = c("GLSRF"),
                                      distribution = NULL, mtry = NULL) {
  # Packages (you can also load these once outside the function)
  library(future.apply)
  library(RandomForestsGLS)
  
  cat("\n -------------------", "\n Starting Model Evaluation for:", response_var, "\n")
  
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
    
    fold_ids <- seq_along(train_test_splits)
    
    # ---- PARALLEL OVER FOLDS ----
    fold_results <- lapply(
      fold_ids,
      function(i) {
        train_set <- train_test_splits[[i]]$train
        test_set  <- train_test_splits[[i]]$test
        
        train_df <- as.data.frame(sf::st_drop_geometry(train_set))
        test_df  <- as.data.frame(sf::st_drop_geometry(test_set))
        
        if (!(response_var %in% colnames(train_df))) {
          stop(paste("Response variable", response_var, "is missing!"))
        }
        
        formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
        
        fold_out <- list()
        
        # --- GLSRF ---
        if ("GLSRF" %in% models) {
          fold_out$GLSRF <- tryCatch({
            
            coords_train <- sf::st_coordinates(train_set)
            coords_test  <- sf::st_coordinates(test_set)
            
            X_train <- as.matrix(train_df[, predictors, drop = FALSE])
            X_test  <- as.matrix(test_df[, predictors, drop = FALSE])
            y_train <- train_df[[response_var]]
            
            glsrf_model <- RandomForestsGLS::RFGLS_estimate_spatial(
              coords = coords_train,
              y      = y_train,
              X      = X_train,
              mtry   = 2
            )
            
            glsrf_pred_mean <- RandomForestsGLS::RFGLS_predict(
              RFGLS_out = glsrf_model,
              Xtest     = X_test
            )
            
            glsrf_pred_spatial <- RandomForestsGLS::RFGLS_predict_spatial(
              RFGLS_out = glsrf_model,
              coords.0  = coords_test,
              Xtest     = X_test
            )
            
            list(
              model            = glsrf_model,
              predictions      = glsrf_pred_spatial$prediction,  # spatial response
              predictions_mean = glsrf_pred_mean$predicted,      # mean function
              AIC              = NA,
              test_set         = test_set
            )
            
          }, error = function(e) {
            cat("\n❌ GLSRF failed for fold", i, "of", cv_name, ":", conditionMessage(e), "\n")
            NULL
          })
        }
        
        fold_out
      }#,
      #future.seed = TRUE
    )
    
    # ---- COLLECT RESULTS BACK INTO fitted_models ----
    for (i in fold_ids) {
      for (model_name in models) {
        if (!is.null(fold_results[[i]][[model_name]])) {
          fitted_models[[model_name]][[cv_name]][[paste0("Fold_", i)]] <-
            fold_results[[i]][[model_name]]
        }
      }
    }
  }
  
  return(fitted_models)
}




# Function : evaluate_models_gls_gaetan ----------------
evaluate_models_gls_gaetan <- function(models_list,
                            response_var = NULL,
                            models       = NULL,
                            cv_names     = NULL,   # names of CV methods to use (e.g. "env_k6", "bloo50")
                            cv_splits    = NULL,   # optional: list of CV splits to compute train sizes
                            output       = c("summary", "fold", "plot", "all")) {
  
  output <- match.arg(output)
  cat("\n🔹 Starting Model Evaluation\n")
  
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
      
      cat("\n🔹 Evaluating model:", model_name, "for response variable:", resp_var, "\n")
      
      available_cvs <- names(models_list[[resp_var]][[model_name]])
      
      cv_methods <- if (!is.null(cv_names)) {
        intersect(cv_names, available_cvs)
      } else {
        available_cvs
      }
      
      for (cv_name in cv_methods) {
        
        cat("\n🔹 Processing CV Method:", cv_name, "\n")
        
        model_folds <- models_list[[resp_var]][[model_name]][[cv_name]]
        # Optional CV splits to get train sizes
        this_cv_split <- if (!is.null(cv_splits)) cv_splits[[cv_name]] else NULL
        
        all_observed    <- c()
        all_predictions <- c()
        aic_values      <- c()
        train_sizes     <- c()
        
        for (i in seq_along(model_folds)) {
          
          cat("   ➡️ Fold:", i, "/", length(model_folds), "\n")
          
          fold_name <- paste0("Fold_", i)
          fold_data <- model_folds[[fold_name]]
          
          
          test_set <- fold_data$test_set
          
          
          test_df <- as.data.frame(st_drop_geometry(test_set))
          
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
            } else if (model_name == "GLSRF") {
              fold_pred <- fold_data$predictions
              
            } else {
              predict(fold_data$model, newdata = test_df, type = "response")
            }
          }, error = function(e) {
            cat("      ❌ Prediction failed for", model_name, "on", fold_name, ":", e$message, "\n")
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
          if (!is.null(this_cv_split) && length(this_cv_split) >= i && !is.null(this_cv_split[[i]]$train)) {
            this_train_size <- nrow(this_cv_split[[i]]$train)
            train_sizes     <- c(train_sizes, this_train_size)
          }
          
          #-----------------------------
          # Per-fold metrics
          #-----------------------------
          if (!all(is.na(fold_pred))) {
            fold_metrics <- compute_metrics(fold_observed, fold_pred,
                                            if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM"))
                                              AIC(fold_data$model) else NULL)
            
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
            cat("      ⚠️ All predictions NA for", fold_name, "- no fold metrics.\n")
          }
        } # end folds
        
        #-----------------------------
        # Global (across folds) metrics
        #-----------------------------
        if (length(all_observed) == 0 || length(all_predictions) == 0) {
          cat("\n⚠️ No valid predictions for model", model_name, "and CV", cv_name, "- skipping global metrics.\n")
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
              
              p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Observed, y = Predicted)) +
                ggplot2::geom_point(alpha = 0.6) +
                ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
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
        }
      }
    }
  }
  
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

#-------------------------- T1.0 GLS RF -----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.4 
# [1] "tpi_mean_log"               "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [6] "temp_mean_1m"               "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [11] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"     "year"                      
# [16] "soft_bottom_clr"            "rock_clr"     

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model
# 1 0.0498418    0.2232528       0.25861 13.97388 16.70465  NA            R GLSRF


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() 
                       %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x, y))), 
              "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)



# [MARBEC-GPU]
#================================= RF =============================================


#-------------------------- T1.0 : RF : bloo50 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry set to default --> thus predictors / 3

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2220545    0.4712266     0.4580655 12.49707 15.15105  NA            R    RF bloo50         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}
beepr::beep()

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.0.rds") 


##  Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.0.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.0_performance.rds")
perf
models_list$R$XGB$bloo50$Fold_1$model







#-------------------------- T1.1 : RF : bloo50 + sqrt(R) ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry set to default --> thus predictors / 3

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# sqrt(R)

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2195448     0.468556     0.4572552 12.50877 15.16944  NA            R    RF bloo50         NA


## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R") %>% mutate(R = sqrt(R))



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}
beepr::beep()

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.1.rds") 


##  Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.1.rds")
perf <- evaluate_models(models_list)
beepr::beep()

# Save perf
saveRDS(perf, "./output/models/rf/T1.1_performance.rds")
perf
models_list$R$XGB$bloo50$Fold_1$model


















#-------------------------- T1.2 : RF : bloo50 + log(R+1) ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry set to default --> thus predictors / 3

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# log(R+1)

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2373872     0.487224     0.4772086 12.50114 15.09188  NA            R    RF bloo50         NA



## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R") %>% mutate(R = log(R+1))



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}
beepr::beep()

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.2.rds") 


## Evaluate performance ----
# models_list <- readRDS("./output/models/rf/T1.2.rds")
perf <- evaluate_models(models_list)
beepr::beep()
# save
saveRDS(perf, "./output/models/rf/T1.2_performance.rds")




















#-------------------------- T1.3 : RF : bloo50, mtry = 2 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2419331     0.491867     0.4800872 12.46602 15.1112  NA            R    RF bloo50         NA


## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.3.rds") 


##  Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.3.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.3_performance.rds")

perf
beepr::beep()







#-------------------------- T1.4 : RF : bloo50, mtry = 3 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 3

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2311374    0.4807675     0.4675131 12.51519 15.13758  NA            R    RF bloo50         NA

## Prep data for model ----
source("./utils/Fct_Modeling.R")
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.4.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.3.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.4_performance.rds")
perf
beepr::beep()






#-------------------------- T1.5 : RF : bloo50, mtry = 1 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2397627    0.4896557     0.4797954 12.66045 15.26152  NA            R    RF bloo50         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.5.rds") 


## Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.5.rds")
perf <- evaluate_models(models_list)

# save
saveRDS(perf, "./output/models/rf/T1.5_performance.rds")

perf
beepr::beep()







#-------------------------- T1.6 : RF : bloo50, mtry = 2, ntree = 1000 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2387987    0.4886703     0.4779052 12.54107 15.13517  NA            R    RF bloo50         NA


## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.6.rds") 


##  Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.6.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.6_performance.rds")

perf
beepr::beep()









#-------------------------- T1.7: RF : bloo50, mtry = 2, ntree = 800 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                      

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2359166    0.4857125     0.4721527 12.54533 15.16021  NA            R    RF bloo50         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.7.rds") 


##  Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.7.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.7_performance.rds")

perf
beepr::beep()









#-------------------------- T1.8 : RF : bloo50, mtry = 2, ntree= 500, bloo20 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 20 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.334184    0.5780865     0.5637794 11.51268 14.24277  NA            R    RF bloo20         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 20

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo20"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.8.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.3.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.3_performance.rds")

perf
beepr::beep()








#-------------------------- T1.9 : RF : bloo50, mtry = 2, ntree= 500, bloo10 ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"                                       <== CHANGE

#---- CV config :
# BLOO CV with buffer size = 10 km

# #---- Perf : 
# R2 Pearson_Corr Spearman_Corr     MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.371679    0.6096548      0.595759 11.1009 13.81871  NA            R    RF bloo10         NA


## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 10

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo10"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.9.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.3.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.9_performance.rds")

perf
beepr::beep()









#-------------------------- T1.10 : RF : bloo50 with coordinates ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" <== CHANGE


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2570363    0.5069875     0.4974259 12.44734 15.0298  NA            R    RF bloo50         NA

## Prep data for model ----
# Extract predictors
predictors <- pred %>% dplyr::select(-c(replicates, grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.10.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.10.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.10_performance.rds")

perf
beepr::beep()












#-------------------------- T1.11 : RF : bloo50 with coordinates + year ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                 <== CHANGE


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.281792    0.5308408     0.5186794 12.23032 14.84524  NA            R    RF bloo50         NA

## Prep data for model ----
tot <- tot %>%
  mutate(year = year(date))



col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat))), "year")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 


# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.11.rds") 


## Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.11.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.11_performance.rds")

perf
beepr::beep()









## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)
col = `%IncMSE`

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()









#-------------------------- T1.12 : RF : bloo50 with coordinates + year + sand + rock  ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1] "northness"                  "eastness"                   "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                 <== CHANGE


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2756541    0.5250277     0.5166787 12.29773 14.89711  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))



col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 


# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.12.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.11.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.12_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./output/models/rf/T1.12_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./output/models/rf/T1.12_variable_importance_IncNodePurity.png")






#-------------------------- T1.13 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
# =================================== BEST RESULTS ON 637 sites =====================================
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"               <== CHANGE
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.286492    0.5352495     0.5244088 12.19134 14.79345  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.13.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.13.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.13_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./output/models/rf/T1.13_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./output/models/rf/T1.13²_variable_importance_IncNodePurity.png")







#-------------------------- T1.14 : RF : predict on surface transect samples only ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of surface_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4912711    0.7009073     0.7085043 11.69162 14.4469  NA            R    RF bloo50         NA


## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "surface_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.14.rds") 


## Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.14_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./output/models/rf/T1.14_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./output/models/rf/T1.14_variable_importance_IncNodePurity.png")





#-------------------------- T1.14.1 : RF : T1.14 again for reproductibily ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of surface_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 


# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4773418    0.6908993      0.701779 11.78027 14.55828  NA            R    RF bloo50         NA


## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "surface_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.14.1.rds") 


## Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.14.1.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.14.1_performance.rds")

perf
beepr::beep()




















#-------------------------- T1.15 : RF : predict on seabed transect samples only ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2946548    0.5428211     0.5062332 10.78129 13.12491  NA            R    RF bloo50         NA


## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.15.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.15_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.15_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.15_variable_importance_IncNodePurity.png")













#-------------------------- T1.15.1 : RF : T1.15 again to check reproductibility---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3041681    0.5515144     0.5185465 10.75682 13.07832  NA            R    RF bloo50         NA

## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.15.1.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.15.1_performance.rds")

perf
beepr::beep()






















#-------------------------- T1.16 : RF : predict on seabed+surface transect samples only ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect + seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2950991    0.5432302      0.541902 11.97284 14.53465  NA            R    RF bloo50         NA
# > 


## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect" | method == "surface_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))

table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.16.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.16_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.16_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.16_variable_importance_IncNodePurity.png")














#-------------------------- T1.17 : RF : predict on seabed+surface transect samples + add method variable ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect +seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect + seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr", "method" <== CHANGE            


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 



## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect" | method == "surface_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))

table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "method")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

unique(predictors$method)


# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.17.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.17_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.17_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.17_variable_importance_IncNodePurity.png")


















#-------------------------- T1.18 : RF : predict on seabed transect samples only, 10 RESPONSE VARIABLES ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect only
# [1] "R"          "Crypto"     "Elasmo"     "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   
# [10] "AngelShark"

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# [1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr         MAE       RMSE AIC Response_Var Model     CV Train_Size
# 1  0.302541819   0.55003802     0.5168938 10.78855561 13.0876425  NA            R    RF bloo50         NA
# 2  0.158950880   0.39868644     0.4534418  3.16379972  3.7601726  NA       Crypto    RF bloo50         NA
# 3  0.029384408   0.17141881     0.1657163  1.21261059  1.5444501  NA       Elasmo    RF bloo50         NA
# 4  0.330589200   0.57496887     0.6030616  0.16297994  0.2608343  NA        DeBRa    RF bloo50         NA
# 5  0.081096820   0.28477503     0.3142297  1.14945787  1.4239121  NA      RedList    RF bloo50         NA
# 6  0.306150507   0.55330869     0.5254283  7.95842153  9.9465155  NA         LRFI    RF bloo50         NA
# 7  0.135010100   0.36743721     0.3270401  2.43256574  2.9766539  NA      TopPred    RF bloo50         NA
# 8  0.320740345   0.56633943     0.4933620  5.75763920  7.1039701  NA   Commercial    RF bloo50         NA
# 9  0.069319536   0.26328603     0.2461465  0.34421391  0.4068560  NA      Grouper    RF bloo50         NA
# 10 0.002987102  -0.05465439    -0.1313480  0.03393542  0.1268795  NA   AngelShark    RF bloo50         NA
# 

## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select(-replicates)
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Fit model ----------
# Initialize list to store models
models_list <- list()

# Iterate over each response variable
for (response_var in colnames(ind)) {   # All resp var
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.18.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.18.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.18_performance.rds")

perf
beepr::beep()


















## Variables importance ----

md     <- models_list
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.18_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.18_variable_importance_IncNodePurity.png")














#-------------------------- T1.19 : RF : predict on surface transect samples only, 10 RESPONSE VARIABLES ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of surface_transect only
# [1] "R"          "Crypto"     "Elasmo"     "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   
# [10] "AngelShark"

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of surface_transect only
# [1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"

#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr        MAE       RMSE AIC Response_Var Model     CV Train_Size
# 1  0.484215472   0.69585593     0.7025530 11.7484758 14.5019672  NA            R    RF bloo50         NA
# 2  0.525075892   0.72462121     0.7337077  3.8385011  4.7280394  NA       Crypto    RF bloo50         NA
# 3  0.002011485   0.04484958     0.1289458  0.8112341  1.0639324  NA       Elasmo    RF bloo50         NA
# 4  0.258669757   0.50859587     0.5070156  0.2665504  0.4397418  NA        DeBRa    RF bloo50         NA
# 5  0.165926003   0.40734016     0.4297370  1.2409074  1.5719045  NA      RedList    RF bloo50         NA
# 6  0.401444344   0.63359636     0.6365691  8.5981903 10.6499926  NA         LRFI    RF bloo50         NA
# 7  0.186106889   0.43140108     0.4715836  2.6195351  3.2732189  NA      TopPred    RF bloo50         NA
# 8  0.312576732   0.55908562     0.5699395  6.2202034  7.7862518  NA   Commercial    RF bloo50         NA
# 9  0.149348948   0.38645692     0.4047496  0.2956961  0.3739465  NA      Grouper    RF bloo50         NA
# 10          NA           NA            NA  0.0000000  0.0000000  NA   AngelShark    RF bloo50         NA


## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "surface_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select(-replicates)
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Fit model ----------
# Initialize list to store models
models_list <- list()

# Iterate over each response variable
for (response_var in colnames(ind)) {   # All resp var
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.19.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.19.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.19_performance.rds")

perf
beepr::beep()






















#-------------------------- T1.20: RF : bloo50 with coordinates + year + sand + rock without eastness & northness + removing site with <5spp. --------------------- 
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 with sites with R<5 removed + removing eastness & northness
# [1]               "tpi_mean_log"               <== CHANGE
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2702516    0.5198573       0.51066 11.92936 14.36478  NA            R    RF bloo50         NA


## Prep data for model ----
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")


# Remove sites with R <5
tot <- tot %>%
  filter(R >=5)

# Remove sites with R <5 in pred and div
pred <- pred %>% filter(replicates %in% tot$replicates)
div <- div %>%  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.20.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.20.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.20_performance.rds")

perf
beepr::beep()






























#-------------------------- T1.20.1: RF : T1.20 to check --------------------- 
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 with sites with R<5 removed + removing eastness & northness
# [1]               "tpi_mean_log"               <== CHANGE
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 


# R2 Pearson_Corr Spearman_Corr     MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.268663    0.5183271     0.5112168 11.9171 14.37421  NA            R    RF bloo50         NA


## Prep data for model ----
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")


# Remove sites with R <5
tot <- tot %>%
  filter(R >=5)

# Remove sites with R <5 in pred and div
pred <- pred %>% filter(replicates %in% tot$replicates)
div <- div %>%  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.20.1.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.20.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.20.1_performance.rds")
perf <- readRDS("./output/models/rf/T1.20.1_performance.rds")
perf
beepr::beep()
































#-------------------------- T1.21 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 3 -------------

## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 3

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"               <== CHANGE
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2830828    0.5320553     0.5232878 12.18647 14.8092  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.21.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.13.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.21_performance.rds")

perf
beepr::beep()




















#-------------------------- T1.22 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2925358    0.5408657     0.5311549 12.13962 14.73589  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.22.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.13.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.22_performance.rds")

perf
beepr::beep()





















#-------------------------- T1.23 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 1, ntree = 500, nodesize = 1 ---------------

## Description ----
#--- Model : RF
# Parameters : mtry = 1, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2879643     0.536623      0.527363 12.29391 14.89635  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.23.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.13.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.23_performance.rds")

perf
beepr::beep()






















#-------------------------- T1.24 : RF : predict on seabed transect samples only ---------------------
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2995532    0.5473146     0.5372111 12.05779 14.64999  NA            R    RF bloo50         NA

## Prep data for model ----

# Keep only surface_transect 
tot <- tot %>%
  filter(method == "seabed_transect")

pred <- pred %>% filter(replicates %in% tot$replicates)

div <- div %>%
  filter(replicates %in% tot$replicates)

# Predictors columns
# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)

tot <- tot %>%
  mutate(year = year(date))
table(tot$year)


col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")



# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)




# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)









## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.24.rds") 


## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.24_performance.rds")

perf
beepr::beep()
























#-------------------------- T1.25 : RF (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add temp_min + area_km2 ----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3034782    0.5508885     0.5428111 12.00968 14.61003  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m", "area_km2")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m", "area_km2")

# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.25.rds") 


## Evaluate performance ----
models_list <- readRDS("./output/models/rf/T1.25.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.25_performance.rds")

perf <- readRDS("./output/models/rf/T1.25_performance.rds")
perf
beepr::beep()
























#-------------------------- T1.26 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add temp_min + surface only ----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# 
# 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4605439    0.6786338     0.6819715 11.78017 14.60515  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m")

# Filter method=surface_transect
tot <- tot %>%
  filter(method == "surface_transect")

div <- div %>%
  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.26.rds") 




## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.26_performance.rds")
perf <- readRDS("./output/models/rf/T1.26_performance.rds")
perf
beepr::beep()
























#-------------------------- T1.27 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- surface only ----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4738402    0.6883605     0.6879786 11.65696 14.48351  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m")

# Filter method=surface_transect
tot <- tot %>%
  filter(method == "surface_transect")

div <- div %>%
  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.27.rds") 




## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.27_performance.rds")

perf
beepr::beep()
























#-------------------------- T1.28 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 5)---------------------------------------------
#-------------------------- surface only | no x y----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 


# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4692602    0.6850257     0.6845877 11.72913 14.54217  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m")

# Filter method=surface_transect
tot <- tot %>%
  filter(method == "surface_transect")

div <- div %>%
  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.28.rds") 




## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.28_performance.rds")

perf
beepr::beep()
























#-------------------------- Compare T1.28 with T1.14 -----------------------------
t14 <- readRDS("./output/models/rf/T1.14.rds")
t28 <- readRDS("./output/models/rf/T1.28.rds")



# Check if same predictors
all.equal(t14$R$RF$bloo50$Fold_1$model$terms, t28$R$RF$bloo50$Fold_1$model$terms)

# Subset Classes 'formula'
t14$R$RF$bloo50$Fold_1$model$terms
t28$R$RF$bloo50$Fold_1$model$terms


# Print differences
setdiff(attr(t14$R$RF$bloo50$Fold_1$model$terms, "variables") %>% as.character(), 
        attr(t28$R$RF$bloo50$Fold_1$model$terms, "variables") %>% as.character())





#-------------------------- T1.29 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 5)---------------------------------------------
#-------------------------- surface only | x and y ----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 



# R2 Pearson_Corr Spearman_Corr      MAE    RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4789338    0.6920504     0.7015957 11.90305 14.6627  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m")

# Filter method=surface_transect
tot <- tot %>%
  filter(method == "surface_transect")

div <- div %>%
  filter(replicates %in% tot$replicates)


# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.29.rds") 




## Evaluate performance ----
#models_list <- readRDS("./output/models/rf/T1.14.rds")
perf <- evaluate_models(models_list)
# save
saveRDS(perf, "./output/models/rf/T1.29_performance.rds")

perf
beepr::beep()

























##--------------------------  T1.29 = T1.28 -----------------------------
#--------------------------- T1.30 : RF (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add temp_min  ----------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv + selection of seabed_transect only
# R

#---- Predictors : 
# from predictors_sel_v.1.4 + selection of seabed_transect only
# 1] "x"                          "y"                          "tpi_mean_log"               "port_dist_m_weight"        
# [5] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# [9] "sal_mean_1m"                "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [17] "year"                       "soft_bottom_clr"            "rock_clr"                 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.2615015    0.5113721     0.4957995 12.31004 14.94447  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "temp_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "temp_min_1m")

# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")






## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Save model ----
saveRDS(models_list, "./output/models/rf/T1.30.rds") 


## Evaluate performance ----
perf <- evaluate_models(models_list)
perf
beepr::beep()
# save
saveRDS(perf, "./output/models/rf/T1.30_performance.rds")

































#-------------------------- T1.31 : RF : bloo20 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 20 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr     MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3588587    0.5990481     0.5931061 11.1705 13.96587  NA            R    RF bloo20         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 20

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo20"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf <- readRDS("./output/models/rf/T1.31_performance.rds")
perf





# Save ----
saveRDS(perf, "./output/models/rf/T1.31_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.31.rds") 

#-------------------------- T1.32 : RF : bloo10 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 
# from predictors_sel_v.1.3=4 
# [1]               "tpi_mean_log"              
# [4] "port_dist_m_weight"         "bathy_mean"                 "wind_mean_1m"              
# [7] "vel_mean_1m"                "temp_mean_1m"               "sal_mean_1m"               
# [10] "grouped_nb_habitat_per_km2" "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"     
# [13] "shore_dist_m_weight_log"    "gravity_mean_log"           "cop_chl_month_mean_log"    
# [16] "Boat_density_month_log"     "x", "y" , "year"                


#---- CV config :
# BLOO CV with buffer size = 20 km

#---- Perf : 

# R2 Pearson_Corr Spearman_Corr     MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3588587    0.5990481     0.5931061 11.1705 13.96587  NA            R    RF bloo20         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
names(pred1.5)

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

names(tot)


tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")



## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 10

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo10"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.32_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.32.rds") 

## Plot observed vs predicted ----
md <- models_list
evaluate_models(md, output = "all")

## Variance importance plot ----
folds  <- md$R$RF$bloo10

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.32_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.32_variable_importance_IncNodePurity.png")




































#-------------------------- T1.33 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

#-------------------------- + area_km2 + estimated_volume_total ------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3969026     0.630002     0.6275899 11.19965 13.80481  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "area_km2")), by = "replicates")

sort(names(tot))

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "area_km2", "estimated_volume_total")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Fit model ----------

# Iterate over each response variable
for (response_var in colnames(ind)) {   # Just R
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.33_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.33.rds") 


## Plot observed vs predicted ----
md <- readRDS("./output/models/rf/T1.33.rds")
#md <- models_list
source("./utils/Fct_Modeling.R")
evaluate_models(md, output = "all")

## Variance importance plot ----
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    tibble::rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.33_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.33_variable_importance_IncNodePurity.png")






































#-------------------------- T1.34 : RF : T33 on RedList & Elasmo ---------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv


#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 500 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.06932629    0.2632989     0.2657861 1.111273 1.439253  NA       Elasmo    RF bloo50         NA
# 2 0.10300356    0.3209417     0.3291147 1.220040 1.524278  NA      RedList    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "area_km2")), by = "replicates")

sort(names(tot))

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "area_km2", "estimated_volume_total")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select(c("Elasmo", "RedList"))
names(ind)

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.34_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.34.rds") 





#-------------------------- T1.35 : RF : bloo50 with coordinates + year + sand + rock without eastness & northness --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

#-------------------------- + area_km2scaled(estimated_volume_total) +scale(x,y,year)  ------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 500 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.3879857     0.622885     0.6261796 11.06587 13.71551  NA            R    RF bloo50         NA

## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "area_km2")), by = "replicates")

sort(names(tot))

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "area_km2", "estimated_volume_total")

tot$estimated_volume_total <- scale(tot$estimated_volume_total)
tot$x <- scale(tot$x)
tot$y <- scale(tot$y)
tot$year <- scale(tot$year)

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


## Fit model ----------
# Initiate models_lists
models_list <- list()


# Iterate over each response variable
for (response_var in colnames(ind)) {   # Just R
  
  # Fit the model using fit_models function
  model <- fit_models(
    dataset = tot, 
    response_var = response_var, 
    predictors = colnames(predictors), 
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(models_list, "./output/models/rf/T1.35.rds") 
saveRDS(perf, "./output/models/rf/T1.35_performance.rds")






#-------------------------- T1.36 : All_fitted_models bloo50/bloo20 (All_fitted_20) [MARBEC GPU GAETAN] ---------------------
md <- readRDS("./output/models/rf/All_fitted_models.rds")

## Description ----
# Predictors 
md$R$RF$bloo50$Fold_1$model$forest$independent.variable.names

# [1] "x"                          "y"                          "northness"                 
# [4] "eastness"                   "tpi_mean_log"               "port_dist_m_weight"        
# [7] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [10] "temp_mean_1m"               "sal_mean_1m"                "grouped_nb_habitat_per_km2"
# [13] "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"      "shore_dist_m_weight_log"   
# [16] "gravity_mean_log"           "cop_chl_month_mean_log"     "Boat_density_month_log"    
# [19] "year"                       "soft_bottom_clr"            "rock_clr"                  
# [22] "area_km2"                   "temp_min_1m"      

# Resp var


# Model + parameters
md$R$RF$bloo50$Fold_1$model$num.trees # ntree = 300
md$R$RF$bloo50$Fold_1$model$mtry # mtry = 3
md$R$RF$bloo50$Fold_1$model$min.node.size # min.node.size = 5
md$R$RF$bloo50$Fold_1$model$r.squared # 0.4464482


## Evaluate perf ----

evaluate_models <- function(models_list,
                            response_var = NULL,
                            models       = NULL,
                            cv_names     = NULL,   # names of CV methods to use (e.g. "env_k6", "bloo50")
                            cv_splits    = NULL,   # optional: list of CV splits to compute train sizes
                            output       = c("summary", "fold", "plot", "all")) {
  
  output <- match.arg(output)
  cat("\n🔹 Starting Model Evaluation\n")
  
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
      
      cat("\n🔹 Evaluating model:", model_name, "for response variable:", resp_var, "\n")
      
      available_cvs <- names(models_list[[resp_var]][[model_name]])
      
      cv_methods <- if (!is.null(cv_names)) {
        intersect(cv_names, available_cvs)
      } else {
        available_cvs
      }
      
      for (cv_name in cv_methods) {
        
        cat("\n🔹 Processing CV Method:", cv_name, "\n")
        
        model_folds <- models_list[[resp_var]][[model_name]][[cv_name]]
        # Optional CV splits to get train sizes
        this_cv_split <- if (!is.null(cv_splits)) cv_splits[[cv_name]] else NULL
        
        all_observed    <- c()
        all_predictions <- c()
        aic_values      <- c()
        train_sizes     <- c()
        
        for (i in seq_along(model_folds)) {
          
          cat("   ➡️ Fold:", i, "/", length(model_folds), "\n")
          
          fold_name <- paste0("Fold_", i)
          fold_data <- model_folds[[fold_name]]
          
          
          test_set <- fold_data$test_set
          
          
          test_df <- as.data.frame(st_drop_geometry(test_set))
          
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
              fold_pred <- fold_data$predictions
            } else if (model_name == "GLSRF") {
              fold_pred <- fold_data$predictions
              
            } else {
              predict(fold_data$model, newdata = test_df, type = "response")
            }
          }, error = function(e) {
            cat("      ❌ Prediction failed for", model_name, "on", fold_name, ":", e$message, "\n")
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
          if (!is.null(this_cv_split) && length(this_cv_split) >= i && !is.null(this_cv_split[[i]]$train)) {
            this_train_size <- nrow(this_cv_split[[i]]$train)
            train_sizes     <- c(train_sizes, this_train_size)
          }
          
          #-----------------------------
          # Per-fold metrics
          #-----------------------------
          if (!all(is.na(fold_pred))) {
            fold_metrics <- compute_metrics(fold_observed, fold_pred,
                                            if (model_name %in% c("GLM", "GAM", "GLM.nb", "spaMM"))
                                              AIC(fold_data$model) else NULL)
            
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
            cat("      ⚠️ All predictions NA for", fold_name, "- no fold metrics.\n")
          }
        } # end folds
        
        
        #-----------------------------
        # Global (across folds) metrics
        #-----------------------------
        if (length(all_observed) == 0 || length(all_predictions) == 0) {
          cat("\n⚠️ No valid predictions for model", model_name, "and CV", cv_name, "- skipping global metrics.\n")
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
              
              p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Observed, y = Predicted)) +
                ggplot2::geom_point(alpha = 0.6) +
                ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
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
        }
      }
    }
  }
  
  #-----------------------------
  # Return according to 'output'
  #-----------------------------
  if (output == "summary") {
    return(performance_summary)
  } else if (output == "fold") {
    return(performance_per_fold)
  } else if (output == "plot") {
    return(plot_list)
  } else {  # "all
    return(list(
      summary  = performance_summary,
      per_fold = performance_per_fold,
      plots    = plot_list
    ))
  }
}

perf <- evaluate_models(md)
perf
# R2 Pearson_Corr Spearman_Corr         MAE       RMSE AIC Response_Var Model     CV
# 1  0.308884015   0.55577335    0.54931249 11.95673879 14.5533653  NA            R    RF bloo50
# 2  0.318806201   0.56462926    0.56204326  3.55933402  4.3646251  NA       Crypto    RF bloo50
# 3  0.050708305   0.22518505    0.22811567  1.14256753  1.4645368  NA       Elasmo    RF bloo50
# 4  0.302181080   0.54971000    0.55333411  0.18424532  0.3145530  NA        DeBRa    RF bloo50
# 5  0.061961983   0.24892164    0.25941153  1.25033416  1.5577066  NA      RedList    RF bloo50
# 6  0.293423673   0.54168595    0.53316755  8.70881603 10.7038085  NA         LRFI    RF bloo50
# 7  0.115327071   0.33959840    0.34073493  2.64682222  3.2316517  NA      TopPred    RF bloo50
# 8  0.278197745   0.52744454    0.51548822  6.22832108  7.6185525  NA   Commercial    RF bloo50
# 9  0.080532690   0.28378282    0.28926803  0.36293335  0.4217820  NA      Grouper    RF bloo50
# 10 0.000138836  -0.01178287    0.02246302  0.02857265  0.1114742  NA   AngelShark    RF bloo50






#==== Surface only : T1.29/T1.28 : pseudo R2 = 0.68 =============================================================
#====  All data : T33 : pseudo R2 = 0.63   (CV50 with volume)         =============================================================#================================================================================================================
#-------------------------- T1.38 = T1.33 on all variables - bloo50 [MARBEC GPU MARIEKE] -------------------------------------------
## Evaluate performance ----
models_list <- readRDS("~/Downloads/T1.37.rds")
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.37_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.32.rds") 

## Plot observed vs predicted ----
md <- models_list
evaluate_models(md, output = "all")

## Variance importance plot ----
folds  <- md$R$RF$bloo10

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.32_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.32_variable_importance_IncNodePurity.png")






































#-------------------------- T1.39 : RF = T1.33 --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

#-------------------------- + area_km2 + estimated_volume_total ------------------------
#-------------------------- + depth_sampling ------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4014184    0.6335759     0.6286837 11.16216 13.75386  NA            R    RF bloo50         NA


## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "area_km2")), by = "replicates")

sort(names(tot))

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "area_km2", "estimated_volume_total", "depth_sampling")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.39_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.39.rds") 


## Plot observed vs predicted ----
md <- readRDS("./output/models/rf/T1.39.rds")
#md <- models_list
source("./utils/Fct_Modeling.R")
evaluate_models(md, output = "all")

## Variance importance plot ----
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    tibble::rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.39_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.39_variable_importance_IncNodePurity.png")














































#-------------------------- T1.40 : RF = T1.33 --------------------- 
#-------------------------- mtry = 2, ntree = 500, nodesize = 1 ---------------

#-------------------------- + area_km2 + estimated_volume_total ------------------------
#-------------------------- + depth_sampling + region ------------------------
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 500 km

#---- Perf :
## Description ----
#--- Model : RF
# Parameters : mtry = 2, ntree = 500, nodesize = 1

#--- Response var : 
# from div_indices_v1.0_sel_v1.1.csv
# R

#---- Predictors : 


#---- CV config :
# BLOO CV with buffer size = 50 km

#---- Perf : 
# R2 Pearson_Corr Spearman_Corr      MAE     RMSE AIC Response_Var Model     CV Train_Size
# 1 0.4014184    0.6335759     0.6286837 11.16216 13.75386  NA            R    RF bloo50         NA




## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "area_km2")), by = "replicates")

sort(names(tot))

# Add year
tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr", "area_km2", "estimated_volume_total", "depth_sampling", "region")

# Extract predictors
predictors <- tot %>% 
  dplyr::select(all_of(col_pred)) 
names(predictors)

# Extract response variables
ind <- div  %>% dplyr::select("R")
names(ind)

## CV config -----
# Convert tot to sf object
tot_sf_4326 <- st_as_sf(tot, coords = c("x", "y"), crs = 4326, remove = FALSE) # make sure you keep "x" and "y" columns with remove = FALSE (the coordinates columns will be needed for spaMM)


# Set coords to numeric
tot_sf_4326$x <- as.numeric(tot_sf_4326$x)
tot_sf_4326$y <- as.numeric(tot_sf_4326$y)




# LOO CV
# Load 
# loo <- loo_cv(tot_sf_4326)


# BLOO CV
# Compute different BLOO-CV configurations 
# Define buffer sizes
#buffer_sizes <- seq(from = 10, to = 200, by = 10)
buffer_sizes = 50

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
# cv_configs <- c(list(loo = loo), bloo_results_list)
cv_configs <- c(bloo_results_list)


rm(tot_sf_4326, loo , bloo_results_list, bloo_results, buffer_sizes, buffer_size)


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
    cv_configs = cv_configs["bloo50"],
    models = c("RF"), 
    distribution = NULL
  )
  
  # Store model in list
  models_list[[response_var]] <- model
}

## Evaluate performance ----
perf <- evaluate_models(models_list)
perf 
beepr::beep()





# Save ----
saveRDS(perf, "./output/models/rf/T1.40_performance.rds")

saveRDS(models_list, "./output/models/rf/T1.40.rds") 


## Plot observed vs predicted ----
md <- readRDS("./output/models/rf/T1.40.rds")
#md <- models_list
source("./utils/Fct_Modeling.R")
evaluate_models(md, output = "all")

## Variance importance plot ----
folds  <- md$R$RF$bloo50

# 1. Extract importance from each RF fold
importance_list <- lapply(folds, function(fold) {
  imp <- randomForest::importance(fold$model)
  as.data.frame(imp) |>
    tibble::rownames_to_column(var = "Feature")
})

# 2. Combine and average
all_importance <- bind_rows(importance_list, .id = "Fold")

# Check columns if you want
# names(all_importance)

mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(desc(IncNodePurity))   # or desc(`%IncMSE`) depending on your RF type



# 3. Plot
# 3. Plot (choose the column IncNodePurity or `%IncMSE`)

ggplot(mean_importance,
       aes(x = reorder(Feature, `%IncMSE`), y = `%IncMSE`)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave("./figures/models/rf/T1.40_variable_importance_percentIncMSE.png")

ggplot(mean_importance,
       aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_col() +
  coord_flip() +
  theme_minimal()
ggsave("./figures/models/rf/T1.40_variable_importance_IncNodePurity.png")









































