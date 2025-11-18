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

# Save model ----
saveRDS(models_list, "./output/models/T0.0.rds") # Response var = R, Crypto and Elasmo

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
















