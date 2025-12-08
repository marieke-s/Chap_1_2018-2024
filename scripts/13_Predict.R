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
library(SHAPforxgboost)
library(spdep)
library(superml)
library(svglite)
library(tidyr)
library(units)
library(virtualspecies)
library(xgboost)




# Load functions
source("./utils/Fct_Modeling.R")

#===============================================================================
#------------- Load and prep data ------------------
# model T.0.0 ----
models_list <- readRDS("./output/models/T0.0.rds")

# predictors_sel_v1.3
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.3.gpkg")



#------------- Plots ------------
md <- models_list$R$XGB$bloo50$Fold_1$model
# plot only the first tree and display the node ID:
xgb.plot.tree(model = md, trees = 0, show_node_id = TRUE)

p <- xgb.plot.multi.trees(model = md, features_keep = 3)
print(p)

#-------------------------- T0.0 : Predict surface on last fold -------------

md <- models_list$R$XGB$bloo50$Fold_637$model

# Load predictors
pred <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors_sel_v1.3.gpkg") %>%
  dplyr::select(-c("x", "y", , "grouped_main_habitat", "sal_mean_40m_month", "temp_mean_40m_month")) %>% # remove 40m 
  rename(wind_mean_1m = ws_mean_surface_month, 
         vel_mean_1m = vel_mean_surface_month, 
         temp_mean_1m = temp_mean_surface_month, 
         sal_mean_1m = sal_mean_surface_month) 
                  


setdiff(names(pred), md$feature_names)
setdiff(md$feature_names, names(pred))

test <- pred

# 1. Drop geometry and keep only the model features, in the right order
X_pred <- test %>%
  sf::st_drop_geometry() %>%
  dplyr::select(all_of(md$feature_names)) %>%
  as.matrix()

# Optional: sanity check
setdiff(md$feature_names, colnames(X_pred))  # should be character(0)

# 2. Predict with your xgboost model
y_hat <- predict(md, newdata = X_pred)

# 3. Attach predictions back to the sf object
test$prediction <- y_hat
sapply(test, function(x) sum(is.na(x)))
summary(test$prediction)


# 4. Plot 
test <- st_read("./output/predictions/T0.0/T0.0_R_XGB_bloo50_Fold637_predictions_2023-07-01_surface.gpkg")
ggplot(test) +
  geom_sf(aes(fill = prediction), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Predicted Species Richness at 1m depth in June 2023",
       fill = "Prediction") +
  theme_minimal()

# Export predictions----
st_write(test, "./output/predictions/T0.0/T0.0_R_XGB_bloo50_Fold637_predictions_2023-07-01_surface.gpkg", delete_dsn = TRUE)

#-------------------------- T0.0 : Predict 40m on last fold -------------

md <- models_list$R$XGB$bloo50$Fold_637$model

# Load predictors
pred <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors_sel_v1.3.gpkg") %>%
  dplyr::select(-c("x", "y", , "grouped_main_habitat", "sal_mean_surface_month", "temp_mean_surface_month")) %>%  
  rename(wind_mean_1m = ws_mean_surface_month, 
         vel_mean_1m = vel_mean_surface_month, 
         temp_mean_1m = temp_mean_40m_month, 
         sal_mean_1m = sal_mean_40m_month) 



setdiff(names(pred), md$feature_names)
setdiff(md$feature_names, names(pred))

test <- pred

# 1. Drop geometry and keep only the model features, in the right order
X_pred <- test %>%
  sf::st_drop_geometry() %>%
  dplyr::select(all_of(md$feature_names)) %>%
  as.matrix()

# Optional: sanity check
setdiff(md$feature_names, colnames(X_pred))  # should be character(0)

# 2. Predict with your xgboost model
y_hat <- predict(md, newdata = X_pred)

# 3. Attach predictions back to the sf object
test$prediction <- y_hat
sapply(test, function(x) sum(is.na(x)))
summary(test$prediction)


# 4. Plot 
test <- st_read("./output/predictions/T0.0/T0.0_R_XGB_bloo50_Fold637_predictions_2023-07-01_40m.gpkg")
ggplot(test) +
  geom_sf(aes(fill = prediction), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Predicted Species Richness at 40m depth in June 2023",
       fill = "Prediction") +
  theme_minimal()

# Export predictions ----
st_write(test, "./output/predictions/T0.0/T0.0_R_XGB_bloo50_Fold637_predictions_2023-07-01_40m.gpkg", delete_dsn = TRUE)


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


#-------------------------- T1.1 : Retrain & Predict on 04-09/2023 ------------------
#--- Retrain model T1.1 : XGBOOST --------------



## Description T1.1 ----
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

# Prep data ----
predictors <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()

# Extract response variables
ind <- div  %>% dplyr::select(-c("replicates"))

# Retrain (all response var) -----
models <- lapply(colnames(ind)[-1], function(response_var) {
  
  print(paste0("==== Training model for response variable: ", response_var))
  
  # Convert character columns to factors
  train_df <- tot %>% mutate(across(where(is.character), as.factor))
  
  # Drop unused levels
  dt <- droplevels(train_df[, names(predictors)])

  # Build Matrix  
  train_x <- model.matrix(~ . - 1, data = dt)
  train_matrix <- xgb.DMatrix(data = train_x, label = train_df[[response_var]])
  
  # Train
  xgb.train(
    params = list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.1, max_depth = 6),
    data = train_matrix,
    nrounds = 100,
    verbose = 0 )
  
    })

names(models) <- colnames(ind)[-1]


# Save models ----
saveRDS(models, "./output/models/xgboost/T1.1_retrain.rds")

# Predict ----
# --- Prep predict data ---
# Load all prediction grids 
files <- list.files("./data/processed_data/predictors/Prediction_grid_v1.1/v1.3/", pattern = "*.gpkg", full.names = TRUE)
pred_list <- lapply(files, st_read)

# Select appropriate predictors 
vars <- c(
  "northness",
  "eastness",
  "tpi_mean_log",
  "port_dist_m_weight",
  "bathy_mean",
  "wind_mean_1m",
  "vel_mean_1m",
  "temp_mean_1m",
  "sal_mean_1m",
  "grouped_nb_habitat_per_km2",
  "canyon_dist_m_weight_log",
  "mpa_dist_m_weight_log",
  "shore_dist_m_weight_log",
  "gravity_mean_log",
  "cop_chl_month_mean_log", # NOT LOGGED IN GRIDS ! 
  "Boat_density_month_log"
)


# Rename cols : boat_density_month_log = Boat_density_month_log


setdiff(vars, names(pred_list[[3]]))



names(pred_list[[3]])

# Rename predictors 










#===============================================================================
#-------------------------- T.1.33:RF: Load and prep data ------------------
#-------------------------- Load and prep data ------------------
# predictors_sel_v.1.4 
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

# mtdt_7_sel_v1.1
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

# div_indices_sel_v1.1.gpkg
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# tot
tot <- pred %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(mtdt), by = "replicates") %>%
  left_join(div, by = "replicates")

# --- Prep data for model ---

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

# Train df
train_df <- tot %>% dplyr::select(c(names(ind), names(predictors)))
names(train_df)

#-------------------------- Retrain model on R ---------------------------
formula <- as.formula(paste("R", "~", paste(names(predictors), collapse = " + ")))

rf <- randomForest::randomForest(formula = formula, data = train_df, ntree = 500, importance = TRUE, mtry = 2, nodesize = 1)

# Save rf
saveRDS(rf, "./output/models/rf/T1.33-retrain.rds")

#-------------------------- Load retrained model ---------------
md <- readRDS("./output/models/rf/T1.33-retrain.rds")
#-------------------------- Predict ---------------------------
data <- readRDS("./data/processed_data/predictors/Prediction_grid_v1.1/Export_grid_v1.1_2023-Apr-Sept_with_predictors-tr_v1.4.rds")



# Bottom ----
# Remove surface sal temp cols
data <- lapply(data, \(df) df |> dplyr::select(!matches("^(temp|sal).*_surf$")))

# Drop geometry 
data <- lapply(data, \(df) df %>% st_drop_geometry())

# Put all colnames to lowercase
data <- lapply(data, function(df) {
  names(df) <- tolower(names(df))
  df
})

# Replace '_surf' or '_bottom' by '_1m' in colnames
data <- lapply(data, function(df) {
  names(df) <- gsub("_(surf|bottom)$", "_1m", names(df))
  df
})




sort(setdiff(names(predictors), names(data$`2023-05-01`)))


sort(names(data$`2023-05-01`))




















