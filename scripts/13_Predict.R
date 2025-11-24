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

#----------- T0.0 : Predict surface on last fold -------------

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

#----------- T0.0 : Predict 40m on last fold -------------

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

