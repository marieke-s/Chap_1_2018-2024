

#---------------- PRELIMINAR PREDICTORS SELECTION / CLEANING ----------------

# Remove range 
# Remove wind stress mins

# Remove range min, max, range, of spatial predictors ----
# Explanation : We remove range min and max predictors that were meant to represent variability within the buffer for spatial predictors (i.e. non-climatic pred).

# Remove range variable
range_vars <- colnames(pred)[grep("_range", colnames(pred))]
pred <- pred %>% dplyr::select(-range_vars)
rm(range_vars)

# Remove min max variable
pred <- pred %>% dplyr::select(-c("dist_port_min", "dist_port_max", "gravity_min", "gravity_max"))

# -7 predictors











# Remove zero variance predictors ----
caret::nearZeroVar(pred_raw %>% st_drop_geometry(), saveMetrics = TRUE)

# nvz : canyon_type, habitats_artificiels_mean, rock_mean + win_min 7day 1 month and 1year
# zeroVar :  win_min 7day 1 month and 1year

#--------------- PREDICTORS SELECTION PROCEDURE ----------------