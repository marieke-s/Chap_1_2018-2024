



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









