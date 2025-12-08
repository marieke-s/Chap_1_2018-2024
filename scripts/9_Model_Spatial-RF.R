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
library(sf)
library(spatialRF)
library(units)
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
library(spatialRF)
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
#-------------  T.1 : Load and prep data ------------------
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



#------------- Prep data for spatial-RF ------------------------------

tot_sf <- st_as_sf(tot, coords = c("x", "y"), crs = 4326)  # adapt names if needed
tot_m  <- st_transform(tot_sf, 3857)                       # metric CRS (meters)

xy_m <- st_coordinates(tot_m) |> as.data.frame()
colnames(xy_m) <- c("x", "y")

nrow(tot)    # 637
nrow(xy_m)   # must be 637
head(xy_m)

dm  <- st_distance(tot_m)             # units [m]
distance_matrix <- as.matrix(drop_units(dm))  # plain numeric

dim(distance_matrix)   # 637 637
range(distance_matrix)

distance_thresholds <- c(0, 1000, 2000, 5000, 10000, 50000, 100000)


#===============================================================================

#-------------------------- T.0.0 : ---------------------------------------------

## --- I. Data Explo ----
# Moaran I heatmap ------
spatialRF::plot_training_df_moran(
  data                    = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix         = distance_matrix,
  distance.thresholds     = distance_thresholds,
  fill.color              = viridis::viridis(100, option = "F", direction = -1),
  point.color             = "gray40"
)

ggsave("./figures/models/spatialrf/T.0.0_moran_I_heatmap.png")

# Scatter plots of predictors vs R ------
spatialRF::plot_training_df(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  ncol = 3,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

ggsave("./figures/models/spatialrf/T.0.0_plot_R_vs_pred.png")


# Interactions between predictors -----
interactions <- spatialRF::the_feature_engineer(
  data                    = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  xy                      = xy_m,
  importance.threshold    = 0.50,
  cor.threshold           = 0.60,
  repetitions             = 100,
  seed                    = 1,
  verbose                 = TRUE
)

# Save results
saveRDS(interactions, "./output/models/spatialrf/T.0.0_interactions_R.rds")

kableExtra::kbl(
  head(interactions$screening, 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# Save kable as svg


## --- II. Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.0_nonspatial_R.rds")


# Model performance ----
spatialRF::print_performance(model.non.spatial)

# Model performance 
# - R squared (oob):                  0.4860573
# - R squared (cor(obs, pred)^2):     0.941395
# - Pseudo R squared (cor(obs, pred)):0.9702551
# - RMSE (oob):                       12.26903
# - RMSE:                             5.236
# - Normalized RMSE:                  0.1939259


## --- III. Spatial rf ----------------
rf_spatial_moran <- rf_spatial(
  data                    = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix         = distance_matrix,    # 637 x 637
  distance.thresholds     = distance_thresholds,
  method                  = "mem.moran.sequential",
  xy                      = xy_m,               # 637 x 2
  n.cores                 = 1
)

saveRDS(rf_spatial_moran, "./output/models/spatialrf/T.0.0_spatial_R.rds")

# Model performance ----
spatialRF::print_performance(rf_spatial_moran)
# Model performance 
# - R squared (oob):                  0.4860573
# - R squared (cor(obs, pred)^2):     0.941395
# - Pseudo R squared (cor(obs, pred)):0.9702551
# - RMSE (oob):                       12.26903
# - RMSE:                             5.236
# - Normalized RMSE:                  0.1939259


## --- IV. Compare models ----
comparison <- spatialRF::rf_compare(
  models = list(
    `Non-spatial` = model.non.spatial,
    `Spatial` = rf_spatial_moran
  ),
  xy = xy_m,
  repetitions = 30,
  training.fraction = 0.8,
  metrics = "r.squared",
  seed = 1)

ggsave("./figures/models/spatialrf/T.0.0_spatial_vs_non-spatial_comparison_R.png")

## --- V. Variable contributions to transferability & Importance -----
# Contribution to transferability ----
contrib <- spatialRF::rf_importance(
  model        = model.non.spatial,
  xy           = xy_m,
  repetitions  = 5,
  metric       = "r.squared",
  distance.step.x = 50000,    # in meters (same units as xy_m)
  distance.step.y = 50000,
  n.cores      = 1,
  verbose      = TRUE
)

# Check contrib object
contrib$performance
contrib$importance$cv.per.variable.plot
contrib$importance


# Save contrib object
saveRDS(contrib, "./output/models/spatialrf/T.0.0_non-spatial_transferability_contrib_R.rds")

# Save plot
ggsave("./figures/models/spatialrf/T.0.0_non-spatial_transferability_contrib_R.png")

# Response curves ----
spatialRF::plot_response_curves(
  model.non.spatial,, 
  quantiles = 0.5,
  ncol = 3
)

ggsave("./figures/models/spatialrf/T.0.0_non-spatial_response_curves_R.png")


spatialRF::plot_response_curves(
  model.non.spatial,
  quantiles = c(0.1, 0.5, 0.9),
  line.color = viridis::viridis(
    3, #same number of colors as quantiles
    option = "F", 
    end = 0.9
  ),
  ncol = 3,
  show.data = TRUE
)


ggsave("./figures/models/spatialrf/T.0.0_non-spatial_response_curves_with_data_R.png")

# Variable importance ----

spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)

ggsave("./figures/models/spatialrf/T.0.0_non-spatial_variable_importance_R.png")

# Contribution vs importance ----
contrib$importance$per.variable %>% 
  ggplot2::ggplot() +
  ggplot2::aes(
    x = contrib$importance$per.variable$importance.oob, # oob importance
    y = contrib$importance$per.variable$importance.cv # contribution to transferability
  ) + 
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::xlab("Importance (out-of-bag)") + 
  ggplot2::ylab("Contribution to transferability") + 
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red4")

ggsave("./figures/models/spatialrf/T.0.0_non-spatial_contrib_vs_importance_R.png")


## --- VI.Cross validation on non spatial model ----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  # repetitions = 30,         # number of spatial folds
  training.fraction = 0.90, #training data fraction on each fold
  metrics = "r.squared",
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
) 

 

# model.non.spatial_cv results ----
spatialRF::print_performance(model.non.spatial_cv)

# check cross validation results 
model.non.spatial_cv$evaluation

## --- AUTRE: Residuals interpratation ----

#-------------------------- T.0.1 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 1,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.1_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.1_nonspatial_cv_R.rds")

# 4  Testing pseudo.r.squared  0.5410000              0.1897728000  0.3192500  0.6562500  0.4851667 0.083 0.2037365128  0.2150000  0.6740000
# 5  Testing        r.squared  0.3010000              0.2164596000  0.1072500  0.4302500  0.2696667 0.076 0.1857360134  0.0460000  0.4540000
# 6  Testing             rmse 13.4245000              1.3232205000 12.9665000 14.4607500 13.5845000 0.429 1.0513438543 12.1860000 14.8520000

# model.non.spatial_cv results 
model.non.spatial_cv <- readRDS("./output/models/spatialrf/T.0.1_nonspatial_cv_R.rds")
model.non.spatial_cv$evaluation

# Plot folds
pr <- tot[, c("x", "y")]
pr$group.2 <- pr$group.1 <- "Training"
pr[model.non.spatial_cv$evaluation$spatial.folds[[3]]$testing, "group.1"] <- "Testing"
pr[model.non.spatial_cv$evaluation$spatial.folds[[6]]$testing, "group.2"] <- "Testing"

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "white") +
  ggplot2::geom_point(data = pr,
                      ggplot2::aes(
                        x = x,
                        y = y,
                        color = group.1
                      ),
                      size = 2
  ) +
  ggplot2::scale_color_viridis_d(
    direction = -1, 
    end = 0.5, 
    alpha = 0.8, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Group") +
  ggplot2::scale_x_continuous(limits = c(2.5, 9.5)) +
  ggplot2::scale_y_continuous(limits = c(41, 44))  +
  ggplot2::ggtitle("Spatial fold 1") + 
  ggplot2::theme(
    legend.position = "none", 
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "white") +
  ggplot2::geom_point(data = pr,
                      ggplot2::aes(
                        x = x,
                        y = y,
                        color = group.2
                      ),
                      size = 2
  ) +
  ggplot2::scale_color_viridis_d(
    direction = -1, 
    end = 0.5, 
    alpha = 0.8, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Group") +
  ggplot2::scale_x_continuous(limits = c(2.5, 9.5)) +
  ggplot2::scale_y_continuous(limits = c(41, 44)) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) + 
  ggplot2::ggtitle("Spatial fold 6") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("")

p1 | p2


spatialRF::plot_evaluation(model.non.spatial_cv)

#-------------------------- T.0.2 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add sal_min ----------------------------
## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "sal_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y))), "year", "soft_bottom_clr", "rock_clr", "sal_min_1m")

# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 1,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.2_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.2_nonspatial_cv_R.rds")

# model.non.spatial_cv results 
model.non.spatial_cv$evaluation

# 4  Testing pseudo.r.squared  0.5195000              0.2379573000  0.3495000  0.6482500  0.4865000 0.080 0.1958650045  0.2080000  0.6880000
# 5  Testing        r.squared  0.2730000              0.2631615000  0.1260000  0.4222500  0.2685000 0.074 0.1805522085  0.0430000  0.4730000
# 6  Testing             rmse 13.2845000              0.7153545000 12.9630000 14.2157500 13.5686667 0.346 0.8484389587 12.7220000 14.7390000


#-------------------------- T.0.3 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Replace sal by sal_min ----------------------------
## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add pred1.5$soft_bottom_clr and pred1.5$rock_clr to tot
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr", "sal_min_1m")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness, x,y, sal_mean_1m))), "year", "soft_bottom_clr", "rock_clr", "sal_min_1m")

# Extract predictors
predictors <- tot %>%
  dplyr::select(all_of(col_pred))
sort(names(predictors))

# Extract response variables
ind <- div  %>% dplyr::select("R")



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 1,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.3_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.3_nonspatial_cv_R.rds")

# model.non.spatial_cv results 
#model.non.spatial_cv <- readRDS("./output/models/spatialrf/T.0.1_nonspatial_cv_R.rds")
model.non.spatial_cv$evaluation



## Variable importance ----
spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)
ggsave("./figures/models/spatialrf/T.0.3_non-spatial_variable_importance_R.png")


#-------------------------- T.0.4 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add temp_min ----------------------------
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



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 1,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.4_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.4_nonspatial_cv_R.rds")

# model.non.spatial_cv results 
model.non.spatial_cv <- readRDS("./output/models/spatialrf/T.0.4_nonspatial_cv_R.rds")
model.non.spatial_cv$evaluation

# $aggregated
# model           metric     median median_absolute_deviation         q1         q3       mean    se           sd        min        max
# 1     Full pseudo.r.squared         NA                        NA         NA         NA  0.9791859    NA           NA         NA         NA
# 2     Full        r.squared         NA                        NA         NA         NA  0.9588051    NA           NA         NA         NA
# 3     Full             rmse         NA                        NA         NA         NA  4.6410000    NA           NA         NA         NA
# 4  Testing pseudo.r.squared  0.5565000              0.1371405000  0.3475000  0.6455000  0.4958333 0.073 0.1783383488  0.2600000  0.6490000
# 5  Testing        r.squared  0.3160000              0.1556730000  0.1260000  0.4167500  0.2723333 0.068 0.1659983936  0.0680000  0.4210000
# 6  Testing             rmse 13.1615000              0.9592422000 12.9320000 14.6855000 13.6280000 0.523 1.2816487818 12.1400000 15.2560000
# 7 Training pseudo.r.squared  0.9789903              0.0004464250  0.9786433  0.9792045  0.9789211 0.000 0.0003674503  0.9784341  0.9793085
# 8 Training        r.squared  0.9584219              0.0008742262  0.9577428  0.9588414  0.9582867 0.000 0.0007193702  0.9573332  0.9590452
# 9 Training             rmse  4.6894500              0.0418834500  4.6676500  4.7486750  4.7079667 0.021 0.0523742939  4.6604000  4.7797000
## Variable importance ----
spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)
ggsave("./figures/models/spatialrf/T.0.4_non-spatial_variable_importance_R.png")


#-------------------------- T.0.4.1 :T.0.4 + CV à c(40, 30, 20, 10 km)  --------
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



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## CV evaluation -----
model.non.spatial <- readRDS("./output/models/spatialrf/T.0.4_nonspatial_R.rds")

model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 10000, # X km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

model.non.spatial_cv$evaluation
# 50km mean r.squared on testing set : 0.4958333
# 40km mean pseudo.r.squared on testing set : 0.4891111
# 30km mean pseudo.r.squared on testing set : 0.5212857 
# 20km mean pseudo.r.squared on testing set : 0.5015833
# 10km mean pseudo.r.squared on testing set : 0.5232414


#-------------------------- T.0.5 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 1)---------------------------------------------
#-------------------------- Add temp_min + surface only ----------------------------
## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add predictors
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



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## Prep data for spatial-RF ------------------------------

tot_sf <- st_as_sf(tot, coords = c("x", "y"), crs = 4326)  # adapt names if needed
tot_m  <- st_transform(tot_sf, 3857)                       # metric CRS (meters)

xy_m <- st_coordinates(tot_m) |> as.data.frame()
colnames(xy_m) <- c("x", "y")

nrow(tot)    # 637
nrow(xy_m)   # must be 637
head(xy_m)

dm  <- st_distance(tot_m)             # units [m]
distance_matrix <- as.matrix(drop_units(dm))  # plain numeric

dim(distance_matrix)   # 637 637
range(distance_matrix)

distance_thresholds <- c(0, 1000, 2000, 5000, 10000, 50000, 100000)


## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 1,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.5_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.5_nonspatial_cv_R.rds")

# model.non.spatial_cv results 
model.non.spatial_cv$evaluation

# $aggregated
# model           metric     median median_absolute_deviation         q1         q3       mean    se           sd        min        max
# 1     Full pseudo.r.squared         NA                        NA         NA         NA  0.9820110    NA           NA         NA         NA
# 2     Full        r.squared         NA                        NA         NA         NA  0.9643456    NA           NA         NA         NA
# 3     Full             rmse         NA                        NA         NA         NA  4.5390000    NA           NA         NA         NA
# 4  Testing pseudo.r.squared  0.5830000                         0  0.5830000  0.5830000  0.5726471 0.024 0.0969200322  0.3150000  0.7040000
# 5  Testing        r.squared  0.3400000                         0  0.3400000  0.3400000  0.3368235 0.024 0.0986085920  0.0990000  0.4960000
# 6  Testing             rmse 10.7990000                         0 10.7990000 11.3600000 11.0340588 0.103 0.4247928128 10.7990000 12.3000000
# 7 Training pseudo.r.squared  0.9812016                         0  0.9812016  0.9812016  0.9812733 0.000 0.0003585525  0.9808850  0.9821851
# 8 Training        r.squared  0.9627565                         0  0.9627565  0.9627565  0.9628974 0.000 0.0007039173  0.9621354  0.9646875
# 9 Training             rmse  4.6864000                         0  4.6864000  4.6864000  4.6894941 0.010 0.0416623552  4.6039000  4.7844000




#-------------------------- T.0.6 : non spatial model with tuned parm (mtry = 2, ntree =500, nodesize = 5)---------------------------------------------
#-------------------------- surface only + x,y ----------------------------
## Prep data for model ----

pred1.5 <- st_read("./data/processed_data/predictors/predictors_sel_v1.5_month.gpkg")
sort(names(pred1.5))

# Add predictors
tot <- tot %>%
  left_join(pred1.5 %>% st_drop_geometry() %>% dplyr::select(c("replicates", "soft_bottom_clr", "rock_clr")), by = "replicates")

sort(names(tot))

tot <- tot %>%
  mutate(year = year(date))

col_pred <- c(colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, grouped_main_habitat, northness, eastness))), "year", "soft_bottom_clr", "rock_clr")


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



# --- Define resp and predictors ---
response_name <- names(ind)
covariates_names <- names(predictors)




## Prep data for spatial-RF ------------------------------

tot_sf <- st_as_sf(tot, coords = c("x", "y"), crs = 4326)  # adapt names if needed
tot_m  <- st_transform(tot_sf, 3857)                       # metric CRS (meters)

xy_m <- st_coordinates(tot_m) |> as.data.frame()
colnames(xy_m) <- c("x", "y")

nrow(tot)    # 637
nrow(xy_m)   # must be 637
head(xy_m)

dm  <- st_distance(tot_m)             # units [m]
distance_matrix <- as.matrix(drop_units(dm))  # plain numeric

dim(distance_matrix)   # 637 637
range(distance_matrix)

distance_thresholds <- c(0, 1000, 2000, 5000, 10000, 50000, 100000)


## Non spatial rf -----------
model.non.spatial <- spatialRF::rf(
  data = tot,
  dependent.variable.name = response_name,
  predictor.variable.names = covariates_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = xy_m, #not needed by rf, but other functions read it from the model
  seed = 1,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = 2,
    ntree = 500,
    min.node.size = 5,
    max.depth = NULL
  )
)

saveRDS(model.non.spatial, "./output/models/spatialrf/T.0.6_nonspatial_R.rds")


## CV evaluation -----
model.non.spatial_cv <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy_m,                  #data coordinates
  #repetitions = 30,         # number of spatial folds
  training.fraction = 0.95, # training data fraction on each fold
  metrics = c("r.squared", "pseudo.r.squared", "rmse"),
  distance.step = 50000, # 50 km distance for spatial folds
  seed = 1,
  verbose = FALSE
)

saveRDS(model.non.spatial_cv, "./output/models/spatialrf/T.0.6_nonspatial_cv_R.rds")

# model.non.spatial_cv results 
model.non.spatial_cv$evaluation

# $aggregated
# model           metric     median median_absolute_deviation         q1         q3       mean    se           sd        min        max
# 1     Full pseudo.r.squared         NA                        NA         NA         NA  0.9692174    NA           NA         NA         NA
# 2     Full        r.squared         NA                        NA         NA         NA  0.9393824    NA           NA         NA         NA
# 3     Full             rmse         NA                        NA         NA         NA  5.6648000    NA           NA         NA         NA
# 4  Testing pseudo.r.squared  0.6460000                         0  0.6460000  0.6460000  0.6291176 0.024 0.1000255482  0.3300000  0.7410000
# 5  Testing        r.squared  0.4180000                         0  0.4180000  0.4180000  0.4055882 0.026 0.1071191269  0.1090000  0.5480000
# 6  Testing             rmse 10.6520000                         0 10.6520000 10.8280000 10.8357059 0.087 0.3585219248 10.6520000 11.6340000
# 7 Training pseudo.r.squared  0.9695230                         0  0.9695230  0.9695230  0.9694996 0.000 0.0005057549  0.9683549  0.9703716
# 8 Training        r.squared  0.9399748                         0  0.9399748  0.9399748  0.9399296 0.000 0.0009804275  0.9377113  0.9416211
# 9 Training             rmse  5.7380000                         0  5.7380000  5.7380000  5.7403118 0.007 0.0302062055  5.6559000  5.7989000



