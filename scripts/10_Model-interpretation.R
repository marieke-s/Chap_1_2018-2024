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

# model T.0.0 ----
models_list <- readRDS("./output/models/T0.0.rds")



#------------- Prep ----------
df <- pred %>%
  st_drop_geometry() %>%
  left_join(div, by = "replicates")

target_col <- colnames(div %>% dplyr::select(-replicates))
feature_cols <- colnames(pred %>% st_drop_geometry() %>% dplyr::select(-c(replicates, x, y, grouped_main_habitat)))

# Separate numeric and categorical predictors
numeric_features <- df %>% dplyr::select(all_of(feature_cols)) %>% dplyr::select(where(is.numeric)) %>% colnames()
factor_features  <-  df %>% dplyr::select(all_of(feature_cols)) %>% dplyr::select(where(is.character)) %>% colnames() # none here 


# Convert character to factor if need 
if (length(factor_features) > 0){
  df[factor_features] <- lapply(df[factor_features], factor)
  message("Converted character columns to factors.")
} else {
  message("No character columns to convert.")
}


#------------- Response var distribution ----------
# Histogram for target (response variables)

for (t in target_col){
  p <- ggplot(df, aes(x = .data[[t]])) +
    geom_histogram(bins = 50) +
    labs(title = paste("Distribution of", t))
  
  print(p)
}


#------------- Model-independent associations ----------
### Spearman correlations (numeric predictors vs numeric target) ----

target_col <- "R"

spearman_results <- purrr::map_dfr(numeric_features, function(v) {
  x <- df[[v]]
  y <- df[[target_col]]
  
  if (all(is.na(x))) return(NULL)
  
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  
  tibble(
    feature   = v,
    method    = "Spearman",
    estimate  = unname(ct$estimate),
    p_value   = ct$p.value
  )
})


## Plot Spearman correlations ----
spearman_plot <- spearman_results %>%
  ggplot(aes(x = reorder(feature, estimate), y = estimate)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = paste("Spearman Correlations with", target_col),
       x = "Features",
       y = "Spearman Correlation Coefficient") +
  theme_minimal()


## Kruskal-Wallis for categorical predictors vs numeric target ----

kw_results <- purrr::map_dfr(factor_features, function(v) {
  x <- df[[v]]
  y <- df[[target_col]]
  # Need at least 2 levels
  if (nlevels(x) < 2) return(NULL)
  kw <- kruskal.test(y ~ x)
  tibble(
    feature   = v,
    method    = "Kruskal-Wallis",
    estimate  = NA_real_,             # H statistic is kw$statistic
    statistic = unname(kw$statistic),
    p_value   = kw$p.value
  )
})

stat_results <- bind_rows(spearman_results, kw_results)

#------------- XGBoost built-in feature importance (model-based, biased) -------
## For one fold ----
xgb_model <- models_list$R$XGB$bloo50$Fold_1$model
class(xgb_model)

xgb_importance <- xgb.importance(
  model = xgb_model,
  feature_names = feature_cols
)

# Plot
xgb.plot.importance(xgb_importance, top_n = 20)



## Averaging across folds ----
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


# Make the plot -----
ggplot(mean_importance, aes(x = reorder(Feature, mean_gain), y = mean_gain)) +
  geom_col() +
  geom_text(
    aes(label = round(mean_gain, 3)),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  theme_minimal() +
  ylim(0, max(mean_importance$mean_gain) * 1.15)  # Add space for text


# Save the plot 
ggsave(
  filename = "./figures/models/T0.0/Var_imp_T0.0_Crypto.png",
  width    = 8,
  height   = 6
)














# Extract importance from each fold's model
importance_list <- lapply(folds, function(fold) {
  xgb.importance(model = fold$model)  # <-- change 'model' to actual element name
})

# Combine all importances into one table
all_importance <- bind_rows(importance_list, .id = "Fold")

# Average importance per variable
mean_importance <- all_importance |>
  group_by(Feature) |>
  summarise(mean_gain = mean(Gain, na.rm = TRUE), .groups = "drop",
            mean_cover = mean(Cover, na.rm = TRUE), .groups = "drop", 
            mean_frequency = mean(Frequency, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_gain))



# Definition of the metrics:
# Gain : contribution of each feature to the model. For boosted tree model, each gain of each feature of each tree is taken into account, then average per feature to give a vision of the entire model. Highest percentage means important feature to predict the label used for the training (only available for tree models);

# Cover : metric of the number of observation related to this feature (only available for tree models)

# Weight : percentage representing the relative number of times a feature have been taken into trees.

#------------- SHAP values (model-based, inherits XGBoost bias) ----
## For one fold ----
# Build the design matrix X from feature columns
X_train <- as.matrix(df[, feature_cols])

# Use training matrix (or full X) for SHAP
shap_vals <- SHAPforxgboost::shap.values(
  xgb_model = models_list$R$XGB$bloo50$Fold_1$model,
  X_train   = X_train
)

# Mean |SHAP| importance per feature
shap_importance <- data.frame(
  feature        = names(shap_vals$mean_shap_score),
  mean_abs_SHAP  = shap_vals$mean_shap_score,
  row.names      = NULL
)

shap_importance[order(-shap_importance$mean_abs_SHAP), ]

shap_long <- SHAPforxgboost::shap.prep(
  xgb_model = xgb_model,
  X_train   = X_train
)

plot <- SHAPforxgboost::shap.plot.summary(shap_long)

# Save the plot
ggsave(
  filename = "./figures/models/T0.0/SHAP_T.0.0_$R$XGB$bloo50$Fold_1$.png",
  plot     = plot,
  width    = 8,
  height   = 6
)



# SHAP dependence plots ----
library(SHAPforxgboost)
library(gridExtra)

## 1. Get the feature names in the order of mean |SHAP| importance
features <- names(shap_vals$mean_shap_score)  

## 2. Split feature list into chunks of 4
feature_chunks <- split(features, ceiling(seq_along(features) / 4))

## 3. Loop over chunks, make 1 plot per feature, arrange in 2x2 panels, and save
for (i in seq_along(feature_chunks)) {
  
  this_chunk <- feature_chunks[[i]]
  
  # one dependence plot per feature
  fig_list <- lapply(
    this_chunk,
    function(f) {
      SHAPforxgboost::shap.plot.dependence(
        data_long = shap_long,
        x         = f,
        dilute    = 5
      )
    }
  )
  
  # open a graphics device to save the panel
  png(
    filename = sprintf("./fugures/models/T0.0/shap_dependence_panel_%02d.png", i),
    width    = 2000,
    height   = 2000,
    res      = 300
  )
  
  # arrange up to 4 plots (2 x 2)
  gridExtra::grid.arrange(grobs = fig_list, ncol = 2)
  
  dev.off()
}


























##  Build one dependence plot per feature
#    We rely on shap.plot.dependence returning a ggplot object
plot_list <- lapply(feature_cols, function(f) {
  p <- SHAPforxgboost::shap.plot.dependence(
    data_long = shap_long,
    x         = f
    # you can also add: color_feature = f or another feature if you like
  )
  # Add a clearer title
  p + ggtitle(paste("SHAP dependence:", f))
})

names(plot_list) <- features_cols

## 3. Arrange plots in panels of 4 (2x2) and save to a multi-page PDF
n_per_page <- 4
n_pages <- ceiling(length(plot_list) / n_per_page)

pdf("shap_dependence_panels.pdf", width = 10, height = 8)  # adjust size as you wish

for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * n_per_page + 1
  idx_end   <- min(i * n_per_page, length(plot_list))
  plots_this_page <- plot_list[idx_start:idx_end]
  
  # Arrange in 2 x 2 grid; if fewer than 4 plots on last page, it still works
  gridExtra::grid.arrange(grobs = plots_this_page, ncol = 2, nrow = 2)
}

dev.off()


# features vs R + correlation -----

library(ggplot2)
library(gridExtra)

# ---------------------------------------------------------------
# 1. Define target and predictor features
# ---------------------------------------------------------------

target_col <- "R"

# numeric predictors only (you can change if needed)
numeric_features <- setdiff(
  names(df)[sapply(df, is.numeric)],
  target_col
)

# Helper function: compute correlations + label
cor_label <- function(x, y) {
  # Pearson
  p_cor <- suppressWarnings(cor.test(x, y, method = "pearson"))
  pearson_r  <- round(p_cor$estimate, 3)
  pearson_p  <- signif(p_cor$p.value, 3)
  
  # Spearman
  s_cor <- suppressWarnings(cor.test(x, y, method = "spearman"))
  spearman_r <- round(s_cor$estimate, 3)
  spearman_p <- signif(s_cor$p.value, 3)
  
  label <- paste0(
    "Pearson r = ", pearson_r, " (p=", pearson_p, ")\n",
    "Spearman ρ = ", spearman_r, " (p=", spearman_p, ")"
  )
  
  return(label)
}

# ---------------------------------------------------------------
# 2. Create one plot per feature: feature vs R + correlation text
# ---------------------------------------------------------------

make_feature_plot <- function(feature) {
  
  x <- df[[feature]]          
  y <- df[[target_col]]
  
  label_text <- cor_label(x, y)
  
  ggplot(df, aes_string(x = feature, y = target_col)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed") +
    labs(
      title    = paste("Feature:", feature),
      subtitle = label_text,
      x        = feature,
      y        = target_col
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 8),       # smaller title
      plot.subtitle = element_text(size = 6)      # optional: smaller subtitle
    )
}


# ---------------------------------------------------------------
# 3. Split features into chunks of 4 for panel plots
# ---------------------------------------------------------------

feature_chunks <- split(numeric_features, ceiling(seq_along(numeric_features) / 4))

# ---------------------------------------------------------------
# 4. Loop over chunks, make panels, save PNGs
# ---------------------------------------------------------------

for (i in seq_along(feature_chunks)) {
  
  this_chunk <- feature_chunks[[i]]
  
  fig_list <- lapply(this_chunk, make_feature_plot)
  
  png(
    filename = sprintf("./figures/models/T0.0/R-features_correlation_%02d.png", i),
    width    = 2000,
    height   = 2000,
    res      = 300
  )
  
  gridExtra::grid.arrange(grobs = fig_list, ncol = 2)
  
  dev.off()
}

