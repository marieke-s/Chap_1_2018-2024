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
library(grf)




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




############################################################
## Script 1: Run causal forest and export results
## File: run_causal_forest.R
############################################################

## Packages
library(grf)
library(sf)
library(dplyr)

##----------------------------------------------------------
## 0. User inputs: where to save results
##----------------------------------------------------------
results_dir <- "./output/models/causal_forest/gravity/"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

##----------------------------------------------------------
## 1. Data: pred (sf) and div$R must already be in environment
##----------------------------------------------------------
# Check objects
if (!exists("pred")) stop("Object 'pred' not found.")
if (!exists("div"))  stop("Object 'div' not found.")
if (!("R" %in% names(div))) stop("div$R not found.")

Y <- div$R
if (!is.numeric(Y)) stop("div$R must be numeric.")
if (nrow(pred) != length(Y)) stop("nrow(pred) and length(div$R) differ.")

##----------------------------------------------------------
## 2. Build covariate matrix X
##----------------------------------------------------------
pred_nongeom <- sf::st_drop_geometry(pred)

X <- pred_nongeom %>%
  dplyr::select(
    c('cop_chl_month_mean_log', 'shore_dist_m_weight_log', 'canyon_dist_m_weight_log', 'port_dist_m_weight', 'bathy_mean')
  )

names(X)

X_mat <- as.matrix(X)

##----------------------------------------------------------
## 3. Define treatment variable W  (EDIT THIS PART)
##----------------------------------------------------------
## You must define your treatment variable here.
## Example (replace with your real definition):

# Example placeholder: using a dummy binary treatment from some column
# W <- as.numeric(pred_nongeom$some_binary_column)
names(pred_nongeom)
W <- pred_nongeom$gravity_mean

if (!exists("W")) {
  stop("You must define W (numeric treatment variable) before running this script.")
}

if (length(W) != nrow(pred)) stop("W must have same length as rows in pred.")
if (!is.numeric(W)) stop("W must be numeric (binary or continuous).")

##----------------------------------------------------------
## 4. Pre-fit Y.hat and W.hat
##----------------------------------------------------------
set.seed(123)

cat("Fitting regression forest for Y.hat...\n")
forest_Y <- regression_forest(X_mat, Y, tune.parameters = "all", num.trees = 200)
Y_hat <- predict(forest_Y)$predictions

cat("Fitting regression forest for W.hat...\n")
forest_W <- regression_forest(X_mat, W, tune.parameters = "all", num.trees = 200)
W_hat <- predict(forest_W)$predictions

##----------------------------------------------------------
## 5. Fit causal forest
##----------------------------------------------------------
cat("Fitting causal forest...\n")
cf <- causal_forest(
  X            = X_mat,
  Y            = Y,
  W            = W,
  Y.hat        = Y_hat,
  W.hat        = W_hat,
  num.trees    = 2000,
  tune.parameters = "all"
)






#------------------------- With train / test ---------------------

##----------------------------------------------------------
## 1. Data: pred (sf) and div$R must already be in environment
##----------------------------------------------------------
# Check objects
if (!exists("pred")) stop("Object 'pred' not found.")
if (!exists("div"))  stop("Object 'div' not found.")
if (!("R" %in% names(div))) stop("div$R not found.")

Y <- div$R
if (!is.numeric(Y)) stop("div$R must be numeric.")
if (nrow(pred) != length(Y)) stop("nrow(pred) and length(div$R) differ.")

##----------------------------------------------------------
## 2. Build covariate matrix X
##----------------------------------------------------------
pred_nongeom <- sf::st_drop_geometry(pred)

X <- pred_nongeom %>%
  dplyr::select(
    c(
      "cop_chl_month_mean_log",
      "shore_dist_m_weight_log",
      "canyon_dist_m_weight_log",
      "port_dist_m_weight",
      "bathy_mean"
    )
  )

names(X)

X_mat <- as.matrix(X)

##----------------------------------------------------------
## 3. Define treatment variable W (continuous gravity_mean)
##----------------------------------------------------------
names(pred_nongeom)
W <- pred_nongeom$gravity_mean

if (!exists("W")) {
  stop("You must define W (numeric treatment variable) before running this script.")
}

if (length(W) != nrow(pred)) stop("W must have same length as rows in pred.")
if (!is.numeric(W)) stop("W must be numeric (binary or continuous).")


## 4. Train–test split ---------------------

set.seed(0)  # for reproducibility

n <- nrow(X_mat)
test_frac <- 0.2
test_size <- floor(test_frac * n)

test_idx  <- sample(1:n, size = test_size)
train_idx <- setdiff(1:n, test_idx)

X_train <- X_mat[train_idx, , drop = FALSE]
X_test  <- X_mat[test_idx,  , drop = FALSE]

Y_train <- Y[train_idx]
Y_test  <- Y[test_idx]

W_train <- W[train_idx]
W_test  <- W[test_idx]  # not strictly needed, but kept for completeness

##----------------------------------------------------------
## 5. Nuisance models for outcome and treatment (on TRAIN only)
##----------------------------------------------------------
set.seed(123)

cat("Fitting regression forest for Y.hat on TRAIN data...\n")
forest_Y <- regression_forest(
  X_train,
  Y_train,
  tune.parameters = "all",
  num.trees = 200
)
Y_hat_train <- predict(forest_Y)$predictions

cat("Fitting regression forest for W.hat on TRAIN data...\n")
forest_W <- regression_forest(
  X_train,
  W_train,
  tune.parameters = "all",
  num.trees = 200
)
W_hat_train <- predict(forest_W)$predictions

##----------------------------------------------------------
## 6. Fit causal forest on TRAIN data
##----------------------------------------------------------
cat("Fitting causal forest on TRAIN data...\n")
cf <- causal_forest(
  X            = X_train,
  Y            = Y_train,
  W            = W_train,
  Y.hat        = Y_hat_train,
  W.hat        = W_hat_train,
  num.trees    = 2000,
  tune.parameters = "all"
)

print(cf)

## Optional: ATE on train set
ate_train <- average_treatment_effect(cf)
cat("\nATE on training data:\n")
print(ate_train)

##----------------------------------------------------------
## 7. Predict CATEs (treatment effects) on TEST data
##----------------------------------------------------------
# This is analogous to: t_effect_test = cf.effect(X_test) in Python
t_effect_test <- predict(cf, X_test)$predictions  # vector of CATEs on test set

cat("\nFirst 10 estimated treatment effects on TEST data:\n")
print(head(t_effect_test, 10))

## Optionally, inspect distribution on test set
hist(
  t_effect_test,
  breaks = 30,
  main = "Estimated treatment effects on TEST set",
  xlab = "tau(x) on test"
)

## Optionally, compare with Y_test (purely exploratory, not causal)
cat("\nSummary of Y_test:\n")
print(summary(Y_test))



##----------------------------------------------------------
## 6. Extract ATE, tau_hat, variable importance
##----------------------------------------------------------
ate <- average_treatment_effect(cf)   # named vector with 'estimate' and 'std.err'
tau_hat <- as.numeric(predict(cf)$predictions)

vi <- variable_importance(cf)
vi_table <- data.frame(
  variable   = colnames(X_mat),
  importance = vi
) %>%
  arrange(desc(importance))

## Also some summary stats useful for interpretation
summary_stats <- list(
  Y_mean = mean(Y),
  Y_sd   = sd(Y),
  tau_mean = mean(tau_hat),
  tau_sd   = sd(tau_hat),
  n = length(Y)
)

##----------------------------------------------------------
## 7. Export results as RDS/CSV
##----------------------------------------------------------
saveRDS(cf,       file.path(results_dir, "causal_forest_model.rds"))
saveRDS(ate,      file.path(results_dir, "causal_forest_ate.rds"))
saveRDS(tau_hat,  file.path(results_dir, "causal_forest_tau_hat.rds"))
saveRDS(vi_table, file.path(results_dir, "causal_forest_vi_table.rds"))
saveRDS(summary_stats, file.path(results_dir, "causal_forest_summary_stats.rds"))

# Optional: export CSVs
write.csv(vi_table,
          file = file.path(results_dir, "causal_forest_variable_importance.csv"),
          row.names = FALSE)

pred_with_tau <- pred
pred_with_tau$tau_hat <- tau_hat

# Optional: export as GeoPackage (or shapefile)
sf::st_write(pred_with_tau,
             dsn   = file.path(results_dir, "causal_forest_tau_hat.gpkg"),
             layer = "tau_hat",
             delete_layer = TRUE)

cat("\nAll causal forest results saved in directory:", results_dir, "\n")
############################################################
## End of run_causal_forest.R
############################################################

############################################################
## Script 2: Interpret causal forest results
## File: interpret_causal_forest.R
############################################################

results_dir <- "./output/models/causal_forest/bathy_mean/"

##----------------------------------------------------------
## 1. Load saved objects
##----------------------------------------------------------
ate_path      <- file.path(results_dir, "causal_forest_ate.rds")
tau_path      <- file.path(results_dir, "causal_forest_tau_hat.rds")
vi_path       <- file.path(results_dir, "causal_forest_vi_table.rds")
summary_path  <- file.path(results_dir, "causal_forest_summary_stats.rds")

if (!file.exists(ate_path))     stop("ATE file not found: ", ate_path)
if (!file.exists(tau_path))     stop("tau_hat file not found: ", tau_path)
if (!file.exists(vi_path))      stop("vi_table file not found: ", vi_path)
if (!file.exists(summary_path)) stop("summary_stats file not found: ", summary_path)

ate          <- readRDS(ate_path)          # named vector: estimate, std.err
tau_hat      <- readRDS(tau_path)          # numeric vector
vi_table     <- readRDS(vi_path)           # data.frame
summary_stats <- readRDS(summary_path)     # list

##----------------------------------------------------------
## 2. Compute basic quantities
##----------------------------------------------------------
est <- as.numeric(ate["estimate"])
se  <- as.numeric(ate["std.err"])

ci_low  <- est - 1.96 * se
ci_high <- est + 1.96 * se

tau_sd  <- summary_stats$tau_sd
Y_sd    <- summary_stats$Y_sd
n       <- summary_stats$n

##----------------------------------------------------------
## 3. Build interpretation text with if/else logic
##----------------------------------------------------------
interpret_lines <- c()

interpret_lines <- c(
  interpret_lines,
  "Causal forest interpretation",
  "============================",
  sprintf("Number of observations: %d", n),
  "",
  "Average Treatment Effect (ATE)",
  "------------------------------",
  sprintf("ATE estimate: %.3f", est),
  sprintf("Standard error: %.3f", se),
  sprintf("95%% CI: [%.3f, %.3f]", ci_low, ci_high),
  ""
)

## Interpretation of statistical significance
if (ci_low > 0) {
  interpret_lines <- c(
    interpret_lines,
    "Interpretation: The estimated treatment effect is **positive** and",
    "statistically significant at the 95% confidence level (CI entirely > 0).",
    ""
  )
} else if (ci_high < 0) {
  interpret_lines <- c(
    interpret_lines,
    "Interpretation: The estimated treatment effect is **negative** and",
    "statistically significant at the 95% confidence level (CI entirely < 0).",
    ""
  )
} else {
  interpret_lines <- c(
    interpret_lines,
    "Interpretation: The 95% confidence interval includes 0.",
    "There is **no strong evidence** of a non-zero average treatment effect;",
    "the effect is statistically uncertain.",
    ""
  )
}

## Magnitude relative to outcome scale (if Y_sd available)
if (!is.null(Y_sd) && !is.na(Y_sd) && Y_sd > 0) {
  rel_est <- est / Y_sd
  interpret_lines <- c(
    interpret_lines,
    "Magnitude relative to outcome variability",
    "-----------------------------------------",
    sprintf("Standard deviation of outcome (Y): %.3f", Y_sd),
    sprintf("ATE / sd(Y): %.3f", rel_est),
    ""
  )
  
  if (abs(rel_est) < 0.2) {
    interpret_lines <- c(
      interpret_lines,
      "The ATE is small relative to the natural variability in Y (|ATE| < 0.2 * sd(Y)).",
      "Even if statistically significant, the practical effect size appears small.",
      ""
    )
  } else if (abs(rel_est) < 0.5) {
    interpret_lines <- c(
      interpret_lines,
      "The ATE is of moderate magnitude relative to sd(Y).",
      ""
    )
  } else {
    interpret_lines <- c(
      interpret_lines,
      "The ATE is large relative to the natural variability in Y (|ATE| > 0.5 * sd(Y)).",
      ""
    )
  }
}

##----------------------------------------------------------
## 4. Heterogeneity in treatment effects
##----------------------------------------------------------
interpret_lines <- c(
  interpret_lines,
  "Heterogeneity in treatment effects",
  "----------------------------------",
  sprintf("Mean of tau_hat (CATEs): %.3f", summary_stats$tau_mean),
  sprintf("SD of tau_hat: %.3f", tau_sd),
  ""
)

# Heuristic: if sd(tau_hat) is big compared to |ATE|, suggest heterogeneity
if (!is.na(tau_sd) && tau_sd > 0) {
  if (abs(est) < 1e-6) {
    # ATE near zero, judge heterogeneity by raw sd
    if (tau_sd > 0.1 * Y_sd) {
      interpret_lines <- c(
        interpret_lines,
        "The average effect is close to zero, but the standard deviation of the",
        "individual treatment effects (tau_hat) is relatively large. This suggests",
        "potential heterogeneity even though the mean effect is weak.",
        ""
      )
    } else {
      interpret_lines <- c(
        interpret_lines,
        "Both the average effect and the variability in tau_hat are small, suggesting",
        "limited evidence of strong treatment-effect heterogeneity.",
        ""
      )
    }
    
  } else {
    ratio_sd_ate <- tau_sd / abs(est)
    interpret_lines <- c(
      interpret_lines,
      sprintf("SD(tau_hat) / |ATE|: %.3f", ratio_sd_ate),
      ""
    )
    
    if (ratio_sd_ate < 0.5) {
      interpret_lines <- c(
        interpret_lines,
        "The variability of tau_hat is small relative to the mean effect.",
        "There is limited evidence for strong heterogeneity in treatment effects.",
        ""
      )
    } else if (ratio_sd_ate < 1.5) {
      interpret_lines <- c(
        interpret_lines,
        "The variability of tau_hat is comparable to the mean effect.",
        "There is some evidence of heterogeneity in treatment effects.",
        ""
      )
    } else {
      interpret_lines <- c(
        interpret_lines,
        "The variability of tau_hat is much larger than the mean effect.",
        "This indicates substantial heterogeneity: the impact of treatment",
        "may differ a lot between locations / covariate profiles.",
        ""
      )
    }
  }
}

##----------------------------------------------------------
## 5. Variable importance interpretation
##----------------------------------------------------------
interpret_lines <- c(
  interpret_lines,
  "Variables driving heterogeneity (variable importance)",
  "----------------------------------------------------",
  ""
)

# Sort and take top 5 variables
vi_table <- vi_table[order(vi_table$importance, decreasing = TRUE), ]
top_k <- min(5, nrow(vi_table))

interpret_lines <- c(
  interpret_lines,
  sprintf("Top %d variables by importance (for treatment-effect heterogeneity):", top_k),
  paste(
    apply(vi_table[1:top_k, ], 1, function(row) {
      sprintf("  - %s (importance = %.3f)", row[["variable"]], as.numeric(row[["importance"]]))
    }),
    collapse = "\n"
  ),
  "",
  "Note: These variables are important for explaining **where** the treatment effect varies,",
  "not necessarily for predicting Y itself.",
  ""
)

##----------------------------------------------------------
## 6. Print and save report
##----------------------------------------------------------
cat(paste(interpret_lines, collapse = "\n"), "\n")

report_path <- file.path(results_dir, "causal_forest_interpretation.txt")
writeLines(interpret_lines, con = report_path)

cat("\nInterpretation written to:\n", report_path, "\n")
############################################################
## End of interpret_causal_forest.R
############################################################


