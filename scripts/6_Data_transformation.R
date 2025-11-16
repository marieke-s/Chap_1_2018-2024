#------------- Description ----
# The aim of this script is to explore all data (Mtdt, species, predictos, div_indices both raw and transformed)





#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(compositions)
library(dplyr)
library(exactextractr)
library(sf)
library(raster)
library(ncdf4)
library(patchwork)
library(lubridate)
library(terra)
library(stringr)
library(pMEM)
library(ggplot2)




#------------- Load functions ------------------------------------------------
source("./utils/Fct_Data-Prep.R")

#------------- Load data ----
# Mtdt_7
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")

# predictors_raw_v2.0
pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v2.1.gpkg")

# div_indices_v1.0
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

























#----------------------------------------- Predictors transformation ------------------
#--- Remove 0 variance predictors ------------------
nzv <- caret::nearZeroVar(pred_raw %>% st_drop_geometry(), saveMetrics = TRUE)

# zeroVar cols
nzv_cols <- rownames(nzv[nzv$zeroVar == TRUE, ]) # "wind_min_7j" "wind_min_1m" "ws_min_1y"  

# Remove zero variance predictors
pred_tr <- pred_raw %>%
  dplyr::select(-all_of(nzv_cols)) # -3 predictors

rm(nzv, nzv_cols)

#--- Replace negative distance values by 0 ------------------
# Distance between seabed and depth sampling  contains negative values due to the difference in bathymetry sources used to compute this variable : depth seafloor is retrieved from SHOM MNT extraction at the buffer level, depth sampling comes from observation of scientist on field with varying measurement instruments (eg. diving coputer). 

# Check negative dist_seabed_depthsampling
pred_tr %>% st_drop_geometry() %>%
  left_join(mtdt %>% st_drop_geometry() %>% dplyr::select(replicates, method), by = "replicates") %>%
  filter(dist_seabed_depthsampling < 0) %>%
  pull(method) %>%
  unique() # dive_transect, dive_motionless and motionless_descended_from_surface --> In all these methods the sampling can occur very close from seafloor --> negative values can be transformed in 0.

# Set negative dist_seabed_depthsampling to 0
pred_tr <- pred_tr %>%
  mutate(dist_seabed_depthsampling = ifelse(dist_seabed_depthsampling < 0, 0, dist_seabed_depthsampling))

#--- CLT of habitat composition ---------
# For the variables : Surface proportion of each habitat within buffer
# --> we need to apply transformation to avoid compositional data bias (see here for details : https://docs.google.com/document/d/1cN9vJ6I4fHzPXZfXjOm77Klk5hCSFBBkOj6Mhgum_S8/edit?tab=t.0 and https://medium.com/@nextgendatascientist/a-guide-for-data-scientists-log-ratio-transformations-in-machine-learning-a2db44e2a455) 
# We apply Centered log-ratio transformation.
# Aitchison, J. (1986) The Statistical Analysis of Compositional Data, Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). 416p. (https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)

# Habitat columns
hab_cols <- c(
  "habitats_artificiels_mean",
  "matte_morte_p_oceanica_mean",
  "algues_infralittorales_mean",
  "soft_bottom_mean",
  "meadow_mean",
  "rock_mean",
  "coralligenous_mean"
)

# Check data
summary(pred_tr[, hab_cols])

# Format data
X <- as.data.frame(pred_tr[, hab_cols])
X <- data.frame(lapply(X, function(col) as.numeric(as.character(col))))
X <- as.matrix(X) + 1e-6  # pseudocount to avoid zeros

# Centered log-ratio transformation
Z <- compositions::clr(acomp(X))
Z <- as.matrix(Z)  # ensures it's 2D even if only one row

# Set column names
Z <- Z[, colnames(Z) != "geom", drop = FALSE]
colnames(Z) <- sub("_mean$", "_clr", hab_cols)

# remove original mean columns and add CLR-transformed ones
pred_tr <- pred_tr[, setdiff(names(pred_tr), hab_cols)]
pred_tr <- cbind(pred_tr, as.data.frame(Z))


# Plots new hab_cols
for (i in colnames(Z)) {
  hist(Z[, i], main = i, breaks = 100)
}

# Check NA  in all cols 
sapply(pred_tr, function(col) sum(is.na(col))) # no NA --> need to handle habitat NA points (5 of them)

rm(X, Z, hab_cols, i)
dev.off()


#--- Aspect : northness and eastness ----
# Explanation from CHATGPT : https://chatgpt.com/s/t_6914b4ea37a48191807bced590285dc3

pred_tr <- pred_tr %>%
  mutate(northness = cos(aspect_mean), 
         eastness = sin(aspect_mean)) %>%
  dplyr::select(-c(aspect_mean, aspect_min, aspect_max))





#--- Log  ------------------
# Transform variables with outliers using max value >/= to 10*median value rule
# Iterate on all predictors, if max value >/= to 10*median value, apply log(x + 1) transformation, else keep the same
# Create a new data frame to store transformed variables

# Detect numeric columns
pred_num <- pred_tr %>% st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-c("x", "y"))  # Do not log coordinates

pred_num_log <- pred_num  # Make a copy to modify

for (i in colnames(pred_num)) {
  # Compute max and median
  max_val <- max(pred_num[[i]], na.rm = TRUE)
  median_val <- median(pred_num[[i]], na.rm = TRUE)
  
  # Apply log transformation if the condition is met
  if (max_val >= 10 * median_val) {
    print(paste("Transforming:", i))  # Print transformed column name
    new_col_name <- paste0(i, "_log")  # New name with _log suffix
    
    # Replace original column with log-transformed version (renamed)
    pred_num_log <- pred_num_log[, !colnames(pred_num_log) %in% i]  # Drop old column
    pred_num_log[[new_col_name]] <- log(pred_num[[i]] + 1)  # Add new transformed column
  }
}

# Print the nb of transformed columns
print(paste("Number of log-transformed columns:", sum(grepl("_log$", colnames(pred_num_log))))) # 51 
colnames(pred_num_log)

# Replace original numeric columns in pred_tr with transformed ones
pred_tr <- pred_tr %>%
  st_drop_geometry() %>%
  dplyr::select(-colnames(pred_num)) %>%  # Drop the original numeric columns
  dplyr::bind_cols(pred_num_log) %>%       # Add the transformed numeric columns
  dplyr::bind_cols(st_geometry(pred_tr)) %>% # Reattach geometry
  st_as_sf() # Convert back to sf 

pred_tr <- pred_tr %>% dplyr::select(-"...143")
rm(pred_num, pred_num_log, i, max_val, median_val, new_col_name)


#--- Scale -------------------------------
# Explanation
# Why scale and normalize variables?
# First, they can improve the performance and accuracy of many machine learning algorithms, such as linear regression, logistic regression, k-means clustering, and principal component analysis. These algorithms rely on the assumption that the features have similar ranges and distributions, and are sensitive to large differences in scale or magnitude. 
# Second, they can make the data more interpretable and comparable, by eliminating the influence of units or scales. For example, if you want to compare the heights and weights of different people, you need to scale them to a common unit, such as meters and kilograms. 
# Third, they can reduce the risk of numerical errors or overflow, by limiting the range of values and avoiding extremely large or small numbers.

# source : https://www.linkedin.com/advice/3/what-most-effective-techniques-scaling-normalizing-57jlc

# Explo : compare scaling methods ----
# Select numerical predictors 
dt <- pred_tr %>% st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-c("x", "y"))  # Do not scale coordinates

#---- Plot predictors distribution 
plot_distributions(dt, 
                   output_path = "./figures/Predictors/Hist_transfo_scaling/", 
                   version_number = "predictors_raw_v2.1")



# Observations
# Highly skewed data (mostly to the left) with some outliers to the right.

#---- Scale and plot again 
# Standardization
# Explanation : Standardization rescales values to have a mean of 0 and a standard deviation of 1, by subtracting the mean and dividing by the standard deviation. 
# Convert scaled matrix to data frame
dt_standard <- as.data.frame(scale(dt))
plot_distributions(dt_standard, hist_color = "gold",
                   output_path = "./figures/Predictors/Hist_transfo_scaling/", 
                   version_number = "predictors_raw_v2.1_strandard_scaling")


# Min Max scale (values between 0 and 1)
# Explanation : Min-max scaling rescales values to a range between 0 and 1, by subtracting the minimum value and dividing by the range. 
process <- caret::preProcess(dt, method=c("range"))
dt_minmax <- predict(process, dt)

plot_distributions(dt_minmax, hist_color = "pink",
                   output_path = "./figures/Predictors/Hist_transfo_scaling/", 
                   version_number = "predictors_raw_v2.1_minmax_scaling")


# Robust scaling 
# Explanation : Robust scaling rescales values to have a median of 0 and an interquartile range of 1, by subtracting the median and dividing by the interquartile range. The advantage of robust scaling is that it is less sensitive to outliers than min-max scaling or standardization.
robust <- apply(dt, 2, function(x) {
  x <- as.numeric(x)  # ensure numeric
  if (all(is.na(x))) return(rep(NA, length(x)))  # handle all-NA columns
  (x - median(x, na.rm = TRUE)) / IQR(x, na.rm = TRUE)
})
dt_robust <- as.data.frame(robust)
plot_distributions(dt_robust, hist_color = "lightgreen",
                   output_path = "./figures/Predictors/Hist_transfo_scaling/", 
                   version_number = "predictors_raw_v2.1_robust_scaling")


# All
# Define transparent color helper
color_transparent <- function(color, alpha) {
  rgb_col <- col2rgb(color)
  rgb(rgb_col[1]/255, rgb_col[2]/255, rgb_col[3]/255, alpha)
}

# Parameters
colors <- c("lightblue", "gold", "pink", "lightgreen")
alpha <- 0.3
colors <- sapply(colors, color_transparent, alpha = alpha)

# Combine datasets into a list
data_list <- list(dt, dt_standard, dt_minmax, dt_robust)

# Set up parameters
n_predictors <- ncol(dt)
n_per_page <- 8
n_pages <- ceiling(n_predictors / n_per_page)

# Output directory and version number
output_path <- "./figures/Predictors/Hist_transfo_scaling/"
version_number <-  "predictors_raw_v2.1_all_scaling"

# Create output folder if it doesn't exist
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Loop through pages
for (page in 1:n_pages) {
  start <- (page - 1) * n_per_page + 1
  end <- min(page * n_per_page, n_predictors)
  predictors <- colnames(dt)[start:end]
  
  # Define filename for this page
  file_name <- file.path(output_path, paste0("overlay_v", version_number, "_page", page, ".png"))
  
  # Open PNG device (high resolution)
  png(filename = file_name, width = 1600, height = 900)
  par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))  # Adjust margins if needed
  
  # Plot histograms
  for (predictor in predictors) {
    # First dataset
    hist(
      data_list[[1]][[predictor]], 
      main = predictor, 
      xlab = predictor, 
      col = colors[1], 
      border = "white", 
      breaks = 20, 
      freq = FALSE
    )
    
    # Overlay remaining datasets
    for (i in 2:length(data_list)) {
      hist(
        data_list[[i]][[predictor]], 
        col = colors[i], 
        border = NA, 
        breaks = 20, 
        freq = FALSE, 
        add = TRUE
      )
    }
  }
  
  dev.off()
  message("Saved: ", file_name)
}


rm(alpha, colors, end, file_name, i, n_pages, n_per_page, n_predictors, output_path, page, plots, predictor, predictors, start, version_number)
rm(robust, process, dt_robust, dt_standard, dt_minmax, dt, data_list)

# Scale : standardization method  ----
# We choose the standardization method
# Select the columns to scale
num <- pred_tr %>% st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

str(pred_tr[, num])

# Replace the original columns with the scaled values
pred_tr[, num] <- scale(pred_tr[, num] %>% st_drop_geometry())


#--- Export transformed predictors -------------------------------
# predictors_tr_v1.0 ----
st_write(pred_tr, "./data/processed_data/predictors/predictors_tr_v1.0.gpkg", append = FALSE)





























