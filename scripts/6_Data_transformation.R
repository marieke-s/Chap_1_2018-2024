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


#--- Aspect and slope [TODO] ----
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
  st_as_sf() %>%  # Convert back to sf 
  dplyr::select(-"...145")

rm(pred_num, pred_num_log, i, max_val, median_val, new_col_name)


#--- Scale [TODO] -------------------------------
# Explanation
# Why scale and normalize variables?
# First, they can improve the performance and accuracy of many machine learning algorithms, such as linear regression, logistic regression, k-means clustering, and principal component analysis. These algorithms rely on the assumption that the features have similar ranges and distributions, and are sensitive to large differences in scale or magnitude. 
# Second, they can make the data more interpretable and comparable, by eliminating the influence of units or scales. For example, if you want to compare the heights and weights of different people, you need to scale them to a common unit, such as meters and kilograms. 
# Third, they can reduce the risk of numerical errors or overflow, by limiting the range of values and avoiding extremely large or small numbers.

# source : https://www.linkedin.com/advice/3/what-most-effective-techniques-scaling-normalizing-57jlc


















# CODES FROM  CHAPTER 1 -----------------
#----- Log predictors -------------
# Explanation : We log the predictors with high variance to reduce the impact of outliers. We log Fishing Effort and Vessel presence because of that, see their distribution :
hist(pred_pooled$Fishing_Eff, breaks = 100)
hist(pred_pooled$Vessel_Presence_mean, breaks = 100)

# Explanation : We log bathymetry because the effect of pressure are exponential with depth. 
# Log bathymetry
pred_pooled$bathy_min_log <- log(abs(pred_pooled$bathy_min + 1))
pred_pooled$bathy_max_log <- log(abs(pred_pooled$bathy_max + 1))
pred_pooled$bathy_mean_log <- log(abs(pred_pooled$bathy_mean + 1))

# Log Fishing Effort
pred_pooled$Fishing_Eff_log <- log(abs(pred_pooled$Fishing_Eff + 1))

# Log Vessel Presence
pred_pooled$Vessel_Presence_mean_log <- log(abs(pred_pooled$Vessel_Presence_mean + 1))
pred_pooled$Vessel_Presence_min_log <- log(abs(pred_pooled$Vessel_Presence_min + 1))
pred_pooled$Vessel_Presence_max_log <- log(abs(pred_pooled$Vessel_Presence_max + 1))

# Visually check the log transformation
hist(pred_pooled$Fishing_Eff_log, breaks = 100)
hist(pred_pooled$Vessel_Presence_mean_log, breaks = 100)
hist(pred_pooled$bathy_mean_log, breaks = 100)

# Remove unloged columns
pred_pooled <- pred_pooled %>% dplyr::select(-c("bathy_min", "bathy_max", "bathy_mean", "Fishing_Eff", "Vessel_Presence_mean", "Vessel_Presence_min", "Vessel_Presence_max"))

# Rename bathy columns (switch _min and _max)
pred_pooled <- pred_pooled %>% dplyr::rename(bathy_min = bathy_min_log, bathy_max = bathy_max_log, bathy_mean = bathy_mean_log)

# Add back the logs to colnames of bathy
colnames(pred_pooled)[grep("bathy", colnames(pred_pooled))] <- paste0(colnames(pred_pooled)[grep("bathy", colnames(pred_pooled))], "_log")

#----- Scale predictors -----------
# Explanation : we scale our predictor data using the standardization method : Standardization rescales values to have a mean of 0 and a standard deviation of 1, by subtracting the mean and dividing by the standard deviation. 
# Select the columns to scale
num <- pred_pooled %>%
  dplyr::select(-spygen_code, -main_habitat, -date, -reserve, -x, -y, -date)

# Ensure all selected columns are numeric
num <- num %>% mutate(across(everything(), as.numeric))

# Replace the original columns with the scaled values
pred_pooled[, colnames(num)] <- scale(num)


# Visualise distribution of scaled predictors
c <- colnames(num)

for (i in c) {
  # Make a histogram
  hist(pred_pooled[[i]], main = i, breaks = 100)
}

rm(c, i)

str(pred_pooled)














#-------------------------- Scale terrain indices predictors -------------

# Define habitat surface predictors
terrain <- c("TPI_min", "TPI_max", "TPI_mean", "TRI_min", "TRI_max", "TRI_mean",
             "slope_min", "slope_max", "slope_mean", "aspect_min", "aspect_max", "aspect_mean", "roughness_min", "roughness_max", "roughness_mean")

# Visual check of terrain indices distribution
# pred <- dplyr::select_if(predictors[, terrain], is.numeric)
# for (i in colnames(pred)) {
#   print(i)
#   # Make a histogram
#   hist(pred[[i]], main = i, breaks = 100)
# }
# 
# rm(pred, i)

# Check variance-to-mean ratio for each habitat
sapply(terrain, function(h) {
  var_to_mean <- var(tot[[h]]) / mean(tot[[h]])
  cat(h, "Variance/Mean:", var_to_mean, "\n")
})

# Interpreting the variance-to-mean ratio : if the ratio is > 1, the variance is greater than the mean, indicating a high level of dispersion in the data.

#----- Function : Log-transform to reduce outlier impact 
log_transform <- function(df, cols) {
  df[paste0(cols, "_log")] <- log(abs(df[cols]) + 1)
  df <- df %>% dplyr::select(-all_of(cols))  # Remove original columns
  return(df)
}

# Log-transform to reduce outlier impact
tot <- log_transform(tot, terrain)

# Standardize log-transformed predictors (mean = 0, sd = 1)
scaled_terrain <- paste0(terrain, "_log")
tot[, scaled_terrain] <- scale(tot[, scaled_terrain])


# Cleanup
rm(terrain, scaled_terrain, log_transform)



#--------- Save scaled habitat and terrain predictors ---------
# Explanation of file naming : P2 for Processing 2 (from data processing computed in analyses/data_prep/4_Data_Processing.R + adding and scaling habitat and terrain variables + region) _RX for Rarity set to threshold 0 (species with less than X occurences of removed), _IX for the biodiversity index computed, here X = dalongeville for Dalongeville indices.
write_rds(tot, "./data/processed_data/data_prep/Med_TOT_2023_P2_R0_Idalongeville_hababsnolog.rds")


#--------------------------------------- 2. Scale predictors ------------------
# Explanation
# Why scale and normalize variables?
# First, they can improve the performance and accuracy of many machine learning algorithms, such as linear regression, logistic regression, k-means clustering, and principal component analysis. These algorithms rely on the assumption that the features have similar ranges and distributions, and are sensitive to large differences in scale or magnitude. 
# Second, they can make the data more interpretable and comparable, by eliminating the influence of units or scales. For example, if you want to compare the heights and weights of different people, you need to scale them to a common unit, such as meters and kilograms. 
# Third, they can reduce the risk of numerical errors or overflow, by limiting the range of values and avoiding extremely large or small numbers.

# source : https://www.linkedin.com/advice/3/what-most-effective-techniques-scaling-normalizing-57jlc
#---- Load CLEAN_1 datasets ------------------
pred <- read.csv("./data/processed_data/data_prep/predictors/Extracted_Predictors/Med_PRED_full4_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled_CLEAN_1.csv")

mtdt <- read.csv("./data/processed_data/data_prep/eDNA/Subselections/Med_mtdt_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled.csv")
mtdt <- mtdt %>% dplyr::select(-"Unnamed..0")

buff <- st_read("./data/processed_data/data_prep/eDNA/Pts_to_Buffer/Med_mtdt_ADNe_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled_buff_CLEAN_1.shp")
buff <- buff %>% rename(spygen_code = spygn_c) 

#---- Select numerical predictors ----
dt <- pred %>% dplyr::select(-spygen_code, -main_habitat, -date, -reserve, -Index_MPA_num)

# Set all columns to numeric
dt <- as.data.frame(lapply(dt, as.numeric))

#---- Plot predictors distribution ------------------
plot_distributions(dt)




# Observations
# Highly skewed data (mostly to the left) with some outliers to the right.

#---- Scale and plot again ------------------
# Standardization
# Explanation : Standardization rescales values to have a mean of 0 and a standard deviation of 1, by subtracting the mean and dividing by the standard deviation. 
# Convert scaled matrix to data frame
dt_standard <- as.data.frame(scale(dt))
plot_distributions(dt_standard, hist_color = "gold")


# Min Max scale (values between 0 and 1)
# Explanation : Min-max scaling rescales values to a range between 0 and 1, by subtracting the minimum value and dividing by the range. 
process <- preProcess(dt, method=c("range"))
dt_minmax <- predict(process, dt)

plot_distributions(dt_minmax, hist_color = "pink")


# Robust scaling 
# Explanation : Robust scaling rescales values to have a median of 0 and an interquartile range of 1, by subtracting the median and dividing by the interquartile range. The advantage of robust scaling is that it is less sensitive to outliers than min-max scaling or standardization.
robust <- apply(dt, 2, function(x) (x - median(x)) / IQR(x))
dt_robust <- as.data.frame(robust)
rm(robust)
plot_distributions(dt_robust, hist_color = "lightgreen")


# All
# Set up parameters for colors and transparency
colors <- c("lightblue", "gold", "pink", "lightgreen")
alpha <- 0.3

color_transparent <- function(color, alpha) {
  rgb_col <- col2rgb(color)
  rgb(rgb_col[1]/255, rgb_col[2]/255, rgb_col[3]/255, alpha)
}

colors <- sapply(colors, color_transparent, alpha = alpha)

# Combine datasets into a list
data_list <- list(dt, dt_standard, dt_minmax, dt_robust)

# Modify the plotting loop
n_predictors <- ncol(dt)
n_per_page <- 8
n_pages <- ceiling(n_predictors / n_per_page)

for (page in 1:n_pages) {
  start <- (page - 1) * n_per_page + 1
  end <- min(page * n_per_page, n_predictors)
  predictors <- colnames(dt)[start:end]
  
  par(mfrow = c(2, 4)) # Set up a grid layout (2 rows, 4 columns)
  
  for (predictor in predictors) {
    # Plot the histogram for the first dataset
    hist(
      data_list[[1]][[predictor]], 
      main = predictor, 
      xlab = predictor, 
      col = colors[1], 
      border = "white", 
      breaks = 20, 
      freq = FALSE
    )
    
    # Overlay histograms for the other datasets
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
  readline(prompt = "Press [Enter] to see the next set of plots...")
}





#---- Save scaled datasets  ----
# We choose the standardization method
# Select the columns to scale
num <- pred %>%
  dplyr::select(-spygen_code, -main_habitat, -date, -reserve, -Index_MPA_num)

# Ensure all selected columns are numeric
num <- num %>% mutate(across(everything(), as.numeric))

# Replace the original columns with the scaled values
pred[, colnames(num)] <- scale(num)

write.csv(pred, "./data/processed_data/data_prep/predictors/Extracted_Predictors/Med_PRED_full4_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled_CLEAN_1_strd.csv", row.names = FALSE)






# Transform variables with outliers using max value >/= to 10*median value rule --------
# Iterate on all predictors, if max value >/= to 10*median value, apply log(x + 1) transformation, else keep the same
# Create a new data frame to store transformed variables
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
print(paste("Number of log-transformed columns:", sum(grepl("_log$", colnames(pred_num_log)))))

colnames(pred_num_log)
#---- Scale variables ----
pred_num_log_sc <- pred_num_log %>% mutate(across(everything(), as.numeric))