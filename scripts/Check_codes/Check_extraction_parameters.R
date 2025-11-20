#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(dplyr)
library(exactextractr)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(terra)
library(stringr)



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- Load and prep data ------------------
# Load mtdt_5 
buff <- st_read("./data/processed_data/Mtdt/mtdt_5.gpkg")

colnames(buff)

# Keep only replicates and geom 
buff <- buff %>%
  dplyr::select(replicates, geom)


###### METHODOLIGAL CHECKS ####

# Parameters
var <- "gravity"
rast <- terra::rast("./data/raw_data/predictors/Gravity/rastGravity.tif")
poly <- buff 


# terra::extract weights=TRUE vs exact=TRUE ----
spatial_extraction <- function(var, rast, poly, stats = c("mean", "min", "max", "range")) {
  
  # Load necessary libraries
  require(terra)
  require(sf)
  require(dplyr)
  
  # Ensure stats contains only allowed values
  stats <- intersect(stats, c("mean", "min", "max", "range"))
  
  # Reproject polygons if CRS differs
  raster_crs <- sf::st_crs(terra::crs(rast))
  poly_crs <- sf::st_crs(poly)
  
  if (!identical(poly_crs, raster_crs)) {
    poly <- sf::st_transform(poly, crs = raster_crs)
  }
  
  # Initialize vectors to store the statistics
  rast_means <- numeric(length = nrow(poly))
  rast_mins <- numeric(length = nrow(poly))
  rast_maxs <- numeric(length = nrow(poly))
  rast_ranges <- numeric(length = nrow(poly))
  
  for (i in 1:nrow(poly)) {
    
    # Extract raster values within the current polygon with weights
    extracted_values <- terra::extract(x = rast, y = poly[i, ], weights = TRUE) # The documentation says that with "weights", "the approximate fraction of each cell" is used whereas with "exact", "the exact fraction" is used. The reason for having both is in part because the argument "weights" predates the argument "exact". "weights" was kept because it could be faster and close enough in most cases.
    
    # Filter out NA values
    extracted_values <- extracted_values[!is.na(extracted_values[, 2]), ]
    
    values <- extracted_values[, 2]
    weights <- extracted_values[, 3]
    
    if (length(values) > 0) {
      if ("mean" %in% stats) {
        rast_means[i] <- sum(values * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
      }
      if ("min" %in% stats) {
        rast_mins[i] <- min(values, na.rm = TRUE)
      }
      if ("max" %in% stats) {
        rast_maxs[i] <- max(values, na.rm = TRUE)
      }
      if ("range" %in% stats) {
        rast_ranges[i] <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
      }
    } else {
      if ("mean" %in% stats) rast_means[i] <- NA
      if ("min" %in% stats)  rast_mins[i]  <- NA
      if ("max" %in% stats)  rast_maxs[i]  <- NA
      if ("range" %in% stats) rast_ranges[i] <- NA
    }
  }
  
  # Append results to poly
  if ("mean" %in% stats) poly[[paste0(var, "_mean")]] <- rast_means
  if ("min" %in% stats)  poly[[paste0(var, "_min")]]  <- rast_mins
  if ("max" %in% stats)  poly[[paste0(var, "_max")]]  <- rast_maxs
  if ("range" %in% stats) poly[[paste0(var, "_range")]] <- rast_ranges
  
  return(poly)
}
spatial_extraction_2 <- function(var, rast, poly, stats = c("mean", "min", "max", "range")) {
  
  # Load necessary libraries
  require(terra)
  require(sf)
  require(dplyr)
  
  # Ensure stats contains only allowed values
  stats <- intersect(stats, c("mean", "min", "max", "range"))
  
  # Reproject polygons if CRS differs
  raster_crs <- sf::st_crs(terra::crs(rast))
  poly_crs <- sf::st_crs(poly)
  
  if (!identical(poly_crs, raster_crs)) {
    poly <- sf::st_transform(poly, crs = raster_crs)
  }
  
  # Initialize vectors to store the statistics
  rast_means <- numeric(length = nrow(poly))
  rast_mins <- numeric(length = nrow(poly))
  rast_maxs <- numeric(length = nrow(poly))
  rast_ranges <- numeric(length = nrow(poly))
  
  for (i in 1:nrow(poly)) {
    
    # Extract raster values within the current polygon with exact weights
    extracted_values <- terra::extract(x = rast, y = poly[i, ], exact = TRUE) 
    
    # Filter out NA values
    extracted_values <- extracted_values[!is.na(extracted_values[, 2]), ]
    
    values <- extracted_values[, 2]
    weights <- extracted_values[, 3]
    
    if (length(values) > 0) {
      if ("mean" %in% stats) {
        rast_means[i] <- sum(values * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
      }
      if ("min" %in% stats) {
        rast_mins[i] <- min(values, na.rm = TRUE)
      }
      if ("max" %in% stats) {
        rast_maxs[i] <- max(values, na.rm = TRUE)
      }
      if ("range" %in% stats) {
        rast_ranges[i] <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
      }
    } else {
      if ("mean" %in% stats) rast_means[i] <- NA
      if ("min" %in% stats)  rast_mins[i]  <- NA
      if ("max" %in% stats)  rast_maxs[i]  <- NA
      if ("range" %in% stats) rast_ranges[i] <- NA
    }
  }
  
  # Append results to poly
  if ("mean" %in% stats) poly[[paste0(var, "_mean")]] <- rast_means
  if ("min" %in% stats)  poly[[paste0(var, "_min")]]  <- rast_mins
  if ("max" %in% stats)  poly[[paste0(var, "_max")]]  <- rast_maxs
  if ("range" %in% stats) poly[[paste0(var, "_range")]] <- rast_ranges
  
  return(poly)
}






# Extraction
buff <- spatial_extraction(var = var, rast = rast, poly = poly, stat = c("min", "max", "mean", "range"))
buff2 <- spatial_extraction_2(var = var, rast = rast, poly = poly, stat = c("min", "max", "mean", "range"))







# make correlations of c("gravity_mean", "gravity_min","gravity_max","gravity_range") in buff vs in buff2 
cor(buff$gravity_mean, buff2$gravity_mean, use = "complete.obs", method = "pearson")
cor(buff$gravity_min, buff2$gravity_min, use = "complete.obs", method = "pearson")
cor(buff$gravity_max, buff2$gravity_max, use = "complete.obs", method = "pearson")
cor(buff$gravity_range, buff2$gravity_range, use = "complete.obs", method = "pearson")

# Conclusion : 1, 1, 0.9999994, 0.9999989. We change to "exact" because also very quick. 





# terra::extract exact=TRUE vs exactextractr::exactextract ####

# Reproject polygons if CRS differs
raster_crs <- sf::st_crs(terra::crs(rast))
poly_crs <- sf::st_crs(poly)

ee <- exactextractr::exact_extract(x = rast, y = poly) 
te <- terra::extract(x = rast, y = poly, exact = TRUE)

# ee légèrement plus rapide 

# bind ee and ted 
df <- merge(
  do.call(rbind, lapply(seq_along(ee), \(i) transform(ee[[i]], ID = i, row = seq_len(nrow(ee[[i]]))))),
  transform(te, row = ave(ID, ID, FUN = seq_along)),
  by = c("ID", "row"),
  all = TRUE
)[, c("ID", "value", "coverage_fraction", "layer", "fraction")]

names(df) <- c("ID", "value_ee", "fraction_ee", "value_te", "fraction_te")



# correlation df$value_te vs df$value_ee 
cor(df$value_te, df$value_ee, use = "complete.obs", method = "pearson")

# correlation df$fraction_te vs df$fraction_ee
cor(df$fraction_te, df$fraction_ee, use = "complete.obs", method = "pearson")


# Conclusion : 1, 1. We keep exact=TRUE because also very quick.










#----- CHL exactextract vs intersection extract ----
# Compare exactextract vs rio.clip data ----
chl_rio <- read.csv("./data/processed_data/predictors/mtdt_5_CHL.csv")
chl_ee <- read.csv("./data/processed_data/predictors/mtdt_5_CHL_exactextract.csv")
colnames(chl_rio)
colnames(chl_ee)

common_cols <- intersect(colnames(chl_rio), colnames(chl_ee))[-1]  # exclude "replicates"

for (col in common_cols) {
  correlation <- cor(chl_rio[[col]], chl_ee[[col]], use = "complete.obs", method = "pearson")
  cat("Correlation between rio.clip and exactextractr for", col, ":", correlation, "\n")
}

# Results :
# corr : 0.99077330.98144640.98714470.99423190.98154780.98984130.99572940.97133930.98447660.99740910.91652120.98090670.99691270.92759980.9785911
# Conclusion : very high correlation. We can keep exactextract extraction. 


rm(chl_ee, chl_rio, common_cols, col, correlation)
