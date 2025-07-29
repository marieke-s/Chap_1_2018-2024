#------------- Description ---------------------
# Purpose: 
# This script aims to compute covariables associated to eDNA replicates. The covariables in this script come from : 
# vector computation 
# extraction from rasters 

# The expected result is a csv file with all original metedata info and a supplementary covariables columns each eDNA replicate.

# Data source: 
# Raster files of predictor (see data section for details)
# A subselection of metadata of the mediterranean eDNA samples grouped by replicates (see 1_Prep_Mtdt.R for details).

# Author: Marieke Schultz

# Contributors : Marie Orblin (script MPA)

# Date script created: 25/07/2025

#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(dplyr)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(terra)
library(stringr)





#------------- Load mtdt ------------------ 
buff <- st_read("./data/processed_data/eDNA/mtdt_5.gpkg")

#------------- Vector data ----------------
######################## Sampling effort ############
# Explanation : the aim is to compute variables that represent the sampling effort in each replicate group : 
# 1. Total volume filtered in the replicate group = estimated_volume_total -->  Already exists
# 2. Nb of PCR replicates = 12 * nb of pooled samples + 12 nb of unpooled samples --> computed here 
# 3. Area covered = area of replicate group buffer --> computed here


# Nb of PCR replicates ----
# Function to compute PCR replicate count
compute_pcr_replicates <- function(rep_string) {
  if (is.na(rep_string) || rep_string == "") return(NA_integer_)
  
  # Split on `/` to get number of samples
  parts <- strsplit(rep_string, "/")[[1]]
  
  # Count each pooled group (with `_`) or single sample as 1 unit
  replicate_units <- sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) == 1)) +  # unpooled
    sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) > 1))     # pooled
  
  return(12 * replicate_units)
}

# Apply to your data
buff$PCR_replicates <- sapply(buff$replicates, compute_pcr_replicates)

# Clean
rm(compute_pcr_replicates)







# Area covered ----
# Project in CRS=2154
buff <- st_transform(buff, crs = 2154)

# Compute area in km2
buff$area_km2 <- st_area(buff) / 1e6  # Convert from m^2 to km^2
buff$area_km2 <- units::set_units(st_area(buff), km^2)


summary(buff$area_km2)
str(buff$area_km2)



# [Optional Plot] ----

# # Ensure area_km2 is numeric (not units) for plotting
# buff$area_km2_num <- as.numeric(units::set_units(buff$area_km2, km^2))
# 
# # filter if needed
# dt <- buff %>%
#   filter(area_km2_num > 5)  # Adjust threshold as needed
# 
# ggplot(dt) +
#   geom_sf(aes(fill = area_km2_num), color = NA) +
#   viridis::scale_fill_viridis(name = "Area (km²)", option = "plasma", direction = -1) +
#   theme_minimal() +
#   labs(title = "Buffer Area Map", fill = "km²")


######################## Distance to shore ############
# Explanation : distance to shore is computed from the coastline shapefile using 'Distance_shore.ipynb'.
# Important methodological considerations : 
# Distance to shore is computed from replicates group's buffer centroids.
# When a centroid is on land (which happened for ~ 67/792 centroids) we set the distance to shore to 1 meter because it was actually not sampled on land (obviously) and it ended up there only because we simplified transects into straight lines between transect start and end points.



# Load data after dist_shore computed  ----
buff <- st_read("./data/processed_data/predictors/mtdt_5_dist-shore.gpkg")
######################## Vessel presence (Luka) ############
######################## MPA : inside-outside reserve ####
# Data :
# The 'protection' variable is : 0 = not in any MPA, 1 = in MPA not fully protected, 2 = in fully protected MPA.
# The methodological method used is : Taking the highest protection level across all the MPA legislation crossed by start and end points of the transects.
# From Marie Orblin and Lola Romant (MARBEC)
# Here we keep only 1 = fully protected areas and 0 = the rest. Then we match this data with our df.

# Load data
mpa <- read.csv("data/raw_data/predictors/MPA/MetaData_AMP.csv")
colnames(mpa)[1] <- "spygen_code"

# Keep only the rows with in df 
mpa <- mpa[mpa$spygen_code %in% mtdt$spygen_code, ] |> 
  dplyr::select(spygen_code, protection)

# Bind mpa$protection to df
df <- df
df <- merge(df, mpa, by = "spygen_code")

# Change protection values 2 --> 1, 1 --> 0
df$reserve <- ifelse(df$protection == 2, 1, 0)

# Remove protection column
df <- df |> dplyr::select(-protection)




######################## MPA : protection level ####
# Data :  
# The 'Index_MPA_num' variable is : 0-5 according to protection levels defined in the MPA Guide (Grorud-Colvert et al. 2021)
# The methodological method used is : Taking the highest protection level across all the MPA legislation crossed by start and end points of the transects.
# From Marie Orblin and Lola Romant (MARBEC)
# Here we match this data with our df. 

# Load data
mpa <- read.csv("data/raw_data/predictors/MPA/MetaData_AMP.csv")
colnames(mpa)[1] <- "spygen_code"

# Keep only the rows with in df 
mpa <- mpa[mpa$spygen_code %in% mtdt$spygen_code, ] |> 
  dplyr::select(spygen_code, Index_MPA_num)

# Bind mpa$Index_MPA_num to df
df <- df
df <- merge(df, mpa, by = "spygen_code")


######################## MPA : distance to reserve ####
# Data :
# The distance to the nearest reserve is computed by Marie Orblin (cf "./analyses/archived/Sent_scripts/in_out_amp24_Marie_Orblin.R". The computed distance is the minimul distance between the centroid of a reserve (fully protected MPA) and the start point of the transect.
# Here we match this data with our df. 


# Load data
mpa <- read.csv("data/raw_data/predictors/MPA/MetaData_AMP.csv")
colnames(mpa)[1] <- "spygen_code"

# df <- read.csv( "./data/processed_data/data_prep/Med_TOT_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled_rare5_coord_P0O0.csv")

# Keep only the rows with in df 
mpa <- mpa[mpa$spygen_code %in% df$spygen_code, ] |> 
  dplyr::select(spygen_code, dist_min_fully_protected_MPA)

# Bind mpa$dist_min_fully_protected_MPA to df
df <- merge(df, mpa, by = "spygen_code")

# Save df
df_3 <- df
write.csv(df, "./data/processed_data/data_prep/Med_TOT_2023_FR_coastal_sd30m_noVH4_Litto3D_balanced2rep_pooled_rare5_coord_dist_MPA_P0O0.csv")





