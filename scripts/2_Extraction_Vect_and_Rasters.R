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
############## Sampling effort ############
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



# Ensure area_km2 is numeric (not units) for plotting
buff$area_km2_num <- as.numeric(units::set_units(buff$area_km2, km^2))

# filter if needed
dt <- buff %>%
  filter(area_km2_num > 5)  # Adjust threshold as needed

ggplot(dt) +
  geom_sf(aes(fill = area_km2_num), color = NA) +
  viridis::scale_fill_viridis(name = "Area (km²)", option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "Buffer Area Map", fill = "km²")

# square root of area
sqrt(7)
0.5*0.5
