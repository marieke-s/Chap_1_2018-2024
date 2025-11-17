#------------- Description ---------------------
# Purpose: 
# This script aims to compute select predictors. 

# Author: Marieke Schultz

# Date script created: 16/11/2025

#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(corrplot)
library(dplyr)
library(exactextractr)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(terra)
library(stringr)
library(pMEM)
library(sf)



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- Load data ------------------
# Load mtdt_7
mtdt <- read_sf("./data/processed_data/Mtdt/mtdt_7.gpkg")

# Load bathy 50m (shom - 100m resolution)
bathy_50m <- terra::rast("./data/processed_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged-50-0.tif")

# Load div_indices_v1.0
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))



#--------------- SITE SELECTION v.1.0 ---------------- 

# Remove buffers that do not overlap with bathy_50m -----

# Use terra::extract to extract bathymetry values --> if NA --> remove polygon

# Extraction
bat <- spatial_extraction(var = "50mbathy", 
                          rast = bathy_50m, 
                          poly = mtdt, 
                          stats = c("max"))

sum(is.na(bat$"50mbathy_max")) # 7 

removed_polys <- bat %>%
  filter(is.na(`50mbathy_max`)) 

# Plot removed_polys on top of bathy_50m ----
plot(bathy_50m)
plot(st_geometry(removed_polys), add = TRUE, col = "red")

# Remove polygons from mtdt ----
mtdt <- mtdt %>%
  filter(!replicates %in% removed_polys$replicates)

div <- div %>%
  filter(!replicates %in% removed_polys$replicates)


# Export mtdt_7_sel_v1.0.gpkg ----
st_write(mtdt, 
         "./data/processed_data/Mtdt/mtdt_7_sel_v1.0.gpkg", 
         delete_dsn = TRUE)


# Export div_indices_v1.0_sel_v1.0.csv ----
# Based on div_indices_v1.0
# write csv 
write_csv2(div,"./data/processed_data/Traits/div_indices_v1.0_sel_v1.0.csv")

#--------------- AUTRES ---------------- 
# remove lockdown

# volume | method | nb of replicates | depth sampling / seafloor

# Selection based on depth  -----
# fish communities differ between 20 and 40meter : Rotanski et al. 2024
# eDNA detection depends on water stratification (high in summer)
# previous decision taken : depth limit set at 30/35 meters

