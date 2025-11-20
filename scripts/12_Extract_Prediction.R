#------------- Description ----
# The aim of this script is to extract predictor values for the predictions cells. 





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
library(pMEM)



# Load functions
source("./utils/Fct_Data-Prep.R")

#--------------------------------------------------- PREDICTION GRID V.1.0  - JUNE 2023 - ALL PREDICTORS ----------------
#------------- Load and prep data ------------------
# grid_v1.0_20230701 ----
buff <- st_read("./data/processed_data/prediction_extent/grid_v1.0_20230701.gpkg") %>% dplyr::select(-c("value"))


#------------- VECTOR DATA ----------------
#--- Sampling effort : Buffer area ############
# Project in CRS=2154
buff <- st_transform(buff, crs = 2154)

# Compute area in km2
buff$area_km2 <- st_area(buff) / 1e6  # Convert from m^2 to km^2
buff$area_km2 <- units::set_units(st_area(buff), km^2)
buff$area_km2 <- as.numeric(buff$area_km2)  # Convert to numeric for easier handling



#--- Reserve : in or out (Laure) [ TO DO ] #############
#--- Distances : to port, canyon and reserve (Martin) ############
# Load data 
dist <- st_read("./data/raw_data/predictors/Prediction_grid_v1.0/patch_distance_canyon_port_mpa_predict.gpkg") %>% as.data.frame() %>% dplyr::select(-gid)

# Merge 
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(dist), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(dist %>% dplyr::select(id, all_of(extra_cols)), by = "id")

rm(dist, extra_cols)

# Count NA per column
sapply(buff, function(x) sum(is.na(x)))



















#------------- RASTER DATA ----------------

#--- Gravity (1 km) : weighted mean, min, max, range ####
# Parameters
var <- "gravity"
rast <- terra::rast("./data/raw_data/predictors/Gravity/rastGravity.tif")
poly <- buff 

# Extraction
gr <- spatial_extraction(var = var, rast = rast, poly = poly, stat = c("min", "max", "mean", "range"))
beepr::beep()

# Merge
# Detect columns in gr that are not in buff
extra_cols <- setdiff(names(gr), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    gr %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )


# Count NA per column
sapply(buff, function(x) sum(is.na(x)))






rm(gr, extra_cols) # 23 gravity NA 



















#--- Bathymetry (0.001°) ####
# Data : we use the SHOM MNT at 100m resolution

# Parameters
var <- "bathy"
rast <- terra::rast("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")
poly <- buff

# Set values >= to 0 to NA (because it's land)
rast[rast >= 0] <- NA

# Extraction
bat <- spatial_extraction(var = var, rast = rast, poly = poly)
beepr::beep()

# Clean 
bat <- bat %>%
  dplyr::select(c("id", "bathy_min", "bathy_max", "bathy_mean", "bathy_range", ))

bat <- bat %>%
  mutate(bathy_mean = abs(bathy_mean),
         bathy_min2  = abs(bathy_max),
         bathy_max2  = abs(bathy_min)) 

bat <- bat %>%
  dplyr::select(-c("bathy_min", "bathy_max")) %>%
  dplyr::rename(
    bathy_min = bathy_min2,
    bathy_max = bathy_max2
  )





# Merge -----
# Detect columns in bat that are not in buff
extra_cols <- setdiff(names(bat), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    bat %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )




rm(bat, extra_cols, rast, poly, var)

#--- Terrain index ####
# Slope and aspect are computed with 8 neighbors.
# TPI, TRI and roughness are automatically fixed to 8 neighbors.

# Load terrain indices
filelist_temp <- list.files("./data/processed_data/predictors/terrain_indices/SHOM_100m/", pattern = ".tif", full.names = TRUE)

rast <- terra::rast(filelist_temp)

# Iterate through each habitat layer in the raster

# Initialize accumulator
buff_all <- buff 

# Loop through each raster layer
for (layer in names(rast)) {
  tmp <- spatial_extraction(
    var = layer,
    rast = rast[[layer]],
    poly = buff_all,
    stats = c("min", "max", "mean")
  )
  
  # Extract only the columns: replicates + {layer}_min, _max, _mean
  tmp <- tmp %>%
    st_drop_geometry() %>%
    dplyr::select(id,
                  all_of(paste0(layer, "_min")),
                  all_of(paste0(layer, "_max")),
                  all_of(paste0(layer, "_mean")))
  
  # Join to original buff_all
  buff_all <- left_join(buff_all, tmp, by = "id")
}
beepr::beep()

# Merge 
# Detect columns in buff_all that are not in buff
extra_cols <- setdiff(names(buff_all), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    buff_all %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )

sapply(buff, function(x) sum(is.na(x)))
rm(filelist_temp, tmp, layer, poly, rast, i, var)
rm(buff_all, extra_cols)


######  Habitat (100m) ####
# Explanation : we retrieve 3 different variables : 
# 1. main habitat within buffer, 
# 2. number of habitats per km² = number of different habitats / buffer area, 
# 3. surface proportion of each habitat within buffer 
# Data :
# A 100m resolution raster file from Andromed with 14 bands representing the surface of different habitats in m².
# Load raster 
rast <- terra::rast("./data/raw_data/predictors/Habitat/bioc_2023_medfr_100m_2154.tif")

# Create 1 layer of presence/absence for each habitat
rast <- rast %>% as.factor() %>% segregate()

# Rename bands
names(rast) <-  c("Association rhodolithes", "matte_morte_P.oceanica",
                  "Coralligene","P.oceanica","Roche du large", "Algues infralittorales",
                  "Galets infralittoraux","fonds meubles circalittoraux",
                  "Fonds meubles infralittoraux", "Habitats artificiels",
                  "Herbiers Cymodocess", "Zone bathyale","Herbier mixte", "Z.noltei")

# Remove "zone bathyale" and "Cymodocees" from raster
rast <- rast[[!(names(rast) %in% c("Zone bathyale", "Herbiers Cymodocess"))]]

# Group habitats 
# In rast, habitats are separatly named as a function of their location (infralittoral, circalittoral, large) and type (rocky, soft bottom, seagrass, algae, coralligenous) --> we group habitats to keep only types. 

rast_grouped <- terra::rast("./data/processed_data/predictors/Habitat/grouped_habitat_raster.tif")












# Main habitat and Number of habitat / km² 
# Important methodological considerations : "Zone bathyale" and "Cymodocees" are not considered as a habitat in this analysis, thus they are ignored when tehy appear to be the main habitat (we take the next more important habitat) and when counting the number of habitats per km².


# For all habitats
buff1 <- calculate_habitats(buff = buff, 
                            rast = rast, 
                            id_column_name = "id",
                            buff_area = "area_km2")


# For grouped habitats
buff1 <- calculate_habitats(buff = buff1, 
                            rast = rast_grouped, 
                            id_column_name = "id", 
                            name = "grouped_",
                            buff_area = "area_km2")
beepr::beep()








# Merge 
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(buff1), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    buff1 %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )

rm(buff1, extra_cols)





# export temp file
st_write(buff, "./data/processed_data/prediction_extent/grid_v1.0_20230701_with_predictors_temp.gpkg", delete_dsn = TRUE)

# check NA
sapply(buff, function(x) sum(is.na(x))) # 15 NA in main_habitat but not in nb of habitat ?? 


















# Surface proportion of each grouped habitat ----

# Make buff a spatial object
buff <- sf::st_as_sf(buff)

# Initialize df
df <- buff

for (layer in names(rast_grouped)) {
  
  temp <- spatial_extraction(
    var   = layer,
    rast  = rast_grouped[[layer]],
    poly  = df,               
    stats = c("mean")
  )
  
  # Find newly created column names
  new_cols <- setdiff(names(temp), names(df))
  
  if (length(new_cols) > 0) {
    # Drop geometry from temp and bind only the new attribute columns
    df <- bind_cols(df, sf::st_drop_geometry(temp)[, new_cols, drop = FALSE])
  } else {
    message(sprintf("No new columns returned for layer '%'", layer))
  }
}


# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(df), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    df %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(df, extra_cols)

























