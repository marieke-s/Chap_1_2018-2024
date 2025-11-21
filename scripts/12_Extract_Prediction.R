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












#--- Main habitat and Number of habitat / km² ----
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





#==========================
#===== export temp file -----
st_write(buff, "./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_20230701_with_predictors_temp.gpkg", delete_dsn = TRUE)

#===== load temp file -----
buff <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_20230701_with_predictors_temp.gpkg")

# check NA
sapply(buff, function(x) sum(is.na(x))) # 15 NA in main_habitat but not in nb of habitat ?? 

#==========================
#--- Surface proportion of each grouped habitat ----
# Handle habitats proportions (grouped only) -----
# From 2a_Prep_Predictors.R:

## 2. Compute weighted means of each habitat within buffers

buff <- sf::st_as_sf(buff)
df   <- buff

# Extract weighted means per habitat
for (h in names(rast_grouped)) {
  temp <- spatial_extraction(
    var   = h,
    rast  = rast_grouped[[h]],
    poly  = df,
    stats = "mean"
  )
  
  # The new column should be like "<h>_mean"
  new_cols <- setdiff(names(temp), names(df))
  
  if (length(new_cols) > 0) {
    df <- dplyr::bind_cols(df, sf::st_drop_geometry(temp)[, new_cols, drop = FALSE])
  } else {
    message(sprintf("No new columns returned for layer '%s'", h))
  }
}

names(df)
sum(is.na(df$"Habitats artificiels_mean")) # 5 NA

# Determine which polygons overlapped the raster
habitats <- paste0(names(rast_grouped), "_mean")
df$overlap <- apply(df[, habitats], 1, function(x) any(is.na(x)))



## 3. For polygons with NO overlap → compute proportions using the 3 nearest pixels

# Convert raster to points
rast_df <- terra::as.data.frame(rast_grouped, xy = TRUE, na.rm = TRUE)
rast_pts <- sf::st_as_sf(rast_df, coords = c("x", "y"), crs = st_crs(df))


# Compute nearest-pixel sums per habitat

no_overlap_ids <- which(df$overlap == TRUE)

for (i in no_overlap_ids) {
  
  # find nearest 3 raster points
  nn  <- nngeo::st_nn(df[i, ], rast_pts, k = 3)
  idx <- unlist(nn)
  
  if (length(idx) == 0) {
    # nothing found – keep NA for this polygon
    next
  }
  
  # for each habitat, sum the values of the 3 nearest pixels
  for (h in names(rast_grouped)) {
    col_mean <- paste0(h, "_mean")
    
    # make sure the column exists and is numeric
    if (!col_mean %in% names(df)) {
      df[[col_mean]] <- NA_real_
    }
    
    vals <- rast_df[idx, h]
    new_val <- if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
    
    # this is the crucial line: assign as a scalar into sf/tibble
    df[i, col_mean] <- new_val
  }
}

# Check results
df %>% filter(overlap == TRUE) %>% dplyr::select(c(habitats)) %>% print()






## 4. Convert these habitat means into TRUE PROPORTIONS

mean_cols <- paste0(names(rast_grouped), "_mean")

total_means <- rowSums(df[, mean_cols] %>% st_drop_geometry(), na.rm = TRUE)


for (col in mean_cols) {
  df[[col]] <- ifelse(
    total_means > 0,
    df[[col]] / total_means,
    NA_real_
  )
}

# Check 
summary(rowSums(df[, mean_cols] %>% st_drop_geometry(), na.rm = TRUE)) # All 1 = good


## 5. Merge back to buff 
buff <- df %>% dplyr::select(-overlap)


## Summary of what we did : 

# | Case                              | Habitat value source                       | Final proportions      |
#   | --------------------------------- | ------------------------------------------ | ---------------------- |
#   | **Polygon overlaps raster**       | Weighted means from `spatial_extraction()` | mean_i / sum(mean_all) |
#   | **Polygon outside raster extent** | Sum of 3 nearest pixels per habitat        | sum_i / sum(sum_all)   |
#   





# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(df), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    df %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )

rm(df, extra_cols)










#=======================================================================================================================
#======================================================= TO DO ===========================================================================================
#=======================================================================================================================
#------------- VECTOR DATA ----------------
###### Distance between seabed and depth sampling ############
# Explanation : bathy_max - depth_sampling

dsamp <- st_read("./data/processed_data/Mtdt/mtdt_5.gpkg") %>% # Load depth_sampling
  st_drop_geometry() %>% 
  dplyr::left_join(buff %>% st_drop_geometry(), by = "replicates") %>% # merge to buff 
  mutate(dist_seabed_depthsampling = bathy_max - depth_sampling)  # compute distance
#dplyr::select(replicates, dist_seabed_depthsampling, depth_sampling, bathy_min, bathy_max, bathy_mean) # keep only necessary columns


# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(dsamp), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    dsamp %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(dsamp, extra_cols)









##### Centroid coordinates ####
# Prep coords from buff centroids
# modif v1.1 : Convert to decimal degree
buff <- sf::st_read("./data/processed_data/predictors/predictors_raw_v1_0.gpkg")

cent <- sf::st_centroid(buff) %>%
  st_transform(crs = 4326) # 

xy <- st_coordinates(cent)

coords <- buff %>%
  mutate(x = xy[, 1],
         y = xy[, 2])


# v1.1 modif 
buff <- buff %>%
  mutate(x = xy[, 1],
         y = xy[, 2])


# Merge ----
extra_cols <- setdiff(names(coords), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    coords %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(coords, extra_cols, cent, xy)



















