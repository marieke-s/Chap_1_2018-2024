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

#--------------------------------------------------- PREDICTION GRID V.1.0 - ALL PREDICTORS ----------------
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



#--- Reserve : in or out (Laure) [ NOT DONE ] #############
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



















#--- Gravity : handle NAs ------
#--- Data prep

# Select buff with NA in gravity columns
gravity_cols <- c("gravity_min", "gravity_max", "gravity_mean", "gravity_range")
buff_na <- buff %>%
  filter(if_any(all_of(gravity_cols), is.na)) # 23 NA

rast <- raster::raster("./data/raw_data/predictors/Gravity/rastGravity.tif")

# Set buff_na crs to match rast
buff_na <- st_transform(buff_na, crs = crs(rast))

#--- Extract min, max and mean on 3 nearest pixels
buff_na <- nearest_value(var = "gravity", rast = rast, df = buff_na, poly = buff_na, k = 3)
beepr::beep()

# Merge back NA values -----
# Match rows between buff and buff_na by id
idx <- match(buff$id, buff_na$id)

# Which rows in buff have a matching row in buff_na?
has_match <- !is.na(idx)

# For each gravity column, fill NAs in buff with values from buff_na
for (col in gravity_cols) {
  # rows where buff has NA in this column and there is a corresponding buff_na row
  na_in_buff <- has_match & is.na(buff[[col]])
  
  # replace only those NAs
  buff[[col]][na_in_buff] <- buff_na[[col]][idx[na_in_buff]]
}

# Compute gravity range : 
buff <- buff %>%
  mutate(gravity_range = gravity_max - gravity_min)

# Check NA
sapply(buff, function(x) sum(is.na(x))) # 0 NA now


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





# Merge 
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
rast <- terra::rast("./data/processed_data/predictors/Habitat/habitat_raster.tif")

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






#--- Handle NAs  ------
# Data prep
# For habitat with NA --> take the 3 closest pixel values

# Select buff with NA in habitat columns
habitat_cols <- c("main_habitat", "grouped_main_habitat", "grouped_nb_habitat_per_km2", "nb_habitat_per_km2")
buff_na <- buff %>%
  filter(if_any(all_of(habitat_cols), is.na)) # 15 NA 

buff_na <- buff_na %>% 
  dplyr::select(id, area_km2)






# Compute main habitat and number of habitats per km2 for buff_na grouped and not grouped 

# For grouped habitats
buff1 <- calculate_habitats_nn(
  buff           = buff_na,
  rast           = rast_grouped,
  id_column_name = "id",
  buff_area      = "area_km2",
  name           = "grouped_"      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)

# For all habitats
buff2 <- calculate_habitats_nn(
  buff           = buff_na,
  rast           = rast,
  id_column_name = "id",
  buff_area      = "area_km2",
  name           = ""      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)

beepr::beep()

# Merge back NA values -----

# Merge buff1 and buff 2 by id 
buff12 <- buff1 %>%
  dplyr::select(-area_km2) %>%
  left_join(
    buff2 %>% st_drop_geometry() %>% dplyr::select(-area_km2),
    by = "id"
  )

# For replicates in buff12, replace the hab_cols columns in buff with those in buff12
hab_cols <- colnames(buff12[, -1] %>% sf::st_drop_geometry())
hab_cols <- intersect(hab_cols, names(buff))

# match buff$replicates to buff12$replicates
idx <- match(buff$id, buff12$id)

# which rows in buff have a match in buff12?
rows_to_update <- !is.na(idx)

# replace values in buff[rows_to_update, hab_cols] 
# with the corresponding rows from buff12
buff[rows_to_update, hab_cols] <- buff12[idx[rows_to_update], hab_cols] %>% sf::st_drop_geometry()

# Check NA ----
sapply(buff, function(x) sum(is.na(x))) # 2 NA remaining for 'main habitat'

# Clean ups
rm(buff1, buff2, buff12, buff_na, hab_cols, idx, rows_to_update)









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
rm(rast, rast_grouped)
rm(h, habitats, new_cols)



#--- Surface proportion of each grouped habitat ----
# Handle habitats proportions (grouped only)
# From 2a_Prep_Predictors.R:


# Habitats (base names, must match layers in rast_grouped)
habitats <- c(
  "Algues.infralittorales",
  "matte_morte_P.oceanica",
  "coralligenous",
  "rock",
  "soft_bottom",
  "meadow",
  "Habitats.artificiels"
)

mean_cols <- paste0(habitats, "_mean")

## 1. Remove old habitat mean columns if present
buff <- buff %>%
  dplyr::select(-any_of(mean_cols))

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
beepr::beep()

names(df)
sum(is.na(df$"Habitats artificiels_mean")) # 13 NA

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
beepr::beep()


sum(is.na(df$"Habitats artificiels_mean")) # 0 NA

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
summary(rowSums(df[, mean_cols] %>% st_drop_geometry(), na.rm = TRUE)) # All 1 = good -> nop min = 0
boxplot(rowSums(df[, mean_cols] %>% st_drop_geometry(), na.rm = TRUE))


# Print row with rowSums = 0
df %>%
  filter(rowSums(df[, mean_cols] %>% st_drop_geometry(), na.rm = TRUE) == 0) %>%
  dplyr::select(id, all_of(mean_cols)) %>%
  print(n = Inf)

# The 2 values = 0 correspond to the 2 remaining NAs (bathyale habitat) -> we leave it to NAs since we will modify the habitat variable to include bathyale zone later on. 

## 5. Merge back to buff 
buff <- df %>% dplyr::select(-overlap)


# Check NA
sapply(buff, function(x) sum(is.na(x))) # 2 NA remaining for habitat variables

## Summary of what we did : 

# | Case                              | Habitat value source                       | Final proportions      |
#   | --------------------------------- | ------------------------------------------ | ---------------------- |
#   | **Polygon overlaps raster**       | Weighted means from `spatial_extraction()` | mean_i / sum(mean_all) |
#   | **Polygon outside raster extent** | Sum of 3 nearest pixels per habitat        | sum_i / sum(sum_all)   |
#   






































#------------- VECTOR DATA ----------------
#--- Distance Distance between seabed and depth sampling [ NOT DONE ] -----
#--- Centroid coordinates -----
# Prep coords from buff centroids
cent <- sf::st_centroid(buff) %>%
  st_transform(crs = 4326) 

xy <- st_coordinates(cent)

coords <- buff %>%
  mutate(x = xy[, 1],
         y = xy[, 2])




# Merge ----
extra_cols <- setdiff(names(coords), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    coords %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )

rm(coords, extra_cols, cent, xy)


#==========================
#===== export temp file -----
st_write(buff, "./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_20230701_with_predictors_temp.gpkg", delete_dsn = TRUE)

#===== load temp file -----
buff <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_20230701_with_predictors_temp.gpkg")

# check NA
sapply(buff, function(x) sum(is.na(x))) # 15 NA in main_habitat but not in nb of habitat ?? 

#==========================


























#------------- NCDF DATA ----------------
#--- MARS3D (1.2km) -----
# MARS3D data was extracted by Pauline using scripts "...". 

# Load MARS3D and combine data

#--- Chlorophyll (1km) ----
# The extraction is made on the ncdf "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc" with the script 2.2_Extract_NCDF.ipynb in python.

# Load CHL extracted data (.geojson)
chl <- read.csv("./data/processed_data/predictors/CHL/mtdt_5_CHL-FULL.csv")

# Check chl values ----
colnames(chl)

# Hist
for (c in colnames(chl)[-1]) {
  hist(chl[[c]], main = paste("Histogram of", c), xlab = c, breaks = 50)
}

# NA
for (c in colnames(chl)[-1]) {
  na_count <- sum(is.na(chl[[c]]))
  cat("Number of NA in", c, ":", na_count, "\n")
}


# Merge buff and chl by replicates -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(chl), names(buff))

# left_join 
buff <- buff %>%
  left_join(
    chl %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(chl, extra_cols, c, na_count)











#--- SST (1km) ----
# SST data is extracted from Copernicus data :see script 2.2_Extract_NCDF.ipynb for more details.

# Load SST extracted data (.geojson)
sst <- read.csv("./data/processed_data/predictors/SST/mtdt_5_SST-FULL.csv")

# Check sst values ----
colnames(sst)

# Hist
for (c in colnames(sst)[-1]) {
  hist(sst[[c]], main = paste("Histogram of", c), xlab = c, breaks = 50)
}

# NA
for (c in colnames(sst)[-1]) {
  na_count <- sum(is.na(sst[[c]]))
  cat("Number of NA in", c, ":", na_count, "\n")
}


# Merge buff and sst by replicates -----
# Detect columns that are not in buff
extra_cols <- setdiff(names(sst), names(buff))

# left_join 
buff <- buff %>%
  left_join(
    sst %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(sst, extra_cols, c, na_count)






























