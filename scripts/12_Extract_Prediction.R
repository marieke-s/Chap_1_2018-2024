#------------- Description ----
# The aim of this script is to extract predictor values for the predictions cells. 





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
library(lubridate)
library(terra)
library(stringr)
library(pMEM)
library(purrr)



# Load functions
source("./utils/Fct_Data-Prep.R")

#--------------------------------------------------- PREDICTION GRID V.1.0 - ALL PREDICTORS ----------------
# Without Laure in out reserve variable
# Without features associated to port, canyon and MPA distance.
# With 2 NA remaining in habitat variables
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



#--------------------------------------------------- PREDICTION GRID V.1.0 - JUNE 2023 (2023-07-01) ----------------
#--- MARS3D (1.2km) -----
# Load MARS3D and combine data
files <- list.files("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01_CUR-WIND-SAL-TEMP/", pattern = "*.geojson", full.names = TRUE)
mars_list <- lapply(files, st_read, quiet = TRUE)

# Combine all dataframes into one
mars_df <- do.call(rbind, mars_list)
names(mars_df)

mars_df <- mars_df %>%
  # Add "_month" to each column except "id" and "geometry"
  rename_with(~ paste0(., "_month"), -c("id", "geometry"))


# Merge with buff
extra_cols <- setdiff(names(mars_df), names(buff))
extra_cols <- extra_cols[!extra_cols %in% c("fid_month", "value_month", "date_month", "depth_sampling_surface_month", "depth_sampling_40m_month", "geometry")]

buff <- buff %>%
  left_join(
    mars_df %>% 
      st_drop_geometry() %>% 
      dplyr::select(id, dplyr::all_of(extra_cols)),
    by = "id"
  )

# Check NA
sapply(buff, function(x) sum(is.na(x)))

# Clean up
rm(mars_df, mars_list, files, extra_cols)
rm(t, temp, rast_pts, rast_grouped, nn, df, buff_na, col, col_mean, cols_to_drop, csv_files, gravity_cols, h, habitat_cols, habitats, i, idx, has_match, mean_cols, na_in_buff, new_cols, new_val, no_overlap_ids, vals, total_means, rast, rast_df, all_tables)


#--- CHL (1 km) -----
chl <- read.csv("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01/Cop_CHL_stats_2023-07-01_FULL.csv")

# Merge with buff
buff <- buff %>% 
  left_join(chl, by = "id")

sapply(buff, function(x) sum(is.na(x)))

rm(chl)


#--- SST (1 km) -----
sst <- read.csv("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01/Cop_SST_stats_2023-07-01_FULL.csv")
sapply(sst, function(x) sum(is.na(x)))

# Merge with buff
buff <- buff %>% 
  left_join(sst, by = "id")

rm(sst)




#--- grid_v1.0_2023-07-01_with_predictors-raw_v1.0 ----
# predictors for grid_v1.0
# fixed predictors : all but : 
# Without Laure in out reserve variable
# Without features associated to port, canyon and MPA distance.
# With 2 NA remaining in habitat variables
# climatic predictors : CHL, SST, WIND, VEL, SAL, TEMP for June 2023; SAL TEMP for surface and 40m depth.  

buff %>% 
  dplyr::select(-value) %>%
  st_write("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_v1.0.gpkg", delete_dsn = TRUE)



#--- grid_v1.0_2023-07-01_with_predictors-raw_v1.1 ----
# same as  grid_v1.0_2023-07-01_with_predictors-raw_v1.1 but filling with habitat NA with "soft_bottom" as nearest pixel value is filled by soft bottom. 
buff_filled <- buff %>%
  mutate(
    main_habitat = ifelse(is.na(main_habitat), "soft_bottom", main_habitat),
    grouped_main_habitat = ifelse(is.na(grouped_main_habitat), "soft_bottom", grouped_main_habitat),
    nb_habitat_per_km2 = ifelse(is.na(nb_habitat_per_km2), 1, nb_habitat_per_km2),
    matte_morte_P.oceanica_mean = ifelse(is.na(matte_morte_P.oceanica_mean), 0, matte_morte_P.oceanica_mean),
    `Algues infralittorales_mean` = ifelse(is.na(`Algues infralittorales_mean`), 0, `Algues infralittorales_mean`),
    coralligenous_mean = ifelse(is.na(coralligenous_mean), 0, coralligenous_mean),
    meadow_mean = ifelse(is.na(meadow_mean), 0, meadow_mean),
    `Habitats artificiels_mean` = ifelse(is.na(`Habitats artificiels_mean`), 0, `Habitats artificiels_mean`),
    rock_mean = ifelse(is.na(rock_mean), 0, rock_mean),
    soft_bottom_mean = ifelse(is.na(soft_bottom_mean), 1, soft_bottom_mean)
    
  )
sapply(buff_filled, function(x) sum(is.na(x))) # 0 NA now
buff_filled %>% 
  dplyr::select(-value) %>%
  st_write("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_v1.1.gpkg", delete_dsn = TRUE)




#--- grid_v1.0_2023-07-01_with_predictors-tr_v1.1 ----
pred_raw <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_v1.1.gpkg")


# As in 6_Data-transformation.R : 

# Remove 0 variance predictors ----
nzv <- caret::nearZeroVar(pred_raw %>% st_drop_geometry(), saveMetrics = TRUE)

# Remove ws_min_surface_month 
pred_tr <- pred_raw %>% 
  dplyr::select(-ws_min_surface_month )



# CLT of habitat composition ----

# Remove capital letters in column names
colnames(pred_tr) <- tolower(colnames(pred_tr))

# Replace space and point by "_"
colnames(pred_tr) <- gsub("[ .]", "_", colnames(pred_tr))

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

# Aspect : northness and eastness ----
pred_tr <- pred_tr %>%
  mutate(northness = cos(aspect_mean), 
         eastness = sin(aspect_mean)) %>%
  dplyr::select(-c(aspect_mean, aspect_min, aspect_max))


# Log  ------------------
# Log the same variables than in the training set.

#--- Handle var with negative values ----


# Set the same as for training
var_to_fix <- c(
  "tpi_min",
  "tpi_mean", 
  "tpi_max")

# Function: abs → log(abs+1) → restore sign
log_with_sign <- function(x) {
  abs_x  <- abs(x)
  log_x  <- log(abs_x + 1)
  ifelse(x < 0, -log_x, log_x)
}

# Apply transformation
for (v in var_to_fix) {
  print('----------')
  print(v)
  
  # 0. Take abs value
  a <- abs(pred_tr[[v]])
  
  # 1. Check if variable needs a log (same rule: max ≥ 10 * median)
  if (max(a, na.rm = TRUE) > 10 * median(a, na.rm = TRUE)) {
    message(paste(v, ": log transform applied"))
    
    # 2. Create the log-signed variable in pred_tr
    pred_tr[[paste0(v, "_log")]] <- log_with_sign(pred_tr[[v]])
    
    # 3. Remove original variable
    pred_tr[[v]] <- NULL
  } else {
    message(paste(v, ": no log transform needed"))
  }
}





# Other variables -----
# Log the same predictors than in the training set. 

# load training predictors
p <- st_read("./data/processed_data/predictors/predictors_tr_v1.2.gpkg")

# retain log cols (containing _log)
log_cols <- colnames(p %>% st_drop_geometry() %>% dplyr::select(contains("_log")))
# remove _log suffix to get original names
log_cols <- sub("_log$", "", log_cols)

# Keep cols in pred_tr that are in log_cols
l <- intersect(log_cols, colnames(pred_tr))

# add to l : "vel_min_surface_month"
l <- c(l, "vel_min_surface_month" )

pred_num <- pred_tr %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(l))

colnames(pred_num_log)

pred_num_log <- pred_num  # Initialize new dataframe for transformed columns

for (i in colnames(pred_num)) {
  # Apply log transformation if the condition is met
      print(paste("Transforming:", i))  # Print transformed column name
    new_col_name <- paste0(i, "_log")  # New name with _log suffix
    
    # Replace original column with log-transformed version (renamed)
    pred_num_log <- pred_num_log[, !colnames(pred_num_log) %in% i]  # Drop old column
    pred_num_log[[new_col_name]] <- log(pred_num[[i]] + 1)  # Add new transformed column

}

# Check NA 
sapply(pred_num_log, function(col) sum(is.na(col))) # no NA

# Print the nb of transformed columns
print(paste("Number of log-transformed columns:", sum(grepl("_log$", colnames(pred_num_log))))) # 19 (+2 tpi) 
colnames(pred_num_log)

# Replace original numeric columns in pred_tr with transformed ones
pred_tr <- pred_tr %>%
  st_drop_geometry() %>%
  dplyr::select(-colnames(pred_num)) %>%  # Drop the original numeric columns
  dplyr::bind_cols(pred_num_log) %>%       # Add the transformed numeric columns
  dplyr::bind_cols(st_geometry(pred_tr)) %>% # Reattach geometry
  st_as_sf() # Convert back to sf 

pred_tr <- pred_tr %>% dplyr::select(-"...84")
rm(pred_num, pred_num_log, i, max_val, median_val, new_col_name)


names(pred_tr)


















# Scale ----
# Select the columns to scale
num <- pred_tr %>% st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-c("x", "y"))%>%  # Do not scale coordinates
  colnames()

str(pred_tr[, num])

# Replace the original columns with the scaled values
pred_tr[, num] <- scale(pred_tr[, num] %>% st_drop_geometry())


# Export ----
# based on grid_v1.0_2023-07-01_with_predictors-raw_v1.1
# same transfo as predictors_tr_v1.2
pred_tr %>%  st_write("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-tr_v1.1.gpkg", delete_dsn = TRUE)




#--- grid_v1.0_2023-07-01_with_predictors-sel_v1.3 ----
pred_raw <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_v1.1.gpkg")
# Select predictors as in predictors_sel_v1.3

pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.3.gpkg")

gp <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-tr_v1.1.gpkg")

names(pred)
names(gp)

cols <- intersect(names(pred), names(gp))
setdiff(colnames(pred), cols)
cols <- c(cols, "ws_mean_surface_month", "vel_mean_surface_month", "temp_mean_40m_month", "temp_mean_surface_month", "sal_mean_40m_month", "sal_mean_surface_month")

gp <- gp %>%
  dplyr::select(all_of(c("id", "date", cols, "geom")))
gp <- gp %>%  dplyr::select(-c("date"))
gp <- gp %>%  dplyr::select(-c("date"))

setdiff(names(gp), names(pred))
setdiff(names(pred), names(gp))
sort(names(gp))
sort(names(pred))



# Export ----
st_write(gp, "./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors_sel_v1.3.gpkg" )

#--------------------------------------------------- PREDICTION GRID V.1.0 - FIXED PREDICTORS ----------------
# grid_v1.0_2023-07-01_with_predictors-raw_FIXED_v1.1 ----

buff <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_v1.1.gpkg")
names(buff)

# remove cols containing 'chl' 'temp' 'sal' 'wind' 'vel' 'sst'
keywords <- c("chl", "temp", "sal", "wind", "vel", "sst", "ws")

cols_to_remove <- names(buff)[grepl(paste(keywords, collapse = "|"),
                                    names(buff),
                                    ignore.case = TRUE)]

buff_clean <- buff %>% dplyr::select(-all_of(cols_to_remove))
names(buff_clean)

# Export 
st_write(buff_clean, "./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_FIXED_v1.1.gpkg", delete_dsn = TRUE)



#--------------------------------------------------- PREDICTION GRID V.1.1 - APRIL-SEPTEMBER 2023 ----------------
#--- Load grid_v1.0_2023-07-01_with_predictors-raw_FIXED_v1.1 ----
buff <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/grid_v1.0_2023-07-01_with_predictors-raw_FIXED_v1.1.gpkg")
# CHL -----
# Lister les fichiers CSV
chl_files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_CHL_2018-2024/grid_v1.0",
  pattern = "\\.csv$", 
  full.names = TRUE
)

# Keep only files ending with "FULL.csv"
chl_files <- chl_files[grepl("FULL\\.csv$", chl_files)]

## 1. Extract the date from the file name
chl_file_info <- tibble::tibble(
  path = chl_files,
  date = str_extract(basename(chl_files), "\\d{4}-\\d{2}-\\d{2}")
)

## 2. Group file paths by date
chl_paths_by_date <- split(chl_file_info$path, chl_file_info$date)

## 3. Read files per date (one df per date)
chl_by_date <- lapply(chl_paths_by_date, function(paths) {
  # in case there is >1 file per date in the future, bind rows
  map_dfr(paths, ~ readr::read_csv(.x, show_col_types = FALSE))
})

## 4. Create buff_by_date by joining each CHL df to buff
buff_by_date <- lapply(names(chl_by_date), function(d) {
  chl_df <- chl_by_date[[d]]
  
  # keep only columns that are not already in buff
  extra_cols <- setdiff(names(chl_df), names(buff))
  
  buff %>%
    left_join(
      chl_df %>%
        dplyr::select(id, dplyr::all_of(extra_cols)),
      by = "id"
    ) %>%
    # if any weird X2023 columns appear, drop them (same as CUR-WIND block)
    dplyr::select(-starts_with("X2023"))
})

names(buff_by_date) <- names(chl_by_date)

## Quick check
names(buff_by_date$`2023-10-01`)
dim(buff_by_date$`2023-05-01`)






























# SST -----
# Lister les fichiers CSV
sst_files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_SST_2018-2024/grid_v1.0",
  pattern = "\\.csv$", 
  full.names = TRUE
)

# Keep only files ending with "FULL.csv"
sst_files <- sst_files[grepl("FULL\\.csv$", sst_files)]



## 1. Extract date from file names
sst_file_info <- tibble(
  path = sst_files,
  date = str_extract(basename(sst_files), "\\d{4}-\\d{2}-\\d{2}")
)

## 2. Group file paths by date
sst_paths_by_date <- split(sst_file_info$path, sst_file_info$date)

## 3. Read SST files per date (one df per date)
sst_by_date <- lapply(sst_paths_by_date, function(paths) {
  map_dfr(paths, ~ readr::read_csv(.x, show_col_types = FALSE))
})

## 4. Append SST data into EXISTING buff_by_date
for (d in intersect(names(buff_by_date), names(sst_by_date))) {
  message("Updating buff_by_date for ", d, " with SST data")
  
  sst_df <- sst_by_date[[d]]
  
  # only add columns that are not already present in this date's buff
  extra_cols <- setdiff(names(sst_df), names(buff_by_date[[d]]))
  
  buff_by_date[[d]] <- buff_by_date[[d]] %>%
    left_join(
      sst_df %>%
        dplyr::select(id, dplyr::all_of(extra_cols)),
      by = "id"
    ) %>%
    # optional: same cleanup pattern as before
    dplyr::select(-starts_with("X2023"))
}

# quick sanity check
names(buff_by_date$`2023-05-01`)
dim(buff_by_date$`2023-05-01`)






# CUR-WIND -----

## 1. List all your GeoJSON files
files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_MARS3D_2018-2024/grid_v1.1/CUR-WIND",
  pattern = "*\\.geojson$",
  full.names = TRUE
)

## 2. Extract the date from the file name
file_info <- tibble::tibble(
  path = files,
  date = str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}")
)

## 3. Group file paths by date
paths_by_date <- split(file_info$path, file_info$date)

## 4. Read + merge files so you have one sf per date
# merged_sf_list is a named list, names = "2023-05-01", etc.
merged_sf_list <- lapply(paths_by_date, function(paths) {
  map_dfr(paths, ~ st_read(.x, quiet = TRUE))
})

dim(merged_sf_list$`2023-05-01`) # should be 4423


## 5. Loop over merged files and APPEND to existing buff_by_date
for (d in intersect(names(buff_by_date), names(merged_sf_list))) {
  message("Updating buff_by_date for ", d, " with CUR-WIND data")
  
  file <- merged_sf_list[[d]] %>% st_drop_geometry()
  
  # only add columns that are not already present in this date's buff
  extra_cols <- setdiff(names(file), names(buff_by_date[[d]]))
  
  buff_by_date[[d]] <- buff_by_date[[d]] %>%
    left_join(
      file %>%
        dplyr::select(id, dplyr::all_of(extra_cols)),
      by = "id"
    ) %>%
    # Remove cols starting with 'X2023' if they appear
    dplyr::select(-starts_with("X2023"))
}

# Check results
names(buff_by_date$`2023-05-01`)

# Clean up 
rm(file_info, paths_by_date, merged_sf_list)
























# SAL-TEMP-surf -----

## 1. List all your GeoJSON files
files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_MARS3D_2018-2024/grid_v1.1/SAL-TEMP/surf",
  pattern = "*\\.geojson$",
  full.names = TRUE
)

## 2. Extract the date from the file name
file_info <- tibble::tibble(
  path = files,
  date = str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}")
)

## 3. Group file paths by date
paths_by_date <- split(file_info$path, file_info$date)

## 4. Read + merge files so you have one sf per date
# merged_sf_list is a named list, names = "2023-05-01", etc.
merged_sf_list <- lapply(paths_by_date, function(paths) {
  map_dfr(paths, ~ st_read(.x, quiet = TRUE))
})

dim(merged_sf_list$`2023-05-01`) # should be 4423


## 5. Loop over merged files and APPEND to existing buff_by_date
for (d in intersect(names(buff_by_date), names(merged_sf_list))) {
  message("Updating buff_by_date for ", d, " with CUR-WIND data")
  
  file <- merged_sf_list[[d]] %>% st_drop_geometry()
  
  # only add columns that are not already present in this date's buff
  extra_cols <- setdiff(names(file), names(buff_by_date[[d]]))
  
  buff_by_date[[d]] <- buff_by_date[[d]] %>%
    left_join(
      file %>%
        dplyr::select(id, dplyr::all_of(extra_cols)),
      by = "id"
    ) %>%
    # Remove cols starting with 'X2023' if they appear
    dplyr::select(-starts_with("X2023")) %>%
   # Add '_surf' suffix to SAL and TEMP columns
    rename_with(~ paste0(., "_surf"), 
                .cols = c(starts_with("sal_"), starts_with("temp_")))
}

# Check results
names(buff_by_date$`2023-05-01`)

# Clean up 
rm(file_info, paths_by_date, merged_sf_list)




# SAL-TEMP-bottom -----

## 1. List all your GeoJSON files
files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_MARS3D_2018-2024/grid_v1.1/SAL-TEMP/bottom",
  pattern = "*\\.geojson$",
  full.names = TRUE
)

## 2. Extract the date from the file name
file_info <- tibble::tibble(
  path = files,
  date = str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}")
)

## 3. Group file paths by date
paths_by_date <- split(file_info$path, file_info$date)

## 4. Read + merge files so you have one sf per date
# merged_sf_list is a named list, names = "2023-05-01", etc.
merged_sf_list <- lapply(paths_by_date, function(paths) {
  map_dfr(paths, ~ st_read(.x, quiet = TRUE))
})

dim(merged_sf_list$`2023-05-01`) # should be 4423


## 5. Loop over merged files and APPEND to existing buff_by_date
for (d in intersect(names(buff_by_date), names(merged_sf_list))) {
  message("Updating buff_by_date for ", d, " with CUR-WIND data")
  
  file <- merged_sf_list[[d]] %>% st_drop_geometry()
  
  # only add columns that are not already present in this date's buff
  extra_cols <- setdiff(names(file), names(buff_by_date[[d]]))
  
  buff_by_date[[d]] <- buff_by_date[[d]] %>%
    left_join(
      file %>%
        dplyr::select(id, dplyr::all_of(extra_cols)),
      by = "id"
    ) %>%
    # Remove cols starting with 'X2023' if they appear
    dplyr::select(-starts_with("X2023")) %>%
    # Add '_surf' suffix to SAL and TEMP columns
    rename_with(~ paste0(., "_bottom"), 
                .cols = c(starts_with("sal_"), starts_with("temp_")))
}

# Check results
names(buff_by_date$`2023-05-01`)

# Clean up 
rm(file_info, paths_by_date, merged_sf_list)



