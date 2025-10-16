#------------- Description ---------------------
# Purpose: 
# This script aims to compute covariables associated to eDNA replicates. The covariables in this script come from : 
# vector computation 
# extraction from rasters 

# The expected result is a csv file with all original metedata info and a supplementary covariables columns each eDNA replicate.

# Data source: 
# Raster and vector files for predictors (see data section for details)
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



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- VECTOR DATA ----------------

###### Distance to shore (Pauline/Marieke) ############
# Explanation : distance to shore is computed from the coastline shapefile using '2.1_Distance_shore.py' from orignal code of Pauline Viguier. 
# Important methodological considerations : 
# Di(Paulostance to shore is computed from replicates group's buffer centroids.
# When a centroid is on land (which happened for ~ 67/792 centroids) we set the distance to shore to 1 meter because it was actually not sampled on land (obviously) and it ended up there only because we simplified transects into straight lines between transect start and end points.

# Run 2.1_Distance_shore.py 
# Bash call
system("python3 ./scripts/2.1_Distance_shore.py")


# Load buff with dist_shore  
buff <- st_read("./data/processed_data/predictors/mtdt_5_dist-shore.gpkg")






###### Distance to port, canyon and reserve (Martin) ############
# Explanation : distances to several entities (port, canyons, MPA reserve) were computed by Martin Paquet in 07/2025 with the distance pipeline explained in Methods.txt.

# Load Martin's distance 
dist <- st_read("./data/raw_data/predictors/Distances/buffer_with_closest_feats.gpkg")
dist <- as.data.frame(dist)

# Merge buff and dist by replicates
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(dist), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(dist %>% dplyr::select(replicates, all_of(extra_cols)), by = "replicates")

rm(dist, extra_cols)

###### In or out reserve (Laure) #############
# Explanation : variable that states whether or not the replicates were done within a fully protected MPA (i.e. a reserve).

# MPA data used : The MPA data used was made in several steps :
# 1) 2024 : Lola Romant assigned protection levels to the French Mediterranean MPA based on MPA Guide Categories during her internship. # --> protectionMed_fr_modif.gpkg

# 2) This data was then modified by Laure Velez to keep only the MPA with the highest protection level when several were overlapping. 
# --> protectionMed_fr.shp

# 3) 07/2025 : Marieke Schultz manually (QGIS) modified the data to correct errors based on current legislation : 
# "Passe de la Réserve naturelle de Cerbicale (Bouches De Bonifacio)" : from level 3 to 4
# "Passe de l'Archipel des îles de Lavezzi (Bouches De Bonifacio)" : from level 3 to 4
# I created the column ‘Fully’ in which levels 1, 2 and 3 are marked “YES” and the rest ‘NO’.
# --> protectionMed_fr_modif.gpkg

# 4) 07/2025 : Laure Velez assigned to each sample whether or not it belonged to a fully protected MPA, using the protectionMed_fr_modif.gpkg data. She checked for each eDNA samples whether or not they belonged to a fully protected MPA, cheking both start point, end point, and transect (eg when start and end points were outside but transect crossing the reserve). 
# --> mtd_mpafully_2018-2024_complete.csv with column "mpa_fully". 


# Load data
in_out <- read.csv("./data/raw_data/predictors/MPA/mtd_mpafully_2018-2024_complete.csv")
unique(in_out$mpa_fully)


# Here we assign the variable "mpa_fully" to the buff data. 
# When all samples of the replicate group have mpa_fully = 1 we set buff$mpa_fully = 1, when all samples of the replicate group have mpa_fully = 0 we set buff$mpa_fully = 0, else we set buff$mpa_fully = NA.

# Step 1: Create a named vector of spygen_code -> mpa_fully
mpa_lookup <- setNames(in_out$mpa_fully, in_out$spygen_code)

# Step 2: Function to calculate mpa_fully for one buff row
get_mpa_fully_status <- function(replicate_str) {
  # Extract SPY codes from the string (handles both '/' and '_')
  codes <- unlist(strsplit(replicate_str, "[/_]"))
  
  # Lookup corresponding mpa_fully values
  values <- mpa_lookup[codes]
  
  # Remove any missing spygen_codes
  values <- values[!is.na(values)]
  
  if (length(values) == 0) {
    return(NA)
  } else if (all(values == 1)) {
    return(1)
  } else if (all(values == 0)) {
    return(0)
  } else {
    return(NA)
  }
}

# Step 3: Apply to buff
buff$mpa_fully <- vapply(buff$replicates, get_mpa_fully_status, FUN.VALUE = numeric(1))

# Clean
rm(get_mpa_fully_status, mpa_lookup, in_out)

# Check results
buff %>% 
  group_by(mpa_fully) %>%
  summarise(count = n()) %>%
  print()

# print the replicates column of buff where mpa_fully is NA
buff %>% filter(is.na(mpa_fully)) %>% dplyr::select(replicates, mpa_fully) %>%
  print()


# No NA found 

###### Comparison distance to reserve Martin with in_out Laure ############

# Make column reserve_martin : if $mpa_dist_m_min = 0 then reserve_martin = 1 else reserve_martin = 0
buff$reserve_martin <- ifelse(buff$mpa_dist_m_min == 0, 1, 0)

# Compare buff$mpa_fully with buff$reserve_martin
comparison <- buff %>%
  dplyr::select(replicates, mpa_fully, reserve_martin) %>%
  mutate(match = ifelse(mpa_fully == reserve_martin, TRUE, FALSE))

# Count matches and mismatches
comparison %>%
  group_by(match) %>%
  summarise(count = n()) %>%
  print()

61/788 *100 # 7.7% mismatch


# Map the mismatches to see where they are located
mismatches <- buff %>%
  filter(mpa_fully != reserve_martin | is.na(mpa_fully) | is.na(reserve_martin))

# Plot mismatches
plot(st_geometry(mismatches), col = ifelse(mismatches$mpa_fully == 1, "blue", "red"), main = "Mismatches between mpa_fully and reserve_martin")


rm(comparison, mismatches)

# remove buff$reserve_martin
buff <- buff %>% dplyr::select(-reserve_martin)


###### Sampling effort ############
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
buff$area_km2 <- as.numeric(buff$area_km2)  # Convert to numeric for easier handling

######################## Vessel presence (Luka) ############


#------------- RASTER DATA ----------------
###### Gravity (1 km) : weighted mean, min, max, range ####
# Data : From Laure Velez (MARBEC)

# Parameters
var <- "gravity"
rast <- terra::rast("./data/raw_data/predictors/Gravity/rastGravity.tif")
poly <- buff 

# Extraction
buff <- spatial_extraction(var = var, rast = rast, poly = poly, stat = c("min", "max", "mean", "range"))




###### Bathymetry (0.001°) : weighted mean, min, max, range ####
# Data : we use the SHOM MNT at 100m resolution

# Parameters
var <- "bathy"
rast <- terra::rast("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")
poly <- buff

# Set values >= to 0 to NA (because it's land)
rast[rast >= 0] <- NA

# Extraction
buff <- spatial_extraction(var = var, rast = rast, poly = poly)


###### Terrain index ####
# Explanation : we compute terrain indices from the SHOM MNT at 100m resolution.
# Methodological choices : 
# Lecours et al. 2017; Rees et al. (2014)
# We compute the following terrain indices: 
# slope, aspect, TPI, TRI and roughness.
# Slope and aspect are computed with 8 neighbors.
# TPI, TRI and roughness are automatically fixed to 8 neighbors.


# Compute terrain indices ----
# Data : we use the SHOM MNT at 100m resolution (loaded in previous section)

# Reproject the raster to WGS84
rast <- terra::project(rast, "EPSG:4326") # default method is "bilinear" if the first layer of rast is numeric (not categorical).

# Compute terrain indices 
compute_selected_terrain_ind(rast = rast, 
                             folder_path = "./data/processed_data/predictors/terrain_indices/SHOM_100m/", 
                             neighbors = 8, 
                             name = "SHOM100m_merged", # string to add in the filname before ".tif"
                             indices = c("slope", "aspect", "TRI", "TPI", "roughness"))

# [OPTIONAL] Check the results ----
# tif_files <- list.files("./data/processed_data/predictors/terrain_indices/SHOM_100m/", pattern = ".tif", full.names = TRUE)
# 
# rasters <- lapply(tif_files, rast)
# 
# names(rasters) <- tools::file_path_sans_ext(basename(tif_files))
# 
# # Plot each raster
# for (i in seq_along(rasters)) {
#   plot(rasters[[i]], main = names(rasters)[i])
# }
# 
# rm(tif_files, rasters)




# Extract terrain indices ----
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
    poly = poly,
    stats = c("min", "max", "mean")
  )
  
  # Extract only the columns: replicates + {layer}_min, _max, _mean
  tmp <- tmp %>%
    st_drop_geometry() %>%
    dplyr::select(replicates,
           all_of(paste0(layer, "_min")),
           all_of(paste0(layer, "_max")),
           all_of(paste0(layer, "_mean")))
  
  # Join to original buff_all
  buff_all <- left_join(buff_all, tmp, by = "replicates")
}
beepr::beep() # Beep to indicate completion


buff <- buff_all

rm(filelist_temp, tmp, layer, buff_all, poly, rast, i, var)
colnames(buff)


######################## Habitat #########################
# Explanation : we retrieve 3 different variables : 
# 1. main habitat within buffer, 
# 2. number of habitats per km² = number of different habitats / buffer area, 
# 3. surface proportion of each habitat within buffer --> To do so we need to apply transformation to avoid compositional data bias (see here for details : https://docs.google.com/document/d/1cN9vJ6I4fHzPXZfXjOm77Klk5hCSFBBkOj6Mhgum_S8/edit?tab=t.0 and https://medium.com/@nextgendatascientist/a-guide-for-data-scientists-log-ratio-transformations-in-machine-learning-a2db44e2a455) 

# ATTENTION : Important methodological considerations :
# "la carto de la cymodocée est tres incomplete. je déconseille de l'utiliser comme predicteur" Julie Deter --> we remove it from the analysis
# "Zone bathyale" is not considered as a habitat in this analysis. 
# Both are ignored when they appear to be the main habitat (we take the next more important habitat) and when counting the number of habitats per km².

# We compute these variables twice : once for detailed habitats and once for the main habitats categories (ie grouped habitats).

# Data :
# A 100m resolution raster file from Andromed with 14 bands representing the surface of different habitats in m² :
# 1	Association a rhodolithes
# 2	Association de la matte morte de Posidonia oceanica
# 3	Biocenose Coralligene
# 4	Biocenose de l herbier a Posidonia oceanica
# 5	Biocenose de la roche du large
# 6	Biocenose des algues infralittorales
# 7	Biocenose des galets infralittoraux
# 8	Fonds meubles circalittoraux
# 9	Fonds meubles infralittoraux
# 10	Habitats artificiels
# 11	Herbiers a Cymodocees
# 12	Zone bathyale
# 13	Herbier mixte a Zostera noltii, Zostera marina, Cymodocea nodosa et Ruppia cirrhosa
# 14	Herbiers a Zostera noltei

# Load raster ----
rast <- terra::rast("./data/raw_data/predictors/Habitat/bioc_2023_medfr_100m_2154.tif")

# Create 1 layer of presence/absence for each habitat
rast <- rast %>% as.factor() %>% segregate()

# Rename bands
names(rast) <-  c("Association rhodolithes", "matte_morte_P.oceanica",
                     "Coralligene","P.oceanica","Roche du large", "Algues infralittorales",
                     "Galets infralittoraux","fonds meubles circalittoraux",
                     "Fonds meubles infralittoraux", "Habitats artificiels",
                     "Herbiers Cymodocess", "Zone bathyale","Herbier mixte", "Z.noltei")


plot(rast)






# Group habitats ----
# Here habitats are separatly named as a function of their location (infralittoral, circalittoral, large) and type (rocky, soft bottom, seagrass, algae, coralligenous) --> we group habitats to keep only types. 

rast_grouped <- rast

# Put "fonds meubles" together
rast_grouped$soft_bottom <- app(rast_grouped[[names(rast_grouped) %in% c("fonds meubles circalittoraux","Fonds meubles infralittoraux")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("fonds meubles circalittoraux","Fonds meubles infralittoraux"))]]

# Put "P.oceanica" "Herbier mixte" "Z.noltei" together
rast_grouped$meadow <- app(rast_grouped[[names(rast_grouped) %in% c("P.oceanica", "Z.noltei", "Herbier mixte" )]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("P.oceanica", "Z.noltei", "Herbier mixte"))]]

# Put "Roche du large" and "Galets infralittoraux" together 
rast_grouped$rock <- app(rast_grouped[[names(rast_grouped) %in% c("Roche du large","Galets infralittoraux")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("Roche du large","Galets infralittoraux"))]]

# Put "Association rhodolithes" and "Coralligene" together
rast_grouped$coralligenous <- app(rast_grouped[[names(rast_grouped) %in% c("Association rhodolithes","Coralligene")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("Association rhodolithes","Coralligene"))]]












# Main habitat and Number of habitat / km² ----
# Important methodological considerations : "Zone bathyale" and "Cymodocees" are not considered as a habitat in this analysis, thus they are ignored when tehy appear to be the main habitat (we take the next more important habitat) and when counting the number of habitats per km².

calculate_habitats <- function(buff, rast, id_column_name, name = "", buff_area) {
  # Check inputs
  if (!inherits(buff, "sf") && !inherits(buff, "SpatVector")) stop("buff must be an sf or SpatVector object")
  if (!inherits(rast, "SpatRaster")) stop("rast must be a SpatRaster object")
  if (!id_column_name %in% colnames(buff)) stop(paste("Column", id_column_name, "not found in buff"))
  if (!buff_area %in% colnames(buff)) stop(paste("Area column", buff_area, "not found in buff"))
  if (!is.character(name)) stop("name must be a character string")
  
  # Define layers to exclude
  exclude_layers <- c("Zone bathyale", "Herbiers Cymodocess")
  
  # Filter out excluded layers
  included_bands <- setdiff(names(rast), exclude_layers)
  rast <- rast[[included_bands]]  # Subset raster to included layers only
  
  # Define dynamic column names
  main_habitat_col <- paste0(name, "main_habitat")
  nb_habitat_col <- paste0(name, "nb_habitat_per_km2")
  
  # Initialize output columns
  buff[[main_habitat_col]] <- NA_character_
  buff[[nb_habitat_col]] <- NA_real_
  
  # Total number of raster bands (after exclusion)
  total_bands <- terra::nlyr(rast)
  band_names <- names(rast)
  
  for (i in seq_len(nrow(buff))) {
    polygon_id <- buff[[id_column_name]][i]
    area_km2 <- buff[[buff_area]][i]
    
    message(
      "--------------------------------------------\n",
      "Processed polygon: ", polygon_id, "\n",
      "Polygon number: ", i, "\n",
      "--------------------------------------------"
    )
    
    polygon_sums <- list()
    
    for (band in seq_len(total_bands)) {
      band_name <- band_names[band]
      raster_band <- terra::subset(rast, band)
      
      val <- terra::extract(raster_band, buff[i, ], weights = TRUE, na.rm = TRUE)
      
      if (!is.null(val) && nrow(val) > 0) {
        weighted_values <- val[, 2] * val[, 3]
        weighted_sum <- if (all(is.na(weighted_values))) 0 else sum(weighted_values, na.rm = TRUE)
      } else {
        weighted_sum <- 0
      }
      
      polygon_sums[[band_name]] <- weighted_sum
      message(band_name, " weighted sum: ", weighted_sum)
    }
    
    sums_vector <- unlist(polygon_sums)
    max_band <- if (all(sums_vector == 0)) NA else names(sums_vector)[which.max(sums_vector)]
    number_habitat <- sum(sums_vector != 0)
    
    habitat_density <- if (is.na(area_km2) || area_km2 == 0) NA else number_habitat / area_km2
    
    buff[[main_habitat_col]][i] <- max_band
    buff[[nb_habitat_col]][i] <- habitat_density
    
    message(" - Main habitat for ", polygon_id, ": ", ifelse(is.na(max_band), "None", max_band))
    message(" - Habitat density (/km2) for ", polygon_id, ": ", habitat_density)
  }
  
  return(buff)
}

# For all habitats
buff1 <- calculate_habitats(buff = buff, 
                            rast = rast, 
                            id_column_name = "replicates",
                           buff_area = "area_km2")


# For grouped habitats
buff1 <- calculate_habitats(buff = buff1, 
                            rast = rast_grouped, 
                            id_column_name = "replicates", 
                           name = "grouped_",
                           buff_area = "area_km2")


































# Surface proportion of each habitat  (ALR transformed) #### TO DO #### ----


######################## BROUILLON / Habitat : surface of each habitat ####
# Load data
rast <- terra::rast("./data/raw_data/predictors/Habitat_and_Anchoring/bioc_landscape_indices_bathy_anchorings_2023_medfr_1000m_2154.tif", lyrs = c(1:7))
rast <- terra::project(rast, "EPSG:4326")

# Rename bands
names(rast) <- c("Posidonia", "Coralligeneous", "Rocks", "Sand", "Dead_Matte", "Other_Seagrass", "Infralitoral_Algae")
plot(rast)

# df <- readRDS("./data/processed_data/data_prep/Med_TOT_2023_P0_R0_Idalongeville.rds")
# poly <- readRDS('./data/processed_data/data_prep/Med_TOT_2023_P0.rds')
# poly <- poly %>% dplyr::select(spygen_code, geometry)

# Convert poly df to SpatVector
# poly <- st_as_sf(poly)
# poly <- terra::vect(poly)
# dim(poly)
# plot(poly)
# class(poly)

# Iterate through each habitat layer in the raster
for (layer in names(rast)) {
  # Call spatial_extraction function for the current habitat layer
  df <- spatial_extraction(
    var = layer,  # Current habitat name
    rast = rast[[layer]],  
    df = df,  
    poly = poly, 
    stats = c("mean")  # Averages the pixel values within the buffer
  )
}

# Check the updated data frame
print(colnames(df))

























######################## GFW : Vessel Presence (0.01°) : weighted mean, min, max, range ####
# Parameters
var <- "Vessel_Presence"
rast <- terra::rast("./data/raw_data/predictors/Vessels/public-global-presence-v20231026_202301-202401_Vessel_Presence_0.01_degres.tif.tif")
rast <- terra::project(rast, "EPSG:4326")
df <- df
poly <- buff

# Extraction
df <- spatial_extraction(var = var, rast = rast, df = df, poly = poly)

# [For Updates---]
write.csv(df, "./data/processed_data/data_prep/predictors/Extracted_Predictors/mtdt_101220242_COV_spatial_cov_no_hab_no_fishing_eff.csv")


######################## GFW : Apparent Fishing Effort (0.01°) : weighted mean ####
df_1 <- df
#  /!\ WARNING : RUN INDEPENDENTLY 
# This means that you need to save your previous df results and start from a new df with no other covariables saved in it. 
# This is because, when binding the df_NA and df, column index number is used. Thus, if you have a different number of column it won't work.

### 1. Extraction of weighted mean
# Parameters 
rast <- terra::rast("./data/raw_data/predictors/Vessels/public-global-fishing-effort-v20231026_Apparent_fishing_effort_202301-202401_0.01_degres.tif")
rast <- terra::project(rast, "EPSG:4326")
poly <- buff
df <- df
var <- "Fishing_Eff"

# Extraction
df <- spatial_extraction(var = var, rast = rast, df = df, poly = poly)



### 2. Extraction of nearest value for NA 
# Parameters
rast <- raster::raster("./data/raw_data/predictors/Vessels/public-global-fishing-effort-v20231026_Apparent_fishing_effort_202301-202401_0.01_degres.tif")
rast <- raster::projectRaster(rast, crs = "EPSG:4326")
df <- df
df_NA <- df[is.na(df$Fishing_Eff_mean), ] # Keep df rows in which Fishing_Eff is NA
poly_NA <- poly |> dplyr::filter(spygen_code %in% df_NA$spygen_code) # Keep only poly that have their spygen_code in df_NA

# Extraction
df_NA <- nearest_value(var = var, rast = rast, df = df_NA, poly = poly_NA)



### 3. Binding both results
df_NA <- df_NA |> dplyr::select(spygen_code, Fishing_Eff_value)
df_binded <- merge(df, df_NA, by = "spygen_code", all = TRUE)

df_binded$Fishing_Eff <- ifelse(is.na(df_binded$Fishing_Eff_mean), 
                                as.numeric(df_binded$Fishing_Eff_value), 
                                as.numeric(df_binded$Fishing_Eff_mean))

df <- df_binded |> dplyr::select(spygen_code, date, Fishing_Eff)

df_2 <- df
write.csv(df, "./data/processed_data/data_prep/predictors/Extracted_Predictors/mtdt_101220242_COV_fishing_eff.csv")


################ Check data for homogenous sampling effort ##############

buff <- as.data.frame(buff)
buff %>%
  group_by(PCR_replicates) %>%
  summarise(count = n()) %>%
  print()


buff %>%
  group_by(estimated_volume_total) %>%
  summarise(count = n()) %>%
  print()

buff %>%
  filter(estimated_volume_total > 55 & estimated_volume_total < 67 & PCR_replicates == 24) %>%
  dim()


summary(buff$area_km2)
hist(buff$area_km2, breaks = 10)
sd(buff$area_km2, na.rm = TRUE)


buff$area_km2 <- as.numeric(buff$area_km2)
buff %>%
  filter(estimated_volume_total > 55 & estimated_volume_total < 67 & PCR_replicates == 24 & area_km2 < 2) %>%
  dim()

