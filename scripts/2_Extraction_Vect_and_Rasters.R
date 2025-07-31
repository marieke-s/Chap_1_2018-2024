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



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- VECTOR DATA ----------------

###### Distance to shore ############
# Explanation : distance to shore is computed from the coastline shapefile using '2.1_Distance_shore.py'.
# Important methodological considerations : 
# Distance to shore is computed from replicates group's buffer centroids.
# When a centroid is on land (which happened for ~ 67/792 centroids) we set the distance to shore to 1 meter because it was actually not sampled on land (obviously) and it ended up there only because we simplified transects into straight lines between transect start and end points.

# Run 2.1_Distance_shore.py ----
# Basic bash call
system("python3 ./scripts/2.1_Distance_shore.py")


# Load buff with dist_shore  ----
buff <- st_read("./data/processed_data/predictors/mtdt_5_dist-shore.gpkg")
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

######################## Vessel presence (Luka) ############
############## MPA : in or out reserve : LAURE #############
# Explanation : variable that states whether or not the replicates were done within a fully protected MPA (ie a reserve).

# MPA data used : protectionMed_fr_modif.gpkg made by Lola Romant during her internship. She assigned protection levels to the French Mediterranean MPA based on MPA Guide Categories. This data was then modified by Laure Velez to keep only the MPA with the highest protection level when several were overlapping (protectionMed_fr.shp). Then it was manually (QGIS) modified by Marieke Schultz to keep correct like this : 
# Passe de la Réserve naturelle de Cerbicale (Bouches De Bonifacio) du niveau 3 au niveau 4
# Passe de l'Archipel des îles de Lavezzi (Bouches De Bonifacio) du niveau 3 au niveau 4
# J'ai créée la colonne "Fully" dans laquelle les niveaux 1, 2 et 3 sont en "YES" et le reste en "NO"
# --> protectionMed_fr_modif.gpkg

# In/Out data used : mtd_2018_2024_ampfully.csv made by Laure Velez in 07/2025. She checked for each eDNA samples whether or not they belonged to a fully protected MPA, cheking both start point, end point, and transect (eg when start and end points were outside but transect crossing the reserve). 

# Load In/Out data
in_out <- read.csv("./data/raw_data/predictors/MPA/mtd_2018_2024_ampfully.csv")
unique(in_out$mpa_fully)


# Here we assign the variable "mpa_fully" to the buff data. 
# When all samples of the replicates have mpa_fully = 1 we set buff$mpa_fully = 1, when all samples of the replicates have mpa_fully = 0 we set buff$mpa_fully = 0, else we set buff$mpa_fully = NA.

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






##############  /!\ MISSING IPOCOM IN LAURE'S DATA ######################## 

######################## Canyon : distance to canyon : MARTIN ####


######################## Port : distance to port : MARTIN ####

######################## MPA : distance to reserve : MARTIN ####

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


######################## Terrain index : slope, aspect, TPI, TRI, roughness ####
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

rm(filelist_temp, tmp, layer, buff_all)



######################## Habitat : LOIC #########################
# Explanation : we retreive different variables : 
# main habitat within buffer, number of habitats within buffer, surface proportion of each habitat within buffer. 
# We do this twice : once for detailed habitats and once for the main habitats categories (ie grouped habitats)
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

# Group habitats ----
rast_grouped <- rast

# Put "fonds meubles" together
rast_grouped$soft_bottom <- app(rast_grouped[[names(rast_grouped) %in% c("fonds meubles circalittoraux","Fonds meubles infralittoraux")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("fonds meubles circalittoraux","Fonds meubles infralittoraux"))]]

# Put "P.oceanica" "Herbier mixte" "Z.noltei" and "Herbiers Cymodocess" together
rast_grouped$meadow <- app(rast_grouped[[names(rast_grouped) %in% c("P.oceanica","Herbiers Cymodocess", "Z.noltei", "Herbier mixte" )]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("P.oceanica","Herbiers Cymodocess", "Z.noltei", "Herbier mixte"))]]

# Put "Roche du large" and "Galets infralittoraux" together 
rast_grouped$rock <- app(rast_grouped[[names(rast_grouped) %in% c("Roche du large","Galets infralittoraux")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("Roche du large","Galets infralittoraux"))]]

# Put "Association rhodolithes" and "Coralligene" together
rast_grouped$coralligenous <- app(rast_grouped[[names(rast_grouped) %in% c("Association rhodolithes","Coralligene")]], fun = "sum")
rast_grouped <- rast_grouped[[!(names(rast_grouped) %in% c("Association rhodolithes","Coralligene"))]]

plot(rast_grouped)


## main biocenose in a cell ----

fun <- function(i) { i / sum(i) } # function to compute proportion of each habitats in the new cells (with lower resolution)

r_bioce <- app(r_bioce, fun = fun)

r_max_bioce <- which.max(r_bioce)
names(r_max_bioce) <- "max_bioce"

values(r_max_bioce) <- case_when(values(r_max_bioce)== 1 ~ "Association matte morte P.oceanica",
                                 values(r_max_bioce)== 3 ~ "Algues infralittorales",
                                 values(r_max_bioce) == 6 ~ "zone bathyale",
                                 values(r_max_bioce) == 7 ~ "meadow",
                                 values(r_max_bioce) == 8 ~ "soft bottom",
                                 values(r_max_bioce) == 9 ~ "corralligenous")

# save into a raster file (class SpatRaster)
writeRaster(r_max_bioce, file = "Outputs/01_read_cleaning/r_max_bioce.tif", overwrite = TRUE)

## number of biocenoses per cell ----

fun1 <- function(i) {sum(i !=0)}

r_n_bioce <- app(r_bioce, fun1)

names(r_n_bioce) <- "n_bioce"

# save into a raster file (class SpatRaster)
writeRaster(r_n_bioce, file = "Outputs/01_read_cleaning/r_n_bioce.tif", overwrite = TRUE)



######################## Habitat : main habitat and number of habitat ####
# Data :
# A raster file with 7 bands representing the surface of different habitats in m² :
# 1 posidonia_seagrass	surface of posidonia in m²	
# 2 coralligeneous	surface of coralligeneous in m²	
# 3 rocks	surface of rocks in m²
# 4 sand	surface of sand in m²	
# 5 dead_matte	surface of dead_matte in m²	
# 5 other_seagrass	surface of other_seagrass in m²	
# 6 infralitoral_algae	surface of infralitoral_algae in m²	
# 7 infralitoral_algae	surface of infralitoral_algae in m²
# from Andromede

# Load necessary data
rast <- terra::rast("./data/raw_data/predictors/Habitat_and_Anchoring/bioc_landscape_indices_bathy_anchorings_2023_medfr_1000m_2154.tif", lyrs = c(1:7))
rast <- terra::project(rast, "EPSG:4326")
df <- df
poly <- as(buff, "SpatVector")

# Rename bands
names(rast) <- c("Posidonia", "Coralligeneous", "Rocks", "Sand", "Dead_Matte", "Other_Seagrass", "Infralitoral_Algae")

# Initialize results columns in df
df$main_habitat <- NA
df$number_habitat <- NA

# Find main habitat and number of habitats for each polygon
for (i in 1:nrow(poly)) {  
  spygn_code <- poly$spygen_code[i]
  print(paste(
    "--------------------------------------------",
    "Processed polygon:", spygn_code,
    "Polygon number:", i,
    "--------------------------------------------"
  ))
  
  polygon_sums <- list()
  
  for (band in 1:7) {
    band_name <- names(rast)[band]
    raster_band <- terra::subset(rast, band)
    
    # Extraction with weights
    val <- terra::extract(raster_band, poly[i, ], weights = TRUE, na.rm = TRUE)
    
    # Check if extraction returned any values
    if (!is.null(val) && nrow(val) > 0) {
      # Compute weighted sum, ignoring NA values
      weighted_values <- val[, 2] * val[, 3]
      weighted_sum <- sum(weighted_values, na.rm = TRUE)
      
      # If all values were NA, set weighted_sum to 0
      if (all(is.na(weighted_values))) {
        weighted_sum <- 0
      }
    } else {
      weighted_sum <- 0
    }
    
    polygon_sums[[band_name]] <- weighted_sum
    
    print(paste(band_name, "weighted sum:", weighted_sum))
  }
  
  # Convert polygon_sums to a named numeric vector
  sums_vector <- unlist(polygon_sums)
  
  # Determine the main habitat (band with the highest weighted sum)
  # If all sums are zero, set main_habitat to NA or a specific value indicating no habitat
  if (all(sums_vector == 0)) {
    max_band <- NA
  } else {
    max_band <- names(sums_vector)[which.max(sums_vector)]
  }
  
  # Determine the number of habitats with non-zero sums
  number_habitat <- sum(sums_vector != 0)
  
  # Store results in the corresponding row in df
  df$main_habitat[df$spygen_code == spygn_code] <- max_band
  df$number_habitat[df$spygen_code == spygn_code] <- number_habitat
  
  print(paste(" - Main habitat for", spygn_code, ":", ifelse(is.na(max_band), "None", max_band)))
  print(paste(" - Number of habitats for", spygn_code, ":", number_habitat))
}


write.csv(df, "./data/processed_data/data_prep/predictors/Extracted_Predictors/mtdt_101220242_COV_protection_habitat.csv")



######################## Habitat : surface of each habitat ####
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

buff$area_km2 <- as.numeric(buff$area_km2)
buff %>%
  filter(estimated_volume_total > 55 & estimated_volume_total < 67 & PCR_replicates == 24 & area_km2 < 2) %>%
  dim()

