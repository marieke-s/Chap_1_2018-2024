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

# Contributors : Martin Paquet, Pauline Viguier, Laure Velez, Marie Orblin

# Date script created: 25/07/2025

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

#------------- Load and prep data ------------------
# Load mtdt_5 
buff <- st_read("./data/processed_data/Mtdt/mtdt_5.gpkg")

colnames(buff)

# Keep only replicates and geom 
buff <- buff %>%
  dplyr::select(replicates, geom)






#------------- VECTOR DATA ----------------
###### Sampling effort : Buffer area ############
# Explanation : the aim is to compute variables that represent the sampling effort in each replicate group : 
# 1. Total volume filtered in the replicate group = estimated_volume_total -->  Already exists
# 2. Nb of PCR replicates = 12 * nb of pooled samples + 12 nb of unpooled samples --> computed here 
# 3. Area covered = area of replicate group buffer --> computed here


# Project in CRS=2154
buff <- st_transform(buff, crs = 2154)

# Compute area in km2
buff$area_km2 <- st_area(buff) / 1e6  # Convert from m^2 to km^2
buff$area_km2 <- units::set_units(st_area(buff), km^2)
buff$area_km2 <- as.numeric(buff$area_km2)  # Convert to numeric for easier handling



###### Reserve : in or out (Laure) #############
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
  as.data.frame() %>%
  group_by(mpa_fully) %>%
  summarise(count = n()) %>%
  print()

# print the replicates column of buff where mpa_fully is NA
any(is.nan(buff$mpa_fully))
buff %>% filter(is.na(mpa_fully)) %>% dplyr::select(replicates, mpa_fully) %>%
  print()


###### Distances : to port, canyon and reserve (Martin) ############
# Explanation : distances to several entities (port, canyons, MPA reserve) were computed by Martin Paquet in 07/2025 with the distance pipeline explained in Methods.txt.
# Documentation port : https://services.sandre.eaufrance.fr/telechargement/geo/PTS/sandre_dictionnaire_PTS_2.pdf
# Documentation canyon : https://opdgig.dos.ny.gov/datasets/250f2b9496854bd098be8154bab04a3a_4/about
# Load data ----
dist <- st_read("./data/raw_data/predictors/Distances/buffer_with_closest_feats_search_outlier_treshold_shore_50m.gpkg")
dist <- as.data.frame(dist)

# Clean ----

dist <- dist %>%
  
  dplyr::select(-c ("ID","date","time_start","depth_sampling","depth_seafloor","lockdown","BiodivMed2023","method","country","region","site","subsite","component","habitat","protection","project", "Tele01","Pleo","Mamm01","Vert01",  "Tele01" ,"Pleo","Mamm01","Vert01","X16s_Metazoa","Bact02" ,"Euka02" ,"duration_total","comments","estimated_volume_total", "mpa_name...17"  )) %>%  # Remove mtdt cols
  
  dplyr::select(-"canyon_Geomorphic") %>% # because gives no info : unique(dist$canyon_Geomorphic) # "Canyon"
  dplyr::select(-"canyon_Ocean") %>% # because all "Mediterranean Sea"
  dplyr::select(-c("port_NumTexteRe", "port_DatePubliT", "port_CdTypeText", "port_MnTypeText", "port_URLTexteRe")) %>% # full NA columns
  dplyr::select(-c("port_CdZonePort", "port_MnTypeZone", "port_DtCreatZon", "port_DtMajZoneP", "port_ComZonePor", "port_CdPort", "port_StatutZone", "port_CdTypeZone", "port_gid", "port_aire_m2")) %>% # We delete metadata associated to ports because incomplete
  rename("mpa_name" = "mpa_name...81") %>%
  dplyr::select(-c("canyon_dist_min_from_buff_m", "port_dist_min_from_buff_m", "mpa_dist_min_from_buff_m")) %>% # These cols were used for methodological checks only ([MARTIN CHECK])calcul de la distance du accCost depuis le buffer à la géométrie de variable la plus proche) 
  dplyr::select(-c("canyon_OBJECTID", "canyon_Canyon_ID", "port_NomZonePor", "port_NomPort", "mpa_name")) # we don't need these metadata cols



  
  
  
  
  
  
# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(dist), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(dist %>% dplyr::select(replicates, all_of(extra_cols)), by = "replicates")

rm(dist, extra_cols)


# Count NA per column
sapply(buff, function(x) sum(is.na(x)))














####### Boats (Luka) ############
# Load csv file
boats <- read.csv("./data/processed_data/predictors/Boats/mtdt_7_boats_month_ALL.csv")

# Merge to buff 

buff <- buff %>% left_join(boats, by = "replicates")
  
names(buff)






#------------- RASTER DATA ----------------

###### Gravity (1 km) : weighted mean, min, max, range ####
# Data : From Laure Velez (MARBEC)


# Parameters
var <- "gravity"
rast <- terra::rast("./data/raw_data/predictors/Gravity/rastGravity.tif")
poly <- buff 

# Extraction
gr <- spatial_extraction(var = var, rast = rast, poly = poly, stat = c("min", "max", "mean", "range"))


# Merge -----
# Detect columns in gr that are not in buff
extra_cols <- setdiff(names(gr), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'

buff <- buff %>%
  left_join(
    gr %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )







rm(gr, extra_cols)



















###### Bathymetry (0.001°) : weighted mean, min, max, range ####
# Data : we use the SHOM MNT at 100m resolution

# Parameters
var <- "bathy"
rast <- terra::rast("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")
poly <- buff

# Set values >= to 0 to NA (because it's land)
rast[rast >= 0] <- NA

# Extraction
bat <- spatial_extraction(var = var, rast = rast, poly = poly)

# Clean -----
bat <- bat %>%
  dplyr::select(c("replicates", "bathy_min", "bathy_max", "bathy_mean", "bathy_range", ))

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
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )




rm(bat, extra_cols, rats, poly)

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

rm(filelist_temp, tmp, layer, poly, rast, i, var)







# Merge -----
# Detect columns in buff_all that are not in buff
extra_cols <- setdiff(names(buff_all), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    buff_all %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )


rm(buff_all, extra_cols)






















######  Habitat (100m) ####
# Explanation : we retrieve 3 different variables : 
# 1. main habitat within buffer, 
# 2. number of habitats per km² = number of different habitats / buffer area, 
# 3. surface proportion of each habitat within buffer 

# WARNING : Important methodological considerations :
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

# Remove "zone bathyale" and "Cymodocees" from raster
rast <- rast[[!(names(rast) %in% c("Zone bathyale", "Herbiers Cymodocess"))]]



# Export 
terra::writeRaster(rast_grouped, "./data/processed_data/predictors/Habitat/habitat_raster.tif")







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


plot(rast_grouped)



# Export 
terra::writeRaster(rast_grouped, "./data/processed_data/predictors/Habitat/grouped_habitat_raster.tif")










# Main habitat and Number of habitat / km² ----
# Important methodological considerations : "Zone bathyale" and "Cymodocees" are not considered as a habitat in this analysis, thus they are ignored when tehy appear to be the main habitat (we take the next more important habitat) and when counting the number of habitats per km².


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








# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(buff1), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    buff1 %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  )

rm(buff1, extra_cols)


























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


  
  
  


#--- making predictors_raw_v2.2 ----
# Load predictors_raw_v2.1
buff <- sf::st_read("./data/processed_data/predictors/predictors_raw_v2.1.gpkg")

# Handle NA ----
# For habitat with NA --> take the 3 closest pixel values

# Select buff with NA in habitat columns
habitat_cols <- c("main_habitat", "grouped_main_habitat", "grouped_nb_habitat_per_km2")
buff_na <- buff %>%
  filter(if_any(all_of(habitat_cols), is.na)) # 5 NA 

buff_na <- buff_na %>% 
  dplyr::select(replicates, area_km2)

# Export 
st_write(buff_na, "./data/processed_data/predictors/Habitat/habitat_NA.gpkg", delete_dsn = TRUE)

## QGIS check : far away from raster but considered coastal because near a small island.





# Compute main habitat and number of habitats per km2 for buff_na grouped and not grouped ----
source("./utils/Fct_Data-Prep.R")

buff1 <- calculate_habitats_nn(
  buff           = buff_na,
  rast           = rast_grouped,
  id_column_name = "replicates",
  buff_area      = "area_km2",
  name           = "grouped_"      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)

buff2 <- calculate_habitats_nn(
  buff           = buff_na,
  rast           = rast,
  id_column_name = "replicates",
  buff_area      = "area_km2",
  name           = ""      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)


# Merge back NA values -----

# Merge buff1 and buff 2 by replicates 
buff12 <- buff1 %>%
  dplyr::select(-area_km2) %>%
  left_join(
    buff2 %>% st_drop_geometry() %>% dplyr::select(-area_km2),
    by = "replicates"
  )


# For replicates in buff12, replace the hab_cols columns in buff with those in buff12

hab_cols <- colnames(buff12[, -1] %>% sf::st_drop_geometry())
hab_cols <- intersect(hab_cols, names(buff))


# match buff$replicates to buff12$replicates
idx <- match(buff$replicates, buff12$replicates)

# which rows in buff have a match in buff12?
rows_to_update <- !is.na(idx)

# replace values in buff[rows_to_update, hab_cols] 
# with the corresponding rows from buff12
buff[rows_to_update, hab_cols] <- buff12[idx[rows_to_update], hab_cols] %>% sf::st_drop_geometry()












# Check and Clean up -----
# Check na again and clean 

habitat_cols <- c("main_habitat", "grouped_main_habitat", "grouped_main_habitat", "grouped_nb_habitat_per_km2")
buff_na <- buff %>%
  filter(if_any(all_of(habitat_cols), is.na)) 

# if buff_na = 0 rows then clean up 

if(nrow(buff_na) == 0){
  message("No more NA in main habitat and nb hab per km2")
  rm(buff_na, buff1, buff2, buff12, hab_cols, habitat_cols, idx, rows_to_update)
} else {
  message( "/!\ Remaining NA !")
}


# Handle habitats proportions (grouped only) -----

# Habitats (base names, must match layers in rast_grouped)
habitats <- c(
  "algues_infralittorales",
  "matte_morte_p_oceanica",
  "coralligenous",
  "rock",
  "soft_bottom",
  "meadow",
  "habitats_artificiels"
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




# Clean up ----
# Check na and clean up


rm(h, habitats, i, idx, mean_cols, new_cols, no_overlap, vals, df, nn, rast_df, rast_pts, temp, col, col_mean, new_val, no_overlap_ids, total_means)
rm(rast, rast_grouped)


######  Habitat (100m) v2 ####
# Explanation : we retrieve 3 different variables : 
# 1. main habitat within buffer, 
# 2. number of habitats per km² = number of different habitats / buffer area, 
# 3. surface proportion of each habitat within buffer 

# WARNING : Important methodological considerations :
# "la carto de la cymodocée est tres incomplete. je déconseille de l'utiliser comme predicteur" Julie Deter --> we remove it from the analysis
# "Zone bathyale" is kept. 

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

# Load & prep predictors_raw_v2.2 ----
buff <- st_read("./data/processed_data/predictors/predictors_raw_v2.2.gpkg")

# Remove habitat cols from buff 
habitat_cols <- c("main_habitat", "nb_habitat_per_km2", 
                  "grouped_main_habitat", "grouped_nb_habitat_per_km2", "algues_infralittorales_mean","habitats_artificiels_mean","soft_bottom_mean","meadow_mean","rock_mean","coralligenous_mean", "matte_morte_p_oceanica_mean")

names(buff)

buff <- buff %>%
  dplyr::select(-any_of(habitat_cols))

# Load and prep rasters ----
rast <- terra::rast("./data/raw_data/predictors/Habitat/bioc_2023_medfr_100m_2154.tif")
plot(rast)

# Create 1 layer of presence/absence for each habitat
rast <- rast %>% as.factor() %>% segregate()

# Rename bands
names(rast) <-  c("Association rhodolithes", "matte_morte_P.oceanica",
                  "Coralligene","P.oceanica","Roche du large", "Algues infralittorales",
                  "Galets infralittoraux","fonds meubles circalittoraux",
                  "Fonds meubles infralittoraux", "Habitats artificiels",
                  "Herbiers Cymodocess", "Zone bathyale","Herbier mixte", "Z.noltei")

# Remove "Cymodocees" from raster
rast <- rast[[!(names(rast) %in% c("Herbiers Cymodocess"))]]



# Export 
terra::writeRaster(rast, "./data/processed_data/predictors/Habitat/habitat_raster_2.tif", overwrite=TRUE)







#   --- Group habitats ---
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

# Keep "Zone bathyale" alone 


plot(rast_grouped)



# Export 
terra::writeRaster(rast_grouped, "./data/processed_data/predictors/Habitat/grouped_habitat_raster_v2.tif", overwrite=TRUE)



# Load rasters directly ----
rast <- terra::rast("./data/processed_data/predictors/Habitat/habitat_raster_2.tif")
rast_grouped <- terra::rast("./data/processed_data/predictors/Habitat/grouped_habitat_raster_v2.tif")







# Main habitat & Nb of habitats/km2 ----
source("./utils/Fct_Data-Prep.R")

buff1 <- calculate_habitats_nn(
  buff           = buff,
  rast           = rast_grouped,
  id_column_name = "replicates",
  buff_area      = "area_km2",
  name           = "grouped_"      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)
st_write(buff1, "./data/processed_data/predictors/predictors_raw_v2.0_temp-buff1.gpkg", delete_dsn = TRUE)
rm(grouped_rast)
buff1 <- st_read("./data/processed_data/predictors/predictors_raw_v2.0_temp-buff1.gpkg")

buff2 <- calculate_habitats_nn(
  buff           = buff,
  rast           = rast,
  id_column_name = "replicates",
  buff_area      = "area_km2",
  name           = ""      # or e.g. "hab_"
  # k_nearest    = 3       # default is 3, can change if needed
)
st_write(buff2, "./data/processed_data/predictors/predictors_raw_v2.0_temp1.gpkg", delete_dsn = TRUE)
rm(rast)
buff2 <- st_read("./data/processed_data/predictors/predictors_raw_v2.0_temp1.gpkg")

names(buff2)

# Merge ----
extra_cols <- setdiff(names(buff1), names(buff))
extra_cols2 <- setdiff(names(buff2), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(
    buff1 %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols)),
    by = "replicates"
  ) %>%
  left_join(
    buff2 %>% 
      st_drop_geometry() %>% 
      dplyr::select(replicates, dplyr::all_of(extra_cols2)),
    by = "replicates"
  )

names(buff)
rm(buff1, buff2, extra_cols, extra_cols2)


# Habitats proportions (grouped only) -----

# Habitats (base names, must match layers in rast_grouped)
## 1. Compute weighted means of each habitat within buffers

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

#   | Case                              | Habitat value source                       | Final proportions      |
#   | --------------------------------- | ------------------------------------------ | ---------------------- |
#   | **Polygon overlaps raster**       | Weighted means from `spatial_extraction()` | mean_i / sum(mean_all) |
#   | **Polygon outside raster extent** | Sum of 3 nearest pixels per habitat        | sum_i / sum(sum_all)   |





# Clean up ----
# Check na and clean up


rm(h, habitats, i, idx, mean_cols, new_cols, no_overlap, vals, df, nn, rast_df, rast_pts, temp, col, col_mean, new_val, no_overlap_ids, total_means)
rm(rast, rast_grouped)





















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




#------------- NCDF DATA ----------------
###### MARS3D (1.2km) #####################################
# MARS3D data was extracted by Pauline using scripts "...". 

# Load MARS3D and combine data

# Current and Wind ----

# List all .geojson files 
geojson_files <- list.files("./data/raw_data/predictors/MARS3D/adne_extract_current_wind/",
                            pattern = "\\.geojson$", full.names = TRUE)


# Open all GeoJSONs and store in a list 
sf_list <- lapply(geojson_files, function(f) {
  sf::st_read(f, quiet = TRUE)
})

# Check name and dimensions of each geojson
for (i in seq(1:3)){
  print(sf_list[[i]][1,]$layer)
  print(dim(sf_list[[i]]))
}

# Combine 
cur_wind <- dplyr::bind_rows(sf_list)






# Temperature and Salinity ----
# List all .geojson files 
geojson_files <- list.files("./data/raw_data/predictors/MARS3D/adne_extract_sal_temp", pattern = "\\.geojson$", 
                            full.names = TRUE)

# Open all GeoJSONs and store in a list 
sf_list <- lapply(geojson_files[1:3], function(f) {
  sf::st_read(f, quiet = TRUE)
})

# Check name and dimensions of each geojson
for (i in seq(1:3)){
  print(sf_list[[i]][1,]$layer)
  print(dim(sf_list[[i]]))
}

# Combine 
sal_temp <- dplyr::bind_rows(sf_list)





# Combine sal_temp and cur_wind -----

# Detect columns in sal_temp that are not in cur_wind
extra_cols <- setdiff(names(sal_temp), names(cur_wind))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
mars3d <- cur_wind %>%
  as.data.frame() %>%
  left_join(as.data.frame(sal_temp) %>% dplyr::select(replicates, all_of(extra_cols)), by = "replicates")

# Keep only necessary columns
mars3d <- mars3d %>%
  dplyr::select(-c("date", "depth_sampling", "depth_seafloor", "datetime", "day_24h", "date_7j", "date_1mois", "date_1an", "path", "geometry"))





# Add mards3d to buff ----
buff <- left_join(buff, mars3d, by = "replicates")
colnames(buff)

# Clean
rm(cur_wind, sal_temp, mars3d, geojson_files, sf_list, extra_cols, i)
buff <- buff %>% dplyr::select(-geometry)




###### Chlorophyll (1km) ----
# Chlorophylle data is extracted from Copernicus data : cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D (DOI : https://doi.org/10.48670/moi-00300). It has a daily temporal resolution and a 1km spatial resolution. It is a L4 product meaning :
# - Level 4 data result from analyses of L3 data (e.g., variables derived from multiple measurements).
# - L4 are those products for which a temporal averaging method or an interpolation procedure is applied to fill in missing data values. Temporal averaging is performed on a monthly basis. The L4 daily products is also called “interpolated” or “cloud free” products. 
# Product user Manual (https://documentation.marine.copernicus.eu/PUM/CMEMS-OC-PUM.pdf)

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











###### SST (1km) ----
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



























#------------- pMEM [not saved in predictors] ------------- 
----------------
# Explanation : check Spatially explicit predictions using spatial eigenvector maps ( https://doi.org/10.1111/2041-210X.14413)

# Generate pMEM variables ------
## 1. Génération des descripteurs spatiaux (pMEM)
# Générer une métrique de distance (euclidienne ici)
dist_metric <- pMEM::genDistMetric()  # ou ajouter delta/angle pour asymétrie

# Choisir une fonction de pondération spatiale, ex : "Gaussian"
dwf <- pMEM::genDWF(fun = "Gaussian", range = 3)  # adapter le range selon l’échelle spatiale

# Prep coords from buff centroids
coords <- st_centroid(buff) %>%
  st_coordinates()

coords <- as.matrix(coords[, c("x", "y")])


## 2. Générer les eigenvecteurs spatiaux prédictifs (SEMap object)
sef <- genSEF(x = coords, m = dist_metric, f = dwf)


## 3. Extraire les vecteurs propres spatiaux (pMEMs) 

# Convertir en data.frame 
pMEM_vars <- as.data.frame(sef) # 681 variables added /!\

# Merge ----
# Detect columns in pMEM_vars that are not in buff
extra_cols <- setdiff(names(pMEM_vars), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff1 <- buff %>%
  cbind(
    pMEM_vars %>% 
      dplyr::select(dplyr::all_of(extra_cols)) 
  )




rm(pMEM_vars, sef, coords)

















#------------- EXPORT PREDICTORS DATA ----------------

# CLEAN UP BUFF -----

x <- buff  

# Clean and uniquify column names in a way that matches GDAL/SQLite expectations
nm <- names(x)
nm <- enc2utf8(nm)
nm <- gsub("[^A-Za-z0-9_]", "_", nm)  # replace punctuation (e.g., dots) with underscore
nm <- tolower(nm)                     # SQLite is case-insensitive; normalize to lower
nm <- make.unique(nm, sep = "_")      # ensure uniqueness after normalization
names(x) <- nm

x <- x %>% mutate(across(c(mpa_fully), as.character))

# Make sure there are no list-cols
is_listcol <- vapply(x, inherits, logical(1), "list")
if (any(is_listcol)) {
  stop("Some columns are list-columns: ", paste(names(x)[is_listcol], collapse=", "),
       ". Unnest or convert before writing.")
}









# Export predictors_raw_v1_0 ----

# predictors_raw_v1_0
# 06/11/2025. 
# Based on mtdt_5.gpkg (788 replicates buffers).
# Predictors (185) :   area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_objectid, canyon_area_km2, canyon_delta_d, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_canyon_id, canyon_shape_length, canyon_shape_area, port_nomzonepor, port_nomport, mpa_protection, mpa_name, mpa_fully_1, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, bathy_mean, bathy_range, bathy_min, bathy_max, date_x, time_start, depth_sampling_x, depth_seafloor_x, lockdown, biodivmed2023, method, country, region, site, subsite, component, habitat, protection, mpa_name_x, project, tele01, pleo, mamm01, vert01, x16s_metazoa, bact02, euka02, estimated_volume_total, duration_total, comments, mpa_name_y, dist_seabed_depthsampling, date_y, depth_sampling_y, depth_seafloor_y, datetime, day_24h, date_7j, date_1mois, date_1an, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y
cat(colnames(x), sep = ", ")


# Ensure file is not locked; remove it before writing
gpkg_path <- "./data/processed_data/predictors/predictors_raw_v1_0.gpkg"
layer_name <- "predictors_raw_v1_0"
if (file.exists(gpkg_path)) unlink(gpkg_path)

st_write(x, dsn = gpkg_path, layer = layer_name, driver = "GPKG", append = FALSE)

























# Export predictors_raw_v1.1 ----
# predictors_raw_v1.1
# 08/11/2025. 
# Based on mtdt_5.gpkg (788 replicates buffers).
# x, y converted to decimal degree
# Predictors (185) : area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_objectid, canyon_area_km2, canyon_delta_d, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_canyon_id, canyon_shape_length, canyon_shape_area, port_nomzonepor, port_nomport, mpa_protection, mpa_name, mpa_fully_1, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, bathy_mean, bathy_range, bathy_min, bathy_max, date_x, time_start, depth_sampling_x, depth_seafloor_x, lockdown, biodivmed2023, method, country, region, site, subsite, component, habitat, protection, mpa_name_x, project, tele01, pleo, mamm01, vert01, x16s_metazoa, bact02, euka02, estimated_volume_total, duration_total, comments, mpa_name_y, dist_seabed_depthsampling, date_y, depth_sampling_y, depth_seafloor_y, datetime, day_24h, date_7j, date_1mois, date_1an, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y

st_write(buff, "./data/processed_data/predictors/predictors_raw_v1.1.gpkg", append = FALSE)



# Export predictors_raw_v1.2 ----
# predictors_raw_v1.2
# 11/11/2025. 
# Based on mtdt_5.gpkg (788 replicates buffers).
# same as redictors_raw_v1.1 but removing columns linked to metadata + some other useless cols

# Predictors (152) : replicates, area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_area_km2, canyon_delta_d, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_canyon_id, canyon_shape_length, canyon_shape_area, port_nomzonepor, port_nomport, mpa_protection, mpa_name, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, bathy_mean, bathy_range, bathy_min, bathy_max, dist_seabed_depthsampling, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y, geom
cat(colnames(pred), sep = ", ")


pred <- st_read("./data/processed_data/predictors/predictors_raw_v1.1.gpkg") %>%
  dplyr::select(-c(date_x, time_start, depth_sampling_x, depth_seafloor_x, lockdown, biodivmed2023, method, country, region, site, subsite, component, habitat, protection, mpa_name_x, project, tele01, pleo, mamm01, vert01, x16s_metazoa, bact02, euka02, estimated_volume_total, duration_total, comments, mpa_name_y, date_y, depth_sampling_y, depth_seafloor_y, datetime)) %>% # mtdt columns
  dplyr::select(-c(day_24h, date_7j, date_1mois, date_1an))%>% # mars3d useless cols
  dplyr::select(-mpa_fully_1)%>% # duplicate of mpa_fully
  dplyr::select(-c(canyon_objectid)) # useless canyon column



st_write(pred, "./data/processed_data/predictors/predictors_raw_v1.2.gpkg", append = FALSE)


# Export predictors_raw_v2.0 ----
# 11/11/2025. 
# ROWS
# Based on mtdt_7.gpkg 
# MODIFICATION : rows from Ange2Mer and no detection are removed
# PREDICTORS : 
# same as predictors_raw_v1.2 regarding the predictors
# replicates : predictor's replicates changed to incorporate removed samples (ie to correspond to mtdt$replicates instead of mtdt$replicates_old)

# Predictors (152) : replicates, area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_area_km2, canyon_delta_d, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_canyon_id, canyon_shape_length, canyon_shape_area, port_nomzonepor, port_nomport, mpa_protection, mpa_name, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, bathy_mean, bathy_range, bathy_min, bathy_max, dist_seabed_depthsampling, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y, geom
cat(colnames(pred), sep = ", ")


## 1. Same columns as predictors_raw_v1.2 : 
pred <- st_read("./data/processed_data/predictors/predictors_raw_v1.2.gpkg")

## 2. Removing rows not in mtdt_7$replicates_old (deleted because either no detection nor ange2mer) :
mtdt_7 <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")
pred <- pred %>% filter(replicates %in% mtdt_7$replicates_old)

setdiff(pred$replicates, mtdt_7$replicates_old) # 0 : now we have the same nb of rows for predictors and mtdt_7 (and thus occ_pooled_v1.1 and div_indices_v1.0)


## 3. Clean predictor's replicate column
## 3.1. Add mtdt_7$replicates correct replicates column to predictors_raw (needed for later)
mtdt_7$replicates_c <- mtdt_7$replicates

pred <- pred %>%
  left_join(
    mtdt_7 %>%
      st_drop_geometry() %>%
      dplyr::select(replicates_c, replicates_old),
    by = c("replicates" = "replicates_old")
  )

## 3.2. Reorder predictor's replicates (as in mtdt_7$replicates) :
pred$replicates <- sapply(pred$replicates, reorder_replicates)

length(setdiff(pred$replicates, mtdt_7$replicates_old)) # 4 : now we differ a bit from replicates_old because we reodered the samples
length(setdiff(pred$replicates, mtdt_7$replicates)) # 17 : but we still differ from replicates because we didnt take into account missing samples. 
length(setdiff(mtdt_7$replicates, pred$replicates)) # 17


## 3.3  Replace predictor's replicates by 'new' replicates (mtdt_7$replicates_c) to account for removed no-detection-samples  within replicate groups and for SPY201174/SPY201181/SPY202491SPY202495 correction : 

pred <- pred %>% 
  dplyr::select(-replicates) %>%
  rename("replicates" = "replicates_c")
  
length(setdiff(pred$replicates, mtdt_7$replicates)) # 0
length(setdiff(mtdt_7$replicates, pred$replicates)) # 0

# Check with occ_pooled_v1.1 and div_indices_v1.0
occ_pooled_v1.1 <- st_read("./data/processed_data/eDNA/occ_pooled_v1.1.csv")
div_indices_v1.0  <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv")

length(setdiff(pred$replicates, occ_pooled_v1.1$replicates)) # 0
length(setdiff(pred$replicates, div_indices_v1.0$replicates)) # 0

# All good !

st_write(pred, "./data/processed_data/predictors/predictors_raw_v2.0.gpkg", append = FALSE)


# Export predictors_raw_v2.1 ----
# 12/11/2025. 
# ROWS
# Based on mtdt_7.gpkg, same as predictors_raw_v2.0
# PREDICTORS : 
# same as predictors_raw_v1.2 - name cols : "canyon_canyon_id"     "port_nomzonepor"      "port_nomport"  "mpa_name" and - "canyon_delta_d"
# predictors (147) : area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_area_km2, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_shape_length, canyon_shape_area, mpa_protection, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, bathy_mean, bathy_range, bathy_min, bathy_max, dist_seabed_depthsampling, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y, replicates, geom
cat(colnames(pred), sep = ", ")

pred <- st_read("./data/processed_data/predictors/predictors_raw_v2.0.gpkg") %>%
  dplyr::select(-c(canyon_canyon_id, port_nomzonepor, port_nomport, mpa_name)) %>%
  dplyr::select(-canyon_delta_d)

# export
st_write(pred, "./data/processed_data/predictors/predictors_raw_v2.1.gpkg", append = FALSE)

# Export predictors_raw_v2.2 ----
# 12/11/2025. 
# ROWS
# Based on mtdt_7.gpkg, same as predictors_raw_v2.1
# PREDICTORS : 
# same as predictors_raw_v2.1 with mofifications on Habitat columns : 
# habitat NA replaced 
# proprotions computed differently --> see Habitat section 
# predictors (147) : area_km2, mpa_fully, canyon_dist_m_min, canyon_dist_m_max, canyon_dist_m_mean, canyon_dist_m_range, canyon_dist_m_weight, port_dist_m_min, port_dist_m_max, port_dist_m_mean, port_dist_m_range, port_dist_m_weight, mpa_dist_m_min, mpa_dist_m_max, mpa_dist_m_mean, mpa_dist_m_range, mpa_dist_m_weight, shore_dist_m_min, shore_dist_m_max, shore_dist_m_mean, shore_dist_m_range, shore_dist_m_weight, canyon_area_km2, canyon_type, canyon_mean_depth, canyon_length, canyon_width, canyon_shape_length, canyon_shape_area, mpa_protection, gravity_mean, gravity_min, gravity_max, gravity_range, aspect_min, aspect_max, aspect_mean, roughness_min, roughness_max, roughness_mean, slope_min, slope_max, slope_mean, tpi_min, tpi_max, tpi_mean, tri_min, tri_max, tri_mean, main_habitat, nb_habitat_per_km2, grouped_main_habitat, grouped_nb_habitat_per_km2, bathy_mean, bathy_range, bathy_min, bathy_max, dist_seabed_depthsampling, wind_min_24h, wind_max_24h, wind_mean_24h, vel_min_24h, vel_max_24h, vel_mean_24h, wind_max_7j, wind_min_7j, wind_mean_7j, vel_max_7j, vel_min_7j, vel_mean_7j, wind_max_1m, wind_min_1m, wind_mean_1m, vel_max_1m, vel_min_1m, vel_mean_1m, ws_max_1y, ws_min_1y, ws_mean_1y, vel_max_1y, vel_min_1y, vel_mean_1y, temp_min_24h, temp_max_24h, temp_mean_24h, sal_min_24h, sal_max_24h, sal_mean_24h, temp_max_7j, temp_min_7j, temp_mean_7j, sal_max_7j, sal_min_7j, sal_mean_7j, temp_max_1m, temp_min_1m, temp_mean_1m, sal_max_1m, sal_min_1m, sal_mean_1m, temp_max_1y, temp_min_1y, temp_mean_1y, sal_max_1y, sal_min_1y, sal_mean_1y, cop_chl_day_mean, cop_chl_day_min, cop_chl_day_max, cop_chl_week_mean, cop_chl_week_min, cop_chl_week_max, cop_chl_month_mean, cop_chl_month_min, cop_chl_month_max, cop_chl_year_mean, cop_chl_year_min, cop_chl_year_max, cop_chl_5years_mean, cop_chl_5years_min, cop_chl_5years_max, cop_analysed_sst_day_mean, cop_analysed_sst_day_min, cop_analysed_sst_day_max, cop_analysed_sst_week_mean, cop_analysed_sst_week_min, cop_analysed_sst_week_max, cop_analysed_sst_month_mean, cop_analysed_sst_month_min, cop_analysed_sst_month_max, cop_analysed_sst_year_mean, cop_analysed_sst_year_min, cop_analysed_sst_year_max, cop_analysed_sst_5years_mean, cop_analysed_sst_5years_min, cop_analysed_sst_5years_max, x, y, replicates, matte_morte_p_oceanica_mean, algues_infralittorales_mean, habitats_artificiels_mean, soft_bottom_mean, meadow_mean, rock_mean, coralligenous_mean, geom
cat(colnames(x), sep = ", ")

pred <- x # (after modifications in the habitat section and in the CLEAN UP BUFF section)

# export
st_write(x, "./data/processed_data/predictors/predictors_raw_v2.2.gpkg", append = FALSE)

# Export predictors_raw_v3.0 -----
# 28/11/2025. 
# ROWS
# Based on mtdt_7.gpkg, same as predictors_raw_v2 : accounting for no detection samples 
# PREDICTORS :  
#  --- Key changes ---
# Boat detection added 
# Habitat including "Zone bathyale"

st_write(x, "./data/processed_data/predictors/predictors_raw_v3.0.gpkg", append = FALSE)



