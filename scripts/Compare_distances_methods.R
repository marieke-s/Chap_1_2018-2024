#------------- Description ---------------------
# Purpose: 
# Compare distances methods as computed by Martin and by other sources.

# Author: Marieke Schultz

# Contributors : Martin Paquet, Pauline Viguier, Laure Velez, Marie Orblin

# Date script created: 6/11/2025

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
# Documentation port columns : https://services.sandre.eaufrance.fr/telechargement/geo/PTS/sandre_dictionnaire_PTS_2.pdf
# Load data ----
dist <- st_read("./data/raw_data/predictors/Distances/buffer_with_closest_feats_search_outlier_treshold_shore_50m.gpkg")
dist <- as.data.frame(dist)

# Clean [ ASK MARTIN !!!!! ]----

dist <- dist %>%
  
  dplyr::select(-c ("ID","date","time_start","depth_sampling","depth_seafloor","lockdown","BiodivMed2023","method","country","region","site","subsite","component","habitat","protection","project", "Tele01","Pleo","Mamm01","Vert01",  "Tele01" ,"Pleo","Mamm01","Vert01","X16s_Metazoa","Bact02" ,"Euka02" ,"duration_total","comments","estimated_volume_total", "mpa_name...17"  )) %>%  # Remove mtdt cols
  
  dplyr::select(-"canyon_Geomorphic") %>% # because gives no info : unique(dist$canyon_Geomorphic) # "Canyon"
  dplyr::select(-"canyon_Ocean") %>% # because all "Mediterranean Sea"
  dplyr::select(-c("port_NumTexteRe", "port_DatePubliT", "port_CdTypeText", "port_MnTypeText", "port_URLTexteRe")) %>% # full NA columns
  dplyr::select(-c("port_CdZonePort", "port_MnTypeZone", "port_DtCreatZon", "port_DtMajZoneP", "port_ComZonePor", "port_CdPort", "port_StatutZone", "port_CdTypeZone", "port_gid", "port_aire_m2")) %>% # We delete metadata associated to ports because incomplete
  rename("mpa_name" = "mpa_name...81") %>%
  dplyr::select(-c("canyon_dist_min_from_buff_m", "port_dist_min_from_buff_m", "mpa_dist_min_from_buff_m")) # These cols were used for methodological checks only ([MARTIN CHECK])calcul de la distance du accCost depuis le buffer à la géométrie de variable la plus proche)


# Count NA per column in dist
sapply(dist, function(x) sum(is.na(x)))




# Q Martin 
# canyon_OBJECTIID
# canynon_Delta_D









# Merge -----
# Detect columns in dist that are not in buff
extra_cols <- setdiff(names(dist), names(buff))

# left_join keeps all rows in buff and brings matching columns from dist by 'replicates'
buff <- buff %>%
  left_join(dist %>% dplyr::select(replicates, all_of(extra_cols)), by = "replicates")

rm(dist, extra_cols)















###### METHODOLIGAL CHECKS : Comparison distance to shore Martin VS Euclidian centroid distance -----
# Load data ############
# Explanation : distance to shore is computed from the coastline shapefile using '2.1_Distance_shore.py' from orignal code of Pauline Viguier. 
# Important methodological considerations : 
# Distance to shore is computed from replicates group's buffer centroids with euclidian distance. 
# When a centroid is on land (which happened for ~ 67/792 centroids) we set the distance to shore to 1 meter because it was actually not sampled on land (obviously) and it ended up there only because we simplified transects into straight lines between transect start and end points.

# Run 2.1_Distance_shore.py to compute distance to shore with euclidian buffer centroid distance.
system("python3 ./scripts/2_Euclidian_Distance_shore.py")

# Load buff with dist_shore (eculidiand) 
dist_shore <- st_read("./data/processed_data/predictors/mtdt_5_dist-shore.gpkg")


# Add dist_shore$dist_shore_m to buff
df <- dist_shore %>%
  as.data.frame() %>%
  dplyr::select(replicates, dist_shore_m) %>%
  left_join(
    buff %>% as.data.frame(),
    by = "replicates"
  )

# Comparison ----



# Correlation df$shore_dist_m_weight (accCost) vs df$dist_shore_m (euclidian)
cor_pearson <- cor(df$shore_dist_m_weight, df$dist_shore_m, use = "complete.obs", method = "pearson")

# Plot 
r2 <- cor_pearson^2
ggplot(data = NULL, aes(x = df$shore_dist_m_weight, y = df$dist_shore_m)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  labs(
    x = "Weighted mean distance to Shore accCost ",
    y = "Distance to Shore euclidian from buffer centroid",
    title = "Relationship between Port Distance Metrics",
    subtitle = sprintf("Pearson r = %.3f | R² = %.3f", cor_pearson, r2)
  ) +
  theme_minimal(base_size = 14)




# Correlation df$shore_dist_m_min (accCost) vs df$dist_shore_m (euclidian)
cor_pearson <- cor(df$shore_dist_m_min, df$dist_shore_m, use = "complete.obs", method = "pearson")

# Plot 
r2 <- cor_pearson^2
ggplot(data = NULL, aes(x = df$shore_dist_m_min, y = df$dist_shore_m)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  labs(
    x = "min distance to Shore accCost ",
    y = "Distance to Shore euclidian from buffer centroid",
    title = "Relationship between Port Distance Metrics",
    subtitle = sprintf("Pearson r = %.3f | R² = %.3f", cor_pearson, r2)
  ) +
  theme_minimal(base_size = 14)




# Clean
rm(df, dist_shore, cor_pearson, r2)
# dev.off() # close plots









###### METHODOLIGAL CHECKS : Comparison distance to reserve Martin VS in_out Laure ############

# Min dist ----

# Make column reserve_martin : if $mpa_dist_m_min = 0 then reserve_martin = 1 else reserve_martin = 0
buff$reserve_martin <- ifelse(buff$mpa_dist_m_min == 0, 1, 0)

# Compare buff$mpa_fully with buff$reserve_martin
comparison <- buff %>%
  dplyr::select(replicates, mpa_fully, reserve_martin) %>%
  mutate(match = ifelse(mpa_fully == reserve_martin, TRUE, FALSE))

# Count matches and mismatches
comparison %>%
  as.data.frame() %>%
  group_by(match) %>%
  summarise(count = n()) %>%
  print()

61/788 *100 # 7.7% mismatch avant jointure ? 
63/788 # 7.9% mismatch après jointure ?
61/788 *100 # 7.7% mismatch à 50m

# Add a column "Who_is_1", when match = TRUE = NA, when match = FALSE and mpa_fully = 1 then "Laure", when match = FALSE and reserve_martin = 1 then "Martin"
comparison <- comparison %>%
  mutate(
    Who_is_1 = case_when(
      match == TRUE ~ NA_character_,
      match == FALSE & mpa_fully == 1 ~ "Laure",
      match == FALSE & reserve_martin == 1 ~ "Martin"
    )
  )

# [OPTIONAL] Export to check in QGIS
st_write(comparison, "./data/processed_data/predictors/MPA/check_reserve_match_50m.gpkg", delete_dsn = TRUE)

#---
# All mismatch are due to Laure = 0 and Martin = 1. This is because Martins's distances to reserve are computed on buffer while Laure is checking each sample point. Thus, close but not inside transects, whn buffer touch the reserve and are thus considered by Martin as "in reserve".


# Max dist ----

# We check if this remains true when using the "max distance to reserve" from Martin instead of the "min distance to reserve".
buff$reserve_martin_max <- ifelse(buff$mpa_dist_m_max == 0, 1, 0)

comparison <- buff %>%
  dplyr::select(replicates, mpa_fully, reserve_martin_max) %>%
  mutate(match = ifelse(mpa_fully == reserve_martin_max, TRUE, FALSE))

comparison %>%
  as.data.frame() %>%
  group_by(match) %>%
  summarise(count = n()) %>%
  print()

105/788 *100 # 13.32% mismatch  

comparison <- comparison %>%
  mutate(
    Who_is_1 = case_when(
      match == TRUE ~ NA_character_,
      match == FALSE & mpa_fully == 1 ~ "Laure",
      match == FALSE & reserve_martin_max == 1 ~ "Martin"
    )
  )

#---
# In thas case mismatches are all due to Laure = 1 and Martin = 0. 


# Mean dist ----

# We check with "weighted mean distance to reserve"
buff$reserve_martin_mean <- ifelse(buff$mpa_dist_m_weight == 0, 1, 0)
comparison <- buff %>%
  dplyr::select(replicates, mpa_fully, reserve_martin_mean, mpa_dist_m_weight) %>%
  mutate(match = ifelse(mpa_fully == reserve_martin_mean, TRUE, FALSE))

comparison %>%
  as.data.frame() %>%
  group_by(match) %>%
  summarise(count = n()) %>%
  print()

105/788 *100 # 13.32% mismatch

comparison <- comparison %>%
  mutate(
    Who_is_1 = case_when(
      match == TRUE ~ NA_character_,
      match == FALSE & mpa_fully == 1 ~ "Laure",
      match == FALSE & reserve_martin_mean == 1 ~ "Martin"
    )
  )


comparison %>%
  as.data.frame() %>%
  dplyr::filter(match == FALSE) %>%
  print()

# [OPTIONAL] Export to check in QGIS
st_write(comparison, "./data/processed_data/predictors/MPA/check_reserve_match_mean_50m.gpkg", delete_dsn = TRUE)

#---
# In thas case mismatches are all due to Laure = 1 and Martin = 0. 




# Clean

rm(comparison)

# remove buff$reserve_martin
buff <- buff %>% dplyr::select(-reserve_martin)


###### METHODOLIGAL CHECKS: Comparison distance to port Martin VS distance to port GFW ############



# Extraction GFW ----

# Data :
# From GFW

# Parameters
var <- "dist_port"
rast <- terra::rast("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1/data/raw_data/predictors/Dist_port/distance-from-port-v1.tiff")
rast <- terra::project(rast, "EPSG:4326")
poly <- buff

# Extraction
df <- spatial_extraction(var = var, rast = rast, poly = buff)  




# Comparison ----

# Compare dff$port_dist_m_weight with df$dist_port_mean with pearson correlation
# Compute Pearson correlation
cor(df$port_dist_m_mean, df$dist_port_mean, use = "complete.obs", method = "pearson") # 50m: 0.6154648
cor(df$port_dist_m_weight, df$dist_port_mean, use = "complete.obs", method = "pearson") # 50m : 0.6153183
cor(df$port_dist_m_min, df$dist_port_min, use = "complete.obs", method = "pearson") # 50m : 0.6142686
cor(df$port_dist_m_max, df$dist_port_max, use = "complete.obs", method = "pearson") # 50m : 0.6254746


# euclid martin MIN vs GFW MIN
cor_pearson <- cor(df$dist_port_min, df$port_dist_euclid_m_min, use = "complete.obs", method = "pearson")
r2 <- cor_pearson^2
ggplot(data = NULL, aes(x = df$dist_port_min, y = df$port_dist_euclid_m_min)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  labs(
    x = "Distance to Port GFW ",
    y = "Distance to Port Martin",
    title = "Relationship between Port Distance Metrics",
    subtitle = sprintf("Pearson r = %.3f | R² = %.3f", cor_pearson, r2)
  ) +
  theme_minimal(base_size = 14)





# euclid martin MIN vs martin MIN
cor_pearson <- cor(df$port_dist_m_min, df$port_dist_euclid_m_min, use = "complete.obs", method = "pearson")
r2 <- cor_pearson^2
ggplot(data = NULL, aes(x = df$port_dist_m_min, y = df$port_dist_euclid_m_min)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  labs(
    x = "Distance to Port Martin ",
    y = "Distance to Port Martin euclid ",
    title = "Relationship between Port Distance Metrics",
    subtitle = sprintf("Pearson r = %.3f | R² = %.3f", cor_pearson, r2)
  ) +
  theme_minimal(base_size = 14)


# Clean
rm(df, cor_pearson, r2, var, rast, poly)

# Potential reasons for differences between GFW and Martin's distances to port :

# raw data port 
# algo distance (euclid or accCost)
# spatial resolution









