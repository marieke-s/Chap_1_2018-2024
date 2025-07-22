#------------- Description ---------------------
# Purpose: 
# This script aims to 
#     - make a subset of the full eDNA metadata file to keep samples needed for our study
#     - pool sampling replicates together within a single point
#     - make a buffer around the pooled sampling points


# The data subset resulting from this script will be used for predictors extraction (xxxx.R) and filtering the occurence dataset (xxx.R)

# Data source: 
# A csv file of the medata of mediterranean eDNA samplings in 2018-2024 produced by Laure Velez and Amandine Avouac (amandine.avouac@umontpellier.fr).
# A csv file for the IPOCOM metadata sent by Celia Bertrand on Slack the 7/07/2025 : eREF_IPOCOM TRANSECT 2024.csv
# Author: Marieke Schultz

# Date script created: 2025-07-07
#------------- Setting up ---------------------
# Remove existing objects
rm(list = ls())

setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(concaveman)
library(dplyr)
library(geosphere)
library(ggplot2)
library(igraph)
library(purrr)
library(leaflet)
library(lubridate)
library(readxl)
library(readr)
library(reticulate)
library(sf)
library(stringr)
library(terra)


# Load functions
source("./utils/Fct_Data-Prep.R")



#------------- Load and clean data ---------------------
# Load the metadata file
mtdtfull <- readr::read_csv("./data/raw_data/eDNA/Med_metadonnees_ADNe - v1.1_2018-2024.csv")

# If mdtdt is version <= 1.1 :
# For row with spygen_code = SPY232658 replace sampling_depth to 80
if (any(mtdtfull$spygen_code == "SPY232658")) {
  mtdtfull <- mtdtfull %>%
    mutate(depth_sampling = ifelse(spygen_code == "SPY232658", 80, depth_sampling)) 
}

# For row with spygen_code = SPY232652 replace sampling_depth to 20
if (any(mtdtfull$spygen_code == "SPY232652")) {
  mtdtfull <- mtdtfull %>%
    mutate(depth_sampling = ifelse(spygen_code == "SPY232652", 20, depth_sampling)) 
}

# Load IPOCOM data
ipocom <- readr::read_delim("./data/raw_data/eDNA/eREF_IPOCOM TRANSECT 2024.csv", delim = ";")
  
# Correct IPOCOM coordinates errors
# SPY233779 IPOCOM_135 Y_DEBUT 43.4785833333333 --> 42.4785833333333
# SPY2401466 IPOCOM_13 Y_DEBUT 43.7057833 --> 43.5057833
# SPY2401440 IPOCOM_12 Y_DEBUT 42.53943333 --> 43.53943333
# SPY2401596 IPOCOM_76 X_FIN 5.44973333333333 --> 5.34973333333333

if (any(ipocom$N_Transect == "IPOCOM_135")) {
  ipocom <- ipocom %>%
    mutate(Y_DEBUT = ifelse(N_Transect == "IPOCOM_135", 42.4785833333333, Y_DEBUT))
}
if (any(ipocom$N_Transect == "IPOCOM_13")) {
  ipocom <- ipocom %>%
    mutate(Y_DEBUT = ifelse(N_Transect == "IPOCOM_13", 43.5057833, Y_DEBUT))
}
if (any(ipocom$N_Transect == "IPOCOM_12")) {
  ipocom <- ipocom %>%
    mutate(Y_DEBUT = ifelse(N_Transect == "IPOCOM_12", 43.53943333, Y_DEBUT))
}
if (any(ipocom$N_Transect == "IPOCOM_76")) {
  ipocom <- ipocom %>%
    mutate(X_FIN = ifelse(N_Transect == "IPOCOM_76", 5.34973333333333, X_FIN))
}


# Error on SPY2401568 IPOCOM_62 ????s

# Correct Mtdtfull coordinates errors
# SPY2401010 : longitude_end_DD 42.902999999999999 --> ???




#------------- 0. Add IPOCOM 2024 data -----------------
# Rename ipocom columns to match mtdtfull
ipocom <- ipocom %>%
  rename(
    spygen_code = N_SPYGEN,
    date = DATE,
    latitude_start_DD = Y_DEBUT,
    longitude_start_DD = X_DEBUT,
    latitude_end_DD = Y_FIN,
    longitude_end_DD = X_FIN,
    subsite = N_Transect,
    time_start = HEURE_DEBUT,
    duration = TEMPS_FILTRATION,
    comments = COMMENTAIRE
  ) 

# Format ipocom date as mtdtfull date
ipocom$date <- as.Date(ipocom$date, format = "%d/%m/%Y")

# Format duration column 
ipocom$duration <- as.numeric(ipocom$duration) / 60

# Remove ipocom columns not in mtdtfull
ipocom <- ipocom %>%
  select(-c("HEURE_FIN", "Temp"))

# Add missing columns to ipocom with NA
missing_cols <- setdiff(names(mtdtfull), names(ipocom))

for (col in missing_cols) {
  ipocom[[col]] <- NA
}

# For column 'country' set "France" for all rows
ipocom$country <- "France"

# For column 'component' set "open_ocean" for all rows
ipocom$pool <- "no"

# Reorder ipocom columns to match mtdtfull
ipocom <- ipocom[, names(ipocom)]

# Combine the two datasets
mtdtcomb <- bind_rows(mtdtfull, ipocom)

# Clean up
rm(mtdtfull, ipocom, col, missing_cols)






#------------- 1. Subset n°1 : France + marine + coastal --------------------
# Explanation:
# Selection of points sampled in France 
# Selection of marine points, i.e. outside harbours, lagoons, ports, rivers, estuaries 
# Filtering out the seamont and open_ocean samples because they are outside the coastal zones of which we are interested in 
# = 377 samples removed

mtdt_1 <- mtdtcomb %>%
  filter(is.na(country) | country == "France") %>%
  filter(
    is.na(component) | !(component %in% c("harbour", 
                                          "lagoon", 
                                          "port", 
                                          "freshwater_river",
                                          "estuary", 
                                          "open_ocean", 
                                          "seamount")))










#------------- 2. Pool 15L pooled samples ------------
# Explanation : There is one protocol that filters simultaneously 2 15L samples on the submersible pump. The 2 samples are then pooled by SpyGen and analysed as an unique 30L sample. There are indicated in the column "pool". 

# Filter samples that need to be pooled
mtdt_pool <- mtdt_1 |>
  dplyr::filter(!(pool == "no"))

# we expect to reduce the number of mtdt_1 rows by 162/2 = 81 rows 

# ---- Automatically pool samples using the `pool` column 
mtdt_2 <- mtdt_1  # initialize

for (p in unique(mtdt_pool$pool)) {
  
  codes <- unlist(strsplit(p, split = "_"))
  new_code <- p
  
  mtdt_2 <- combine_rows(
    data = mtdt_2,  # ✅ keep updating the pooled dataset
    codes_to_combine = codes,
    new_code = new_code,
    proceed = TRUE
  )
}


# Clean up 
rm(mtdt_pool, p, codes, new_code)










#------------- 3. Pool replicates Amandine -------






#------------- 3. Pool replicates Marieke --------------------
# Explanation: In most protocols, two samples were done at the same time and same location : they are field replicates. For the future analyses, we want to pool these replicates together. Problem is, there are no replicate values saved in the metadata file. This section aims to recover which samples go together. 

#---- Extract SHOM bathy -----
# Explanation: 
# Since there are NA values in the mtdt depth column, we use the bathymetry extracted from a digital surface model (MNT bathymétrique du SHOM: https://diffusion.shom.fr/donnees/bathymerie/mnt-facade-gdl-ca-homonim.html) to fill NA values. 


# Dataset 
dt <- mtdt_2

# Load bathy
bathy <- terra::rast("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")

# Extract SHOM bathymetry from start/end coordinates
start_bathy = abs(terra::extract(bathy, terra::vect(dt, geom = c("longitude_start_DD", "latitude_start_DD"), crs = crs(bathy)))[, 2])
end_bathy   = abs(terra::extract(bathy, terra::vect(dt, geom = c("longitude_end_DD", "latitude_end_DD"), crs = crs(bathy)))[, 2])

dt <- dt |>
  mutate(
    max_shom_bathy = pmax(start_bathy, end_bathy, na.rm = TRUE), # The maximum bathymetry between start and end points is retained
    combined_bathy = coalesce(as.numeric(depth_seafloor), max_shom_bathy)
  ) 

# Correlation between depth_seafloor and max_shom_bathy
cor(dt$depth_seafloor, dt$max_shom_bathy, use = "pairwise.complete.obs")

# Assign back dataset name
mtdt_2 <- dt
rm(dt, bathy, end_bathy, start_bathy)











#---- Buffer transect --------------------
buff <- buffer_transect(
  df = mtdt_2,
  start_lon_col = "longitude_start_DD",
  start_lat_col = "latitude_start_DD",
  end_lon_col = "longitude_end_DD",
  end_lat_col = "latitude_end_DD",
  buffer_dist = 500 # buffer distance is meters
)

beepr::beep()

# Plot buffers
plot(buff, max.plot = 1)

# Check if buff are valid
sum(!st_is_valid(buff))

# Save buff
st_write(buff, "./data/processed_data/buff_500m_2023_FR_surface_coastal_30L.gpkg", delete_dsn = TRUE)



#---- Assign replicates with subsite ----

# STEP 1 — Identify samples with known replicate-defining subsites
known_rep_replicates <- buff %>%
  mutate(
    # Normalize (aller/retour) to common subsite key
    replicate_key = case_when(
      str_detect(subsite, "^piaf_\\d+$") ~ subsite,
      str_detect(subsite, "^t_\\d+$") ~ subsite,
      str_detect(subsite, "^IPOCOM_\\d+$") ~ subsite,
      str_detect(subsite, "^\\d+_\\((aller|retour)\\)$") ~ str_replace(subsite, "_\\((aller|retour)\\)", ""),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(replicate_key)) %>%
  group_by(replicate_key) %>%
  mutate(replicate = cur_group_id()) %>%
  ungroup() %>%
  select(-replicate_key)

# STEP 2 — Separate out samples not already assigned to replicates
remaining_to_group <- buff %>%
  filter(!spygen_code %in% known_rep_replicates$spygen_code)






#------------- Assign replicates with buffer MISTRAL ----------------
# Load necessary libraries
library(sf)

# Assuming 'remaining_to_group' is your data frame and it has a column 'wkt_geometry' with WKT strings
# Ensure that 'wkt_geometry' is a character vector of WKT strings
remaining_to_group$wkt_geometry <- as.character(remaining_to_group$wkt_geometry)

# Convert the WKT geometry column to an sfc object
geometry_sfc <- st_as_sfc(remaining_to_group$wkt_geometry, crs = 4326) # Assuming WGS84 coordinate system

# Create a new data frame excluding the original geometry column
remaining_to_group_sf <- remaining_to_group
remaining_to_group_sf$wkt_geometry <- NULL

# Convert the data frame to an sf object with the correct geometry column
sf_data <- st_as_sf(remaining_to_group_sf, sf_column_name = "geometry_sfc")

# Perform spatial overlap analysis
overlaps <- st_intersects(sf_data, by_element = TRUE)

# Create a list to store groups of overlapping geometries
overlap_groups <- list()

# Iterate over each geometry to find overlaps
for (i in 1:nrow(sf_data)) {
  # Find which geometries overlap with the current geometry
  overlapping_indices <- which(overlaps[i, ])
  
  # Add the current geometry's index to the overlapping indices
  overlapping_indices <- c(i, overlapping_indices)
  
  # Sort the indices to have a consistent grouping
  overlapping_indices <- sort(unique(overlapping_indices))
  
  # Create a key for the group
  group_key <- paste(sort(overlapping_indices), collapse = "_")
  
  # Store the group in the list
  overlap_groups[[group_key]] <- overlapping_indices
}

# Create a replicate index based on the overlap groups
replicate_index <- integer(nrow(sf_data))
names(replicate_index) <- 1:nrow(sf_data)

for (group in overlap_groups) {
  replicate_index[as.character(group)] <- match(paste(sort(group), collapse = "_"), names(overlap_groups))
}

# Add the replicate index to the original data frame
remaining_to_group$replicate <- replicate_index

# View the updated data frame with the replicate index
head(remaining_to_group)




# # STEP 3 — Apply buffer-based grouping to only the remaining samples
# grouped_remaining <- remaining_to_group %>%
#   mutate(date = as.Date(date)) %>%
#   {
#     overlaps <- st_intersects(.)
#     g <- graph_from_adj_list(overlaps)
#     .$buffer_group <- components(g)$membership
#     .
#   } %>%
#   group_by(buffer_group) %>%
#   mutate(n_in_group = n()) %>%
#   ungroup() %>%
#   mutate(overlap_date_group = if_else(n_in_group > 2, paste(buffer_group, date, sep = "_"), as.character(buffer_group))) %>%
#   select(-buffer_group, -n_in_group) %>%
#   group_by(overlap_date_group) %>%
#   group_split() %>%
#   lapply(function(df_group) {
#     if (nrow(df_group) > 2) {
#       df_group %>%
#         group_by(max_shom_bathy) %>%
#         mutate(bathy_sub_id = cur_group_id()) %>%
#         ungroup() %>%
#         mutate(
#           bathy_sub_id = dense_rank(bathy_sub_id),
#           bathy_subgroup = paste0(overlap_date_group[1], "_b", bathy_sub_id)
#         ) %>%
#         select(-bathy_sub_id)
#     } else {
#       df_group$bathy_subgroup <- df_group$overlap_date_group
#       df_group
#     }
#   }) %>%
#   bind_rows() %>%
#   group_by(bathy_subgroup) %>%
#   mutate(replicate = max(known_rep_replicates$replicate, na.rm = TRUE) + cur_group_id()) %>%  # shift replicate IDs to avoid overlap
#   ungroup() %>%
#   select(-overlap_date_group, -bathy_subgroup)

# # STEP 3 — Simplified buffer-based grouping using st_intersection
# 
# # 3.1 Ensure valid geometries and CRS consistency
# remaining_to_group <- remaining_to_group %>%
#   st_make_valid() %>%
#   st_buffer(0)
# st_crs(remaining_to_group) <- st_crs(buff)
# 
# # 3.2 Prepare for spatial grouping
# remaining_to_group <- remaining_to_group %>%
#   mutate(date = as.Date(date),
#          id = row_number())
# 
# # 3.3 Compute intersections
# intersection_df <- st_intersection(
#   remaining_to_group %>% select(id, geometry),
#   remaining_to_group %>% select(id, geometry)
# ) %>%
#   filter(id != id.1) %>%
#   select(id, id.1)
# 
# # 3.4 Build adjacency graph and get group memberships
# g <- igraph::graph_from_data_frame(intersection_df, directed = FALSE)
# membership <- igraph::components(g)$membership
# remaining_to_group$buffer_group <- membership[as.character(remaining_to_group$id)]
# 
# # 3.5 Assign replicate IDs with optional bathymetry stratification
# grouped_remaining <- remaining_to_group %>%
#   group_by(buffer_group) %>%
#   mutate(n_in_group = n()) %>%
#   ungroup() %>%
#   mutate(overlap_date_group = if_else(n_in_group > 2, paste(buffer_group, date, sep = "_"), as.character(buffer_group))) %>%
#   select(-n_in_group) %>%
#   group_by(overlap_date_group) %>%
#   group_split() %>%
#   lapply(function(df_group) {
#     if (nrow(df_group) > 2) {
#       df_group %>%
#         group_by(bathy_band = cut(df_group$max_shom_bathy, breaks = pretty(df_group$max_shom_bathy, n = 3), include.lowest = TRUE)) %>%
#         mutate(
#           bathy_sub_id = cur_group_id(),
#           bathy_subgroup = paste0(overlap_date_group[1], "_b", bathy_sub_id)
#         ) %>%
#         ungroup() %>%
#         select(-bathy_sub_id, -bathy_band)
#     } else {
#       df_group$bathy_subgroup <- df_group$overlap_date_group
#       df_group
#     }
#   }) %>%
#   bind_rows() %>%
#   group_by(bathy_subgroup) %>%
#   mutate(replicate = max(known_rep_replicates$replicate, na.rm = TRUE) + dense_rank(cur_group_id())) %>%
#   ungroup() %>%
#   select(-overlap_date_group, -bathy_subgroup, -id)









# # STEP 3 — Simplified buffer-based grouping using st_intersection
# 
# # 3.1 Ensure valid geometries and CRS consistency
# remaining_to_group <- remaining_to_group %>%
#   st_make_valid() %>%
#   st_buffer(0)
# st_crs(remaining_to_group) <- st_crs(buff)
# 
# # 3.2 Prepare for spatial grouping
# remaining_to_group <- remaining_to_group %>%
#   mutate(date = as.Date(date),
#          id = row_number())
# 
# 
# # Detect the geometry column name
# geom_col <- attr(remaining_to_group, "sf_column")
# 
# # 3.3 Compute intersections — robust to geometry column name
# sf1 <- remaining_to_group %>% select(id, geometry = all_of(geom_col))
# sf2 <- remaining_to_group %>% select(id_other = id, geometry = all_of(geom_col))
# 
# # Perform spatial intersection
# intersection_df <- st_intersection(sf1, sf2) %>%
#   filter(id != id_other) %>%
#   st_drop_geometry() %>%
#   select(id, id_other)
# 
# 
# # 3.4 Build adjacency graph and get group memberships
# g <- igraph::graph_from_data_frame(intersection_df, directed = FALSE)
# membership <- igraph::components(g)$membership
# remaining_to_group$buffer_group <- membership[as.character(remaining_to_group$id)]
# 
# # 3.5 Assign replicate IDs with optional bathymetry stratification
# grouped_remaining <- remaining_to_group %>%
#   group_by(buffer_group) %>%
#   mutate(n_in_group = n()) %>%
#   ungroup() %>%
#   mutate(overlap_date_group = if_else(n_in_group > 2, paste(buffer_group, date, sep = "_"), as.character(buffer_group))) %>%
#   select(-n_in_group) %>%
#   group_by(overlap_date_group) %>%
#   group_split() %>%
#   lapply(function(df_group) {
#     if (nrow(df_group) > 2) {
#       df_group %>%
#         group_by(bathy_band = cut(df_group$max_shom_bathy, breaks = pretty(df_group$max_shom_bathy, n = 3), include.lowest = TRUE)) %>%
#         mutate(
#           bathy_sub_id = cur_group_id(),
#           bathy_subgroup = paste0(overlap_date_group[1], "_b", bathy_sub_id)
#         ) %>%
#         ungroup() %>%
#         select(-bathy_sub_id, -bathy_band)
#     } else {
#       df_group$bathy_subgroup <- df_group$overlap_date_group
#       df_group
#     }
#   }) %>%
#   bind_rows() %>%
#   group_by(bathy_subgroup) %>%
#   mutate(replicate = max(known_rep_replicates$replicate, na.rm = TRUE) + dense_rank(cur_group_id())) %>%
#   ungroup() %>%
#   select(-overlap_date_group, -bathy_subgroup, -id)


# STEP 4 — Combine both groups
final_grouped <- bind_rows(known_rep_replicates, grouped_remaining)

# STEP 5 — Manual override
final_grouped <- final_grouped %>%
  mutate(
    replicate = case_when(
      spygen_code %in% c("SPY232834", "SPY232833") ~ min(replicate[spygen_code %in% c("SPY232834", "SPY232833")]),
      TRUE ~ replicate
    )
  )


#---- Check results ----
final_grouped %>%
  st_drop_geometry() %>%
  count(replicate) %>%
  count(n, name = "num_groups") %>%
  arrange(n)

# Add n_in_replicate column
final_grouped <- final_grouped %>%
  group_by(replicate) %>%
  mutate(n_in_replicate = n()) %>%
  ungroup()


# Save as geopackage
st_write(final_grouped, "./data/processed_data/grouped_buffer_1.gpkg")

# #---- Assign replicates with time and date ----
# # 1. Split replicate groups by size
# rep_groups_split <- final_grouped %>%
#   st_drop_geometry() %>%
#   count(replicate, name = "group_size")
# 
# # 2. Join size info to original data
# final_grouped <- final_grouped %>%
#   left_join(rep_groups_split, by = "replicate")
# 
# # 3. Separate groups needing refinement (size > 2)
# to_refine <- final_grouped %>% filter(group_size > 2)
# keep_as_is <- final_grouped %>% filter(group_size <= 2)
# 
# # 4. Refine over-large groups using date + time_start
# refined <- to_refine %>%
#   group_by(replicate, date, time_start) %>%
#   mutate(new_replicate = cur_group_id()) %>%
#   ungroup()
# 
# # 5. Offset new replicate IDs to avoid overlap
# offset <- max(final_grouped$replicate, na.rm = TRUE)
# refined <- refined %>%
#   mutate(replicate = new_replicate + offset) %>%
#   select(-new_replicate, -group_size)
# 
# # 6. Recombine all together
# final_grouped_refined <- bind_rows(
#   keep_as_is %>% select(-group_size),
#   refined
# )

#---- Assign replicates with time and date ----

# 1. Split replicate groups by size
rep_groups_split <- final_grouped %>%
  st_drop_geometry() %>%
  count(replicate, name = "group_size")

# 2. Join size info to original data
final_grouped <- final_grouped %>%
  left_join(rep_groups_split, by = "replicate")

# 3. Separate groups to refine (where group_size ≠ 2)
to_refine <- final_grouped %>% filter(group_size != 2)
keep_as_is <- final_grouped %>% filter(group_size == 2)

# 4. Refine non-standard groups using date + time_start
refined <- to_refine %>%
  group_by(replicate, date, time_start) %>%
  mutate(new_replicate = cur_group_id()) %>%
  ungroup()

# 5. Offset new replicate IDs to avoid overlap
offset <- max(final_grouped$replicate, na.rm = TRUE)
refined <- refined %>%
  mutate(replicate = new_replicate + offset) %>%
  select(-new_replicate, -group_size)

# 6. Recombine all together
final_grouped_refined <- bind_rows(
  keep_as_is %>% select(-group_size),
  refined
)











#---- Check results ----
final_grouped_refined %>%
  st_drop_geometry() %>%
  count(replicate) %>%
  count(n, name = "num_groups") %>%
  arrange(n)


grouped_buffer_2 <- final_grouped_refined %>%
  group_by(replicate) %>%
  mutate(n_in_replicate = n()) %>%
  ungroup()

# save as gpkg
st_write(grouped_buffer_2, "./data/processed_data/grouped_buffer_2.gpkg", delete_dsn = TRUE)



#---- Plot results ----
t <- final_grouped %>% filter(n_in_replicate == 4)



ggplot(t) +
  geom_sf(aes(fill = factor(n_in_replicate)), color = "black", size = 0.2) +
  scale_fill_manual(
    values = c("1" = "red", "2" = "#91cf60", "3" = "#d9ef8b", "4" = "#fee08b", "5" = "pink", "6" = "blue"),
    na.value = "grey80",
    name = "Samples per replicate",
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(title = "Replicate Groups by Sample Count",
       subtitle = "Red = singleton replicates",
       fill = "Samples") +
  theme_minimal()











#---- Check the result ------------
# count number of samples per replicate group
buffered_grouped %>%
  st_drop_geometry() %>%               
  count(replicate) %>%            
  count(n, name = "num_groups") %>%    
  arrange(n)                           


buffered_grouped <- buffered_grouped %>%
  group_by(replicate) %>%
  mutate(n_in_replicate = n()) %>%
  ungroup()

ggplot(buffered_grouped) +
  geom_sf(aes(fill = factor(n_in_replicate)), color = "black", size = 0.2) +
  scale_fill_manual(
    values = c("1" = "red", "2" = "#91cf60", "3" = "#d9ef8b", "4" = "#fee08b", "5" = "pink", "6" = "blue"),
    na.value = "grey80",
    name = "Samples per replicate",
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(title = "Replicate Groups by Sample Count",
       subtitle = "Red = singleton replicates",
       fill = "Samples") +
  theme_minimal()



# Clean up ----
dev.off()
df <- buffered_grouped %>%
  select(-c("wkt_geometry", "n_in_replicate"))
rm(buff, mtdt_2023_FR_surface_coastal_30L, buffered_grouped)



























#------------- 3. Subset n°2 : 40m depth --------------------
#------------- 4. Make merged buffer ---------------------





#--------------------------- Boat detection Luka------------
# Load data
boat <- readr::read_csv("./data/raw_data/predictors/Boats/predictions_Luka_20250714.csv")

# Convert coordinates 
boat <- boat %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 32632) %>%
  st_transform(crs = 4326)

# Filter data of interest
boats <- boat %>% 
  filter(type =="boat") %>% 
  filter(confidence > 0.5) 


#--- Check correlation between length and confidence ---

# Check for normality on a sample 
shapiro.test(sample(boat$length, 5000))
shapiro.test(sample(boat$confidence, 5000))

# Check correlation
cor(boat$length, boat$confidence, method = "pearson")

#--- 



# Select all dates starting with "202308"
boats <- boats %>%
  filter(str_starts(date, "202308"))

# Save as gpkg
st_write(boats, "./data/processed_data/boats_Luka_202308.gpkg", delete_dsn = TRUE)


#--------------------------- Check replicate column Amandine -----------
dt <- readr::read_csv("./data/raw_data/eDNA/Med_metadonnees_ADNe_2018_2024_V4.csv")
dt$time_start <- as.character(dt$time_start)

years <- c(2018,2019,2020,2021,2022,2023,2024)

for (year in years) {
  # Filter for the current year
  dt_year <- dt %>%
    filter(year(date) == year) %>%
    select(-c("lockdown", "BiodivMed2023","Tele01","Pleo","Mamm01", "Vert01","X16s_Metazoa","Bact02","Euka02" ))%>%
    filter(is.na(country) | country == "France") %>%
    filter(
      is.na(component) | !(component %in% c("harbour", 
                                            "lagoon", 
                                            "port", 
                                            "freshwater_river",
                                            "estuary", 
                                            "open_ocean", 
                                            "seamount")))
    
  
  # create buffer
  buff <- buffer_transect(
    df = dt_year,
    start_lon_col = "longitude_start_DD",
    start_lat_col = "latitude_start_DD",
    end_lon_col = "longitude_end_DD",
    end_lat_col = "latitude_end_DD",
    buffer_dist = 500 # buffer distance is meters
  )
  
  # save as gpkg
  st_write(buff, paste0("./data/processed_data/buff_Amandine_", year, ".gpkg"), delete_dsn = TRUE)
}

#---- Buffer transect --------------------
buff <- buffer_transect(
  df = dt,
  start_lon_col = "longitude_start_DD",
  start_lat_col = "latitude_start_DD",
  end_lon_col = "longitude_end_DD",
  end_lat_col = "latitude_end_DD",
  buffer_dist = 500 # buffer distance is meters
)

beepr::beep()

# Plot buffers
plot(buff, max.plot = 1)

# Check if buff are valid
sum(!st_is_valid(buff))

# Save buff
st_write(buff, "./data/processed_data/replicates_Amandine", delete_dsn = TRUE)
