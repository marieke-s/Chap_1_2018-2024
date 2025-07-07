#------------- Description ---------------------
# Purpose: 
# This script aims to 
#     - pool sampling replicates together within a single point
#     - make a buffer around the pooled sampling points
#     - make a subset of the full eDNA metadata file to keep samples needed for our study

# The data subset resulting from this script will be used for predictors extraction (xxxx.R) and filtering the occurence dataset (xxx.R)

# Data source: 
# A csv file of the medata of mediterranean eDNA samplings in 2018-2024 produced by Laure Velez and Amandine Avouac (amandine.avouac@umontpellier.fr).

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

plot = "YES"


#------------- Load data ---------------------
# Load the metadata file
mtdtfull <- readr::read_csv("./data/raw_data/eDNA/Med_metadonnees_ADNe - v1.1_2018-2024.csv")


#------------- 1. Subset n°1 : France + marine + 30L --------------------
# Explanation:
# Selection of points sampled in France (90 samples removed)
# Selection of marine points, i.e. outside harbours, lagoons, ports, rivers, estuaries (213 samples removed)
# Filtering out the seamont and open_ocean samples because they are outside the coastal zones of which we are interested in (= 52 samples removed) 
# Selection of 30L samples to unsure consistence through the protocols (There are 209 15L samples, 214 60L samples, 20 45L samples, 2 40L samples and 104 NA ≈ 550 samples removed)

mtdt_1 <- mtdtfull |>
  filter(country == "France") |>  # Filter France
  filter(!(component == "harbour"|
           component =="lagoon"|
           component == "port"|
           component == "freshwater_river"|
           component == "estuary" |
           component == "open_ocean" |
           component == "seamount")) |> 
  filter(estimated_volume == "15" |
           estimated_volume == "30")


#------------- 2. Pool replicates --------------------


#---- Extract SHOM bathy because some depth_seafloor are NA -----
# Explanation: 
# Since there are NA values in the mtdt depth column, we use the bathymetry extracted from a digital surface model (MNT bathymétrique du SHOM: https://diffusion.shom.fr/donnees/bathymerie/mnt-facade-gdl-ca-homonim.html) to fill NA values. 


# Dataset 
dt <- mtdt_1

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
mtdt_1 <- dt
rm(dt, bathy, end_bathy, start_bathy)











#---- Buffer transect --------------------
buff <- buffer_transect(
  df = mtdt_1,
  start_lon_col = "longitude_start_DD",
  start_lat_col = "latitude_start_DD",
  end_lon_col = "longitude_end_DD",
  end_lat_col = "latitude_end_DD",
  buffer_dist = 500 # buffer distance is meters
)

# beepr::beep()

# Plot buffers
# plot(buff, max.plot = 1)

# Check if buff are valid
sum(!st_is_valid(buff))



#---- Assign replicates --------------------
buffered_grouped <- buff %>%
  # Step 1: Group samples by overlapping buffers
  mutate(date = as.Date(date)) %>%
  {
    overlap_matrix <- st_intersects(.)
    g <- igraph::graph_from_adj_list(overlap_matrix)
    .$buffer_group <- igraph::components(g)$membership
    .
  } %>%
  group_by(buffer_group) %>%
  mutate(n_in_group = n()) %>%
  ungroup() %>%
  mutate(
    overlap_date_group = if_else(
      n_in_group > 2,
      paste(buffer_group, date, sep = "_"),
      as.character(buffer_group)
    )
  ) %>%
  select(-buffer_group, -n_in_group) %>%
  
  # Step 2: Apply bathy-based subgroups to groups > 2 samples
  group_by(overlap_date_group) %>%
  group_split() %>%
  lapply(function(df_group) {
    if (nrow(df_group) > 2) {
      df_group %>%
        group_by(max_shom_bathy) %>%
        mutate(bathy_sub_id = cur_group_id()) %>%
        ungroup() %>%
        mutate(
          bathy_sub_id = dense_rank(bathy_sub_id),
          bathy_subgroup = paste0(overlap_date_group[1], "_b", bathy_sub_id)
        ) %>%
        select(-bathy_sub_id)
    } else {
      df_group$bathy_subgroup <- df_group$overlap_date_group
      df_group
    }
  }) %>%
  bind_rows() %>%
  
  # Step 3: Assign clean replicate ID 
  group_by(bathy_subgroup) %>%
  mutate(replicate = cur_group_id()) %>%
  select(-overlap_date_group, -bathy_subgroup) %>%
  ungroup() %>%
  
  # ---- Step 4: Manual Adjustments 
  # Force SPY232834 and SPY232833 into the same replicate group (assign lowest of their current replicate values)
  mutate(
    replicate = case_when(
      spygen_code %in% c("SPY232834", "SPY232833") ~ min(replicate[spygen_code %in% c("SPY232834", "SPY232833")]),
      TRUE ~ replicate
    )
  ) 


# Check the result ------------
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
    values = c("1" = "red", "2" = "#91cf60", "3" = "#d9ef8b", "4" = "#fee08b", "5" = "#fc8d59", "6" = "#d73027"),
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
