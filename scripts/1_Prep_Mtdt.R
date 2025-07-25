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
library(tidyr)


# Load functions
source("./utils/Fct_Data-Prep.R")



#------------- Load Metadata ---------------------
# Load the metadata file
mtdtfull <- readr::read_csv("./data/raw_data/eDNA/Med_metadonnees_ADNe_2018_2024_V4.csv")

#------------- Clean Metadata ---------------------
# Replicates errors ----
# Replicates are 2 or more samples taken at the same location (within a 500m radius), on the same 1/2 day (≤ 6h difference), at the same sampling depth (+/- 5m if done using the same method).
# This definition allows to group together samples for which the extraction of MARS3D environmental variables will be virtually the same. It doesn't necessarily correspond to the replicates as they were intended in the field.
# The replicates were pooled together by Amandine Avouac (based on various metadata) but some modifications need to be made to correspond to the above definition (cf mail sent on 22/07/2025 to Amandine and Laure). 

# Helper to standardize replicate naming
standardize_replicates <- function(codes) {
  paste(sort(codes), collapse = "/")
}

# Define replicate groups
replicate_groups <- list(
  # 2018
  c("SPY180630", "SPY180635", "SPY181822", "SPY181823"),
  c("SPY180634", "SPY180934", "SPY181811", "SPY181820"),
  c("SPY181810", "SPY181814", "SPY181829", "SPY181830"),
  
  # 2019
  c("SPY182155", "SPY182160", "SPY182543", "SPY182544"),
  c("SPY181975", "SPY182168", "SPY182540", "SPY182542"),
  c("SPY192280", "SPY192289", "SPY192305", "SPY192312"),
  c("SPY180765", "SPY181884", "SPY181894", "SPY181910",
    "SPY180785", "SPY181511", "SPY181890", "SPY181893"),
  
  # 2020
  c("SPY180622", "SPY201155", "SPY201156", "SPY201189"),
  c("SPY200460", "SPY200470", "SPY200455", "SPY200484"),
  c("SPY200456", "SPY200465", "SPY200468", "SPY200474"),
  c("SPY200462", "SPY201154", "SPY200472", "SPY200476"),
  c("SPY200457", "SPY200461", "SPY200489", "SPY200494"),
  c("SPY180629", "SPY190597", "SPY181145", "SPY192300"),
  c("SPY180625", "SPY181700", "SPY181701", "SPY181708"),
  c("SPY181824", "SPY181828", "SPY181699", "SPY192310"),
  c("SPY200463", "SPY200464", "SPY200458", "SPY200487", 
    "SPY200459", "SPY200488", "SPY200485", "SPY200493"),
  c("SPY201175", "SPY201191", "SPY202496", "SPY202498"),
  c("SPY201179", "SPY202492", "SPY201176", "SPY202499"),
  
  # 2021
  c("SPY212716", "SPY212756", "SPY212715", "SPY212717"),
  c("SPY212711", "SPY212759", "SPY212713", "SPY212758"),
  c("SPY212710", "SPY212712", "SPY212719"),
  
  # 2022 - split, so handled separately later
  
  # 2023
  c("SPY232138", "SPY233971", "SPY232132", "SPY232197"),
  c("SPY232130", "SPY232272", "SPY232135", "SPY232154"),
  c("SPY234951", "SPY234961", "SPY234952", "SPY234957"),
  c("SPY233998", "SPY234954", "SPY234956", "SPY234960", "SPY234955", "SPY234959"),
  c("SPY232255", "SPY233721")
)

# Update replicates
mtdtfull <- mtdtfull %>%
  rowwise() %>%
  mutate(replicates = {
    new_repl <- NULL
    for (grp in replicate_groups) {
      if (spygen_code %in% grp) {
        new_repl <- standardize_replicates(grp)
        break
      }
    }
    ifelse(!is.null(new_repl), new_repl, replicates)
  }) %>%
  ungroup()



# 2022: split SPY220281/SPY220285 into separate replicates
mtdtfull <- mtdtfull %>%
  mutate(replicates = case_when(
    spygen_code == "SPY220281" ~ "SPY220281",
    spygen_code == "SPY220285" ~ "SPY220285",
    TRUE ~ replicates
  ))

# Optional: ensure consistency in naming
mtdtfull <- mtdtfull %>%
  mutate(replicates = ifelse(is.na(replicates) | replicates == "", "no", replicates))

rm(replicate_groups, standardize_replicates)

# Other errors ----

# Remove "...1" column
mtdtfull <- mtdtfull %>%
  select(-contains("...1"))



# sampling_depth 
if (any(mtdtfull$spygen_code == "SPY232658")) {
  mtdtfull <- mtdtfull %>%
    mutate(depth_sampling = ifelse(spygen_code == "SPY232658", 80, depth_sampling)) 
}

if (any(mtdtfull$spygen_code == "SPY232652")) {
  mtdtfull <- mtdtfull %>%
    mutate(depth_sampling = ifelse(spygen_code == "SPY232652", 20, depth_sampling)) 
}



# coordinates
if (any(mtdtfull$spygen_code == "SPY2401010")) {
  mtdtfull <- mtdtfull %>%
    mutate(latitude_end_DD = ifelse(spygen_code == "SPY2401010", 42.750149999999998, latitude_end_DD))
}

if (any(mtdtfull$spygen_code == "SPY2400997")) {
  mtdtfull <- mtdtfull %>%
    mutate(longitude_end_DD = ifelse(spygen_code == "SPY2400997", 9.08493333333333, longitude_end_DD),
           latitude_end_DD = ifelse(spygen_code == "SPY2400997", 42.750149999999998, latitude_end_DD))
}

if (any(mtdtfull$spygen_code == "SPY2401011")) {
  mtdtfull <- mtdtfull %>%
    mutate(latitude_end_DD = ifelse(spygen_code == "SPY2401011", 42.753000000000000, latitude_end_DD),
           longitude_end_DD = ifelse(spygen_code == "SPY2401011", 9.176416666666670, longitude_end_DD))
}



# component
mtdtfull <- mtdtfull %>%
  mutate(component = ifelse(spygen_code == "SPY212587", "offshore", component))

mtdtfull <- mtdtfull %>%
  mutate(component = ifelse(spygen_code == "SPY211262", "offshore", component))

mtdtfull <- mtdtfull %>%
  mutate(component = ifelse(spygen_code == "SPY211311", "offshore", component))







#------------- Load IPOCOM ---------------------
ipocom <- readr::read_csv("./data/raw_data/eDNA/AB_eREF_IPOCOM TRANSECT 2024_metadatas template.csv")




#------------- Clean IPOCOM ---------------------
# Coordinates errors ----

if (any(ipocom$subsite_andromede == "IPOCOM_12")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_12", 43.53943333, latitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_13")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_13", 43.5057833, latitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_22")) {
  ipocom <- ipocom %>%
    mutate(longitude_end_DD = ifelse(subsite_andromede == "IPOCOM_22", 6.854, longitude_end_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_34")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_34", 43.246, latitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_40")) {
  ipocom <- ipocom %>%
    mutate(latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_40", 43.148, latitude_end_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_58")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_58", 43.068, latitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_62")) {
  ipocom <- ipocom %>%
    mutate(
      latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_62", 43.04214, latitude_start_DD),
      longitude_start_DD = ifelse(subsite_andromede == "IPOCOM_62", 5.46331, longitude_start_DD),
      latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_62", 43.05471, latitude_end_DD),
      longitude_end_DD = ifelse(subsite_andromede == "IPOCOM_62", 5.46201, longitude_end_DD)
    )
}

if (any(ipocom$subsite_andromede == "IPOCOM_76")) {
  ipocom <- ipocom %>%
    mutate(longitude_end_DD = ifelse(subsite_andromede == "IPOCOM_76", 5.34973333333333, longitude_end_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_82")) {
  ipocom <- ipocom %>%
    mutate(latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_82", 43.317, latitude_end_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_135")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_135", 42.4785833333333, latitude_start_DD))
}


# Depth_sampling missing ----
# Mail Adele Barroil 25/07/25 : "Pour ces deux transects on a eu un problème avec la sonde CTD donc pas de données. Mais pour SPY2401554/SPY2401555 ca devait être entre 20 et 25 m et pour SPY2401674/SPY2401675 vers 10m"

if (any(ipocom$replicates == "SPY2401554/SPY2401555")) {
  ipocom <- ipocom %>%
    mutate(depth_sampling = ifelse(replicates == "SPY2401554/SPY2401555", 22.5, depth_sampling))
}

if (any(ipocom$replicates == "SPY2401674/SPY2401675")) {
  ipocom <- ipocom %>%
    mutate(depth_sampling = ifelse(replicates == "SPY2401674/SPY2401675", 10, depth_sampling))
}



# Format ipocom to match mtdtfull ----
ipocom <- ipocom %>%
  rename(BiodivMed2023 = BiodivMed)

ipocom$date <- as.Date(ipocom$date, format = "%d/%m/%Y")

# Convert hms/difftime to numeric seconds
ipocom$duration <- as.numeric(ipocom$duration, units = "secs")

# Convert to minutes if mtdtfull uses minutes
ipocom$duration <- ipocom$duration / 60







#------------- Merge IPOCOM to Metadata -----------------

# Combine the two datasets
mtdtcomb <- bind_rows(mtdtfull, ipocom)

# Combine the two datasets and keep only the columns present in mtdtfull filling with NA if necessary

# Clean up
rm(mtdtfull, ipocom)








#------------- Subset n°1 : France + marine + coastal --------------------
# Explanation:
# - Selection of points sampled in France 
# - Selection of marine points: outside harbours, lagoons, ports, rivers, estuaries 
# - Filtering out the seamount, canyon and open_ocean samples 
# - BUT: keep any sample in a replicate group if one of them passes the filters
# - AND: if replicates == "no", it must pass the filters individually

# STEP 1: Identify valid replicate groups (excluding "no")
valid_replicates <- mtdtcomb %>%
  filter(replicates != "no") %>%
  filter(is.na(country) | country == "France") %>%
  filter(is.na(component) | !(component %in% c(
    "harbour", "lagoon", "port", 
    "freshwater_river", "estuary", 
    "open_ocean", "seamount", "canyon"
  ))) %>%
  pull(replicates) %>%
  unique()

# STEP 2: Keep:
# - all samples in valid replicate groups
# - singletons (replicates == "no") only if they meet the filters individually
mtdt_1 <- mtdtcomb %>%
  filter(
    (replicates %in% valid_replicates) |
      (
        replicates == "no" &
          (is.na(country) | country == "France") &
          (is.na(component) | !(component %in% c(
            "harbour", "lagoon", "port", 
            "freshwater_river", "estuary", 
            "open_ocean", "seamount", "canyon"
          )))
      )
  )

rm(valid_replicates)

#------------- Extract bathymetry ---------------------
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

#------------- Subset n°2 : remove offshore > 50m depth --------------------
# Explanation : we want to remove samples far away from the coast, yet some offshore samples are close enough to be kept (since offshore was assigned to samples > 6km away from the coast). Thus we remove samples that are BOTH offshore AND > 50m depth, this ensures keeping coastal samples relatively close to shore.

# STEP 1: Identify valid replicate groups (excluding "no")
valid_replicates_component_bathy <- mtdt_1 %>%
  filter(replicates != "no") %>%
  filter(
    is.na(component) | component != "offshore" |
      (!is.na(combined_bathy) & combined_bathy <= 50)
  ) %>%
  pull(replicates) %>%
  unique()

# STEP 2: Keep:
# - all rows in valid replicate groups
# - singletons (replicates == "no") only if they are not deep offshore
mtdt_2 <- mtdt_1 %>%
  filter(
    (replicates %in% valid_replicates_component_bathy) |
      (
        replicates == "no" &
          (
            is.na(component) | component != "offshore" |
              (!is.na(combined_bathy) & combined_bathy <= 50)
          )
      )
  )

# Clean up
rm(valid_replicates_component_bathy)



#------------- Buffer samples ---------------------
mtdt_2 <- buffer_transect(
  df = mtdt_2,
  start_lon_col = "longitude_start_DD",
  start_lat_col = "latitude_start_DD",
  end_lon_col = "longitude_end_DD",
  end_lat_col = "latitude_end_DD",
  buffer_dist = 500
)

# ggplot(mtdt_2) +
#   geom_sf(aes(geometry = wkt_geometry), fill = "lightblue", color = "darkblue", alpha = 0.1) +
#   theme_minimal() +
#   ggtitle("500m samples buffer 2028-2022")

# Save as gpkg
mtdt_2$time_start <- as.character(mtdt_2$time_start)
sf::write_sf(mtdt_2, "./data/processed_data/eDNA/mtdt_2.gpkg", delete_dsn = TRUE)

#------------- Buffer replicates ------------

# Group by replicate and compute concave hulls
replicates_buffer <- mtdt_2 |>
  filter(replicates != "no") |>
  group_split(replicates) |>
  
  # Apply concave hull algorithm to each replicate group
  purrr::map_dfr(function(mtdt_2) {
    hull <- concaveman::concaveman(mtdt_2, concavity = 10, length_threshold = 10)   # length_threshold : length of the shortest edge in the hull
    tibble(
      replicate = unique(mtdt_2$replicates),
      geometry = st_geometry(hull)
    )
  }) |>
  st_as_sf(crs = 4326)


# Drop Z and M dimensions (convert ZM to XY)
replicates_buffer <- st_zm(replicates_buffer, drop = TRUE, what = "ZM")

# Combine replicates_buffer for replicates != "no" and wkt_geometry for replicates == "no" in col named "replicates_geometry"
mtdt_3 <- mtdt_2 |>
  rename(sample_geometry = wkt_geometry) |>
  rowwise() |>
  mutate(replicates_geometry = if (replicates != "no") {
    matched_index <- match(replicates, replicates_buffer$replicate)
    if (!is.na(matched_index)) {
      st_as_text(replicates_buffer$geometry[[matched_index]])
    } else {
      NA_character_
    }
  } else {
    st_as_text(sample_geometry)
  }) |>
  ungroup()

# Convert mtdt_3 to non spatial data frame
mtdt_3 <- st_drop_geometry(mtdt_3)

# Convert mtdt_3 to spatial object with replicates_geometry
mtdt_3 <- st_as_sf(mtdt_3, wkt = "replicates_geometry", crs = 4326)

# Plot
plot(mtdt_3[1])

# Clean 
rm(replicates_buffer)












#------------- Export ------------
# Save as gpkg
mtdt_3$time_start <- as.character(mtdt_3$time_start)
sf::write_sf(mtdt_3, "./data/processed_data/eDNA/mtdt_3.gpkg", delete_dsn = TRUE)

#------------- Group Mtdt by replicates -----------------
# Explanation : make a simplified dataset grouped by replicates for NCDF extraction 


# --- Simplified mtdt for MARS3D extraction ----
# Prepare dataset ----
# Keep "spygen_code", "replicates", "date", "time_start", country", "region", "site", "component", "replicates_geometry", "depth_sampling"
mtdt_4 <- mtdt_3 %>%                
  select(spygen_code, replicates, date, time_start, country, region, replicates_geometry, depth_sampling) 

# Replace time_start = NULL by 7:00
mtdt_4$time_start[is.na(mtdt_4$time_start)] <- "07:00:00"

# Check if any mtdt_4$time_start is NA
any(is.na(mtdt_4$time_start))

# Replace replicates = "no" by "spygen_code"
mtdt_4$replicates[mtdt_4$replicates == "no"] <- mtdt_4$spygen_code[mtdt_4$replicates == "no"]

# Remove spygen_code
mtdt_4 <- mtdt_4 %>%                
  select(-spygen_code)




# Group by replicates ----
# When diff depth_sampling --> prints the rows and make the mean across depth_sampling 
mtdt_4 <- mtdt_4 %>%
  select(replicates, date, time_start, country, region, replicates_geometry, depth_sampling) %>%
  group_by(replicates) %>%
  summarise(
    date = first(date),
    time_start = min(time_start, na.rm = TRUE),
    
    # Use first value if all non-NA values are identical; otherwise compute mean
    depth_sampling = if (n_distinct(na.omit(depth_sampling)) == 1) {
      first(na.omit(depth_sampling))
    } else {
      print(cur_data_all())  # This prints the full data for the current group
      mean(depth_sampling, na.rm = TRUE)
    },
    
    country = first(country),
    region = first(region),
    replicates_geometry = first(replicates_geometry),
    .groups = "drop"
  )





# Export ----
# Date and time for Gaétan (aggregation of MARS3D data) and Martin (extraction Canyon) 
# Sent on 24/07/2025
mtdt_4 %>%
  select(replicates, date, time_start) %>%
  write_csv("./data/processed_data/eDNA/mtdt_4_date_time.csv")










# Geom date, time, depth_sampling for Paule (MARS3D extraction) 
st_write(mtdt_4, "./data/processed_data/eDNA/mtdt_4.gpkg", delete_dsn = TRUE)

























# --- Complete mtdt extraction ----
# Prepare dataset ----
mtdt_5 <- mtdt_3 

# Replace time_start = NULL by 7:00
mtdt_5$time_start[is.na(mtdt_5$time_start)] <- "07:00:00"

# Check if any mtdt_4$time_start is NA
any(is.na(mtdt_5$time_start))

# Replace replicates = "no" by "spygen_code"
mtdt_5$replicates[mtdt_5$replicates == "no"] <- mtdt_5$spygen_code[mtdt_5$replicates == "no"]

# Remove spygen_code
mtdt_5 <- mtdt_5 %>%                
  select(-spygen_code)




# Group by replicates ----
mtdt_5 <- mtdt_5 %>%
  select(-c("latitude_start_DD", "longitude_start_DD", 
                    "latitude_end_DD", "longitude_end_DD", 
                    "pool", "subsite_andromede", "max_shom_bathy", "combined_bathy", "mpa_dist"))


mtdt_5 <- mtdt_5 %>%
  group_by(replicates) %>%
  summarise(
    date = first(na.omit(date)),
    time_start = min(time_start, na.rm = TRUE),
    
    # Use first value if all non-NA values are identical; otherwise compute mean
    depth_sampling = if (n_distinct(na.omit(depth_sampling)) == 1) {
      first(na.omit(depth_sampling))
    } else {
      print(cur_data_all())  # This prints the full data for the current group
      mean(depth_sampling, na.rm = TRUE)
    },
    
    # Use first value if all non-NA values are identical; otherwise take the first non-NA value
    depth_seafloor = if(n_distinct(na.omit(depth_seafloor)) == 1) {
      first(na.omit(depth_seafloor))
    } else {
      print(cur_data_all())  # This prints the full data for the current group
      first(na.omit(depth_seafloor))
    },
    
    lockdown = first(na.omit(lockdown)),
    BiodivMed2023 = first(na.omit(BiodivMed2023)),
    method = first(na.omit(method)),
    country = first(na.omit(country)),
    region = first(na.omit(region)),
    site = first(na.omit(site)),
    subsite = first(na.omit(subsite)),
    component = first(na.omit(component)),
    habitat = first(na.omit(habitat)),
    protection = first(na.omit(protection)),
    mpa_name = first(na.omit(mpa_name)),
    replicates_geometry = first(na.omit(replicates_geometry)),
    Tele01 = first(na.omit(Tele01)),
    Pleo = first(na.omit(Pleo)),  
    Mamm01 = first(na.omit(Mamm01)),
    Vert01 = first(na.omit(Vert01)),
    X16s_Metazoa = first(na.omit(X16s_Metazoa)),
    Bact02 = first(na.omit(Bact02)), 
    Euka02 = first(na.omit(Euka02)),
    
  # if comments differ combine them by copy pasting them with their associated replicates id :
    comments = paste(unique(na.omit(comments)), collapse = "; "),  # Combine unique comments with a semicolon
    .groups = "drop"
  )






