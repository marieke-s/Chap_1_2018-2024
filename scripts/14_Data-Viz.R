#------------- Description ----
# The aim of this script is to explore all data (Mtdt, species, predictos, div_indices both raw and transformed)





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
library(patchwork)
library(lubridate)
library(terra)
library(stringr)
library(pMEM)
library(ggplot2)




#------------- Load functions 
source("./utils/Fct_Data-Prep.R")

#------------------------------------------------- MAP GRID CELLS x sampling x time ------------------------------
# Load data -----
grid <- st_read("./data/processed_data/prediction_extent/grid_v1.1.gpkg")
buff <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

str(buff)
str(grid)



# Prep data -----
buff <- buff %>%
  mutate(year = year(date),
         month = month(date))
str(buff)


# Intersect grid and buff by month_year -----
sf_use_s2(FALSE)  # turn off s2 globally for sf

## 1. Create month_year in buff ("Janvier 2018", "Février 2018", etc.)

mois_fr <- c("Janvier", "Février", "Mars", "Avril", "Mai", "Juin",
             "Juillet", "Août", "Septembre", "Octobre", "Novembre", "Décembre")

buff$month_year <- paste(mois_fr[buff$month], buff$year)

# (Optionally keep them ordered chronologically)
buff$month_year <- factor(buff$month_year,
                          levels = unique(buff[order(buff$year, buff$month), "month_year", drop = TRUE]))

## 2. Make sure buff and grid share the same CRS

if (st_crs(buff) != st_crs(grid)) {
  grid <- st_transform(grid, st_crs(buff))
}

## 3. For each unique month_year, intersect and create sampling_{month_year} in grid

unique_month_year <- levels(buff$month_year)  # or unique(buff$month_year) if not factor

for (my in unique_month_year) {
  # subset buff to that month_year
  buff_my <- buff[buff$month_year == my, ]
  if (nrow(buff_my) == 0) next  # just in case
  
  # build a safe column name: sampling_Janvier_2018, sampling_Fevrier_2018, etc.
  col_name <- paste0("sampling_",
                     gsub("[^A-Za-z0-9]", "_", my))  # replace spaces/accents/punctuation by "_"
  
  # initialize column as FALSE
  grid[[col_name]] <- FALSE
  
  # find intersections grid–buff_my
  # st_intersects() returns a list of integer vectors (indices of buff_my that intersect each grid cell)
  inter <- st_intersects(grid, buff_my)
  
  # cells with at least one intersecting buffer
  hits <- lengths(inter) > 0
  
  # set TRUE where there is at least one intersection
  grid[[col_name]][hits] <- TRUE
}


plot(grid["sampling_Avril_2023"])

## 4. Export the resulting grid (choose whatever format you need)

# Example: GeoPackage
st_write(grid, "./data/processed_data/Data-viz/intersection-grid_v1-mtdt_7.gpkg", layer = "grid_sampling", delete_layer = TRUE)






#----- Make plots ----
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)

grid <- st_read("./data/processed_data/Data-viz/intersection-grid_v1-mtdt_7.gpkg", layer = "grid_sampling")

# 1. Identify the sampling columns
sampling_cols <- grep("^sampling_", names(grid), value = TRUE)

# 2. Separate geometry so we can pivot the attributes
grid_geom <- grid[, "geom"]           # keeps sf geometry + row order
grid_attr <- st_drop_geometry(grid)   # attributes only

# 3. Long format: one row per cell per month_year column
grid_long <- grid_attr %>%
  dplyr::select(id, all_of(sampling_cols)) %>%
  pivot_longer(
    cols      = all_of(sampling_cols),
    names_to  = "month_year_col",
    values_to = "sampling"
  )

# 4. Reconstruct a readable month_year label from the column names
#    e.g. "sampling_Janvier_2018" -> "Janvier 2018"
grid_long <- grid_long %>%
  mutate(
    month_year = gsub("^sampling_", "", month_year_col),
    month_year = gsub("_", " ", month_year)
  )

# Order month_year chronologically 
grid_long$month_year <- factor(grid_long$month_year,
                               levels = unique_month_year)


# 5. Join back geometry
grid_long_sf <- st_as_sf(
  left_join(grid_geom %>% mutate(id = grid_attr$id),
            grid_long,
            by = "id")
)

# 6. Build the plot: sampled cells in yellow, others in light grey
p <- ggplot(grid_long_sf) +
  geom_sf(aes(fill = sampling), color = NA) +
  scale_fill_manual(
    values = c(`FALSE` = "grey80", `TRUE` = "yellow"),
    guide = "none"
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    title    = "Sites échantillonnés",
    subtitle = "{closest_state}",
    x = NULL, y = NULL
  ) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )

# 7. Animate over month_year
anim <- p +
  transition_states(
    states = month_year,
    state_length = 0.7,
    transition_length = 0.3
  ) +
  ease_aes("linear")

# 8. Render and save GIF
n_states <- length(levels(grid_long_sf$month_year))

anim_gif <- animate(
  anim,
  nframes = length(levels(grid_long_sf$month_year)),  # 1 frame per month
  fps = 2,
  width   = 800,
  height  = 600
)


anim_save("./figures/sampling_through_time.gif", animation = anim_gif)













#------------------------------------------------- MAP ALL MED eDNA points FRANCE ------------------------------
# Load the metadata file
mtdtfull <- readr::read_csv("./data/raw_data/Mtdt/Med_metadonnees_ADNe_2018_2024_V4.csv")


# Load the metadata file
mtdtfull <- readr::read_csv("./data/raw_data/Mtdt/Med_metadonnees_ADNe_2018_2024_V4.csv")

#------------- Clean Metadata ---------------------
# Replicates errors ---
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
  c("SPY212710", "SPY212712"),
  
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

# Remove SPY212719 --> it could be associated with SPY212710/SPY212712 but one is the MPA whether the 2 other are not. It doesn't have any replicates but is well to close to SPY212710/SPY212712 to be kept alone. 

mtdtfull <- mtdtfull %>%
  filter(!(spygen_code == "SPY212719"))

# Optional: ensure consistency in naming
mtdtfull <- mtdtfull %>%
  mutate(replicates = ifelse(is.na(replicates) | replicates == "", "no", replicates))

rm(replicate_groups, standardize_replicates)

# Other errors ---

# Remove "...1" column
mtdtfull <- mtdtfull %>%
  dplyr::select(-contains("...1"))



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







# estimated_volume
# Explanation : pump filter seawater at 1L/min, thus estimated_volume can be retrieved using duration 

mtdtfull <- mtdtfull %>%
  mutate(estimated_volume = ifelse(is.na(estimated_volume) & !is.na(duration), duration, estimated_volume)) 

mtdtfull <- mtdtfull %>%
  mutate(estimated_volume = ifelse(replicates == "SPY211147_SPY211148", 16, estimated_volume))







# country 
mtdtfull <- mtdtfull %>%
  mutate(country = ifelse(replicates == "SPY180926/SPY180931", "Spain", country),
         country = ifelse(replicates == "SPY181817/SPY181827/SPY182545/SPY182546", "Spain", country),
         country = ifelse(replicates == "SPY180928/SPY180929", "Spain", country),
         country = ifelse(replicates == "SPY182155/SPY182160/SPY182543/SPY182544", "Spain", country)
  )







#------------- Load IPOCOM 
ipocom <- readr::read_csv("./data/raw_data/Mtdt/AB_eREF_IPOCOM TRANSECT 2024_metadatas template.csv")




#------------- Clean IPOCOM 
# Coordinates errors ---

if (any(ipocom$subsite_andromede == "IPOCOM_12")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_12", 43.53943333, latitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_13")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_13", 43.54591, latitude_start_DD), 
           longitude_start_DD = ifelse(subsite_andromede == "IPOCOM_13", 7.08869, longitude_start_DD))
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
      latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_62", 43.0702333333333, latitude_start_DD),
      longitude_start_DD = ifelse(subsite_andromede == "IPOCOM_62", 5.77218333333333, longitude_start_DD),
      latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_62", 43.0908, latitude_end_DD),
      longitude_end_DD = ifelse(subsite_andromede == "IPOCOM_62", 5.7704, longitude_end_DD)
    )
}

if (any(ipocom$subsite_andromede == "IPOCOM_76")) {
  ipocom <- ipocom %>%
    mutate(longitude_end_DD = ifelse(subsite_andromede == "IPOCOM_76", 5.2837, longitude_end_DD), 
           latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_76", 43.2586, latitude_end_DD)
           
           
    )
}

if (any(ipocom$subsite_andromede == "IPOCOM_82")) {
  ipocom <- ipocom %>%
    mutate(latitude_end_DD = ifelse(subsite_andromede == "IPOCOM_82", 43.317, latitude_end_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_111")) {
  ipocom <- ipocom %>%
    mutate(
      latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_111", 43.3705166666667, latitude_start_DD),
      longitude_start_DD = ifelse(subsite_andromede == "IPOCOM_111", 3.66183333333333, longitude_start_DD))
}

if (any(ipocom$subsite_andromede == "IPOCOM_135")) {
  ipocom <- ipocom %>%
    mutate(latitude_start_DD = ifelse(subsite_andromede == "IPOCOM_135", 42.4785833333333, latitude_start_DD))
}


# Depth_sampling missing ---
# Mail Adele Barroil 25/07/25 : "Pour ces deux transects on a eu un problème avec la sonde CTD donc pas de données. Mais pour SPY2401554/SPY2401555 ca devait être entre 20 et 25 m et pour SPY2401674/SPY2401675 vers 10m"

if (any(ipocom$replicates == "SPY2401554/SPY2401555")) {
  ipocom <- ipocom %>%
    mutate(depth_sampling = ifelse(replicates == "SPY2401554/SPY2401555", 22.5, depth_sampling))
}

if (any(ipocom$replicates == "SPY2401674/SPY2401675")) {
  ipocom <- ipocom %>%
    mutate(depth_sampling = ifelse(replicates == "SPY2401674/SPY2401675", 10, depth_sampling))
}



# Format ipocom to match mtdtfull ---
ipocom <- ipocom %>%
  rename(BiodivMed2023 = BiodivMed)

ipocom$date <- as.Date(ipocom$date, format = "%d/%m/%Y")

# Convert hms/difftime to numeric seconds
ipocom$duration <- as.numeric(ipocom$duration, units = "secs")

# Convert to minutes if mtdtfull uses minutes
ipocom$duration <- ipocom$duration / 60







#------------- Merge IPOCOM to Metadata 

# Combine the two datasets
mtdtcomb <- bind_rows(mtdtfull, ipocom)

# Combine the two datasets and keep only the columns present in mtdtfull filling with NA if necessary

# Clean up
rm(mtdtfull, ipocom)








#------------- Subset France  --------------------

# STEP 1: Identify valid replicate groups (excluding "no")
valid_replicates <- mtdtcomb %>%
  filter(replicates != "no") %>%
  filter(is.na(country) | country == "France") %>%
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
          (is.na(country) | country == "France") )
      )



rm(valid_replicates)


#------------- Group by replicates  --------------------

# Replace replicates = "no" by "spygen_code"
mtdt_1$replicates[mtdt_1$replicates == "no"] <- mtdt_1$spygen_code[mtdt_1$replicates == "no"]

# Remove spygen_code
mtdt_1 <- mtdt_1 %>%                
  dplyr::select(-spygen_code)



mtdt_full_replicates <- mtdt_1 %>%
  group_by(replicates) %>%
  summarise(
    date = first(na.omit(date)),
    time_start = min(time_start, na.rm = TRUE),
    
    # Use first value if all non-NA values are identical; otherwise compute mean
    depth_sampling = if (n_distinct(na.omit(depth_sampling)) == 1) {
      first(na.omit(depth_sampling))
    } else {
      #print(cur_data_all())  # This prints the full data for the current group
      mean(depth_sampling, na.rm = TRUE)
    },
    
    # Use first value if all non-NA values are identical; otherwise take the first non-NA value
    depth_seafloor = if(n_distinct(na.omit(depth_seafloor)) == 1) {
      first(na.omit(depth_seafloor))
    } else {
      #print(cur_data_all())  # This prints the full data for the current group
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
    project = first(na.omit(project)),
    Tele01 = first(na.omit(Tele01)),
    Pleo = first(na.omit(Pleo)),  
    Mamm01 = first(na.omit(Mamm01)),
    Vert01 = first(na.omit(Vert01)),
    X16s_Metazoa = first(na.omit(X16s_Metazoa)),
    Bact02 = first(na.omit(Bact02)), 
    Euka02 = first(na.omit(Euka02)),
    estimated_volume_total = sum(na.omit(estimated_volume)),
    duration_total = sum(na.omit(duration)),
    latitude_start_DD = first(na.omit(latitude_start_DD)),
    longitude_start_DD = first(na.omit(longitude_start_DD)),
    latitude_end_DD = first(na.omit(latitude_end_DD)),
    longitude_end_DD = first(na.omit(longitude_end_DD)),
    
    # if comments differ combine them by copy pasting them with their associated replicates id :
    comments = paste(unique(na.omit(comments)), collapse = "; "),  # Combine unique comments with a semicolon
    .groups = "drop"
  )




# --- subset france again ----
mtdt_full_replicates <- mtdt_full_replicates %>%
  filter(is.na(country) | country == "France")

mtdt_full_replicates <- mtdt_full_replicates %>%
  filter(latitude_start_DD >= 41,
         latitude_start_DD <= 44)

# --- remove blank ----
mtdt_full_replicates <- mtdt_full_replicates %>%
  filter(!component == "blank")


# --- remove glider ----
mtdt_full_replicates <- mtdt_full_replicates %>%
  filter(!method == "glider")

# --- correct coordinates ----
# set end coordinates of replicates "SPY2401047/SPY2401048" to NA
mtdt_full_replicates <- mtdt_full_replicates %>%
  mutate(
    latitude_end_DD = ifelse(replicates == "SPY2401047/SPY2401048", NA, latitude_end_DD),
    longitude_end_DD = ifelse(replicates == "SPY2401047/SPY2401048", NA, longitude_end_DD)
  )


#--------------- Map data prep ---------------------
dt <- mtdt_full_replicates

# Make dt an sf object with POINT geometry from latitude and longitude start
dt <- st_as_sf(
  dt,
  coords = c("longitude_start_DD", "latitude_start_DD"),
  crs = 4326,
  remove = FALSE
)


# Make transect geometry if latitude_end_DD and longitude_end_DD are available
# dt <- dt %>%
#   mutate(
#     longitude_start_DD = as.numeric(longitude_start_DD),
#     latitude_start_DD  = as.numeric(latitude_start_DD),
#     longitude_end_DD   = as.numeric(longitude_end_DD),
#     latitude_end_DD    = as.numeric(latitude_end_DD)
#   )
# 
# geom_list <- purrr::pmap(
#   list(dt$longitude_start_DD, dt$latitude_start_DD,
#        dt$longitude_end_DD,   dt$latitude_end_DD),
#   function(lon_s, lat_s, lon_e, lat_e) {
#     
#     has_end <- !is.na(lon_e) && !is.na(lat_e)
#     
#     if (has_end) {
#       st_linestring(matrix(c(lon_s, lat_s,
#                              lon_e, lat_e), 
#                            ncol = 2, byrow = TRUE))
#     } else {
#       st_point(c(lon_s, lat_s))
#     }
#   }
# )
# 
# dt_sf <- st_sf(
#   dt,
#   geometry = st_sfc(geom_list, crs = 4326)
# )


# plot(dt_sf["geometry"])
# dt <- dt_sf
# Make year col
dt$year <- year(dt$date)


# renames port
dt <- dt %>%
  mutate(component = ifelse(component == "harbour", "port", component))
table(dt$component)

# translate all components in french
dt <- dt %>%
  mutate(component = case_when(
    component == "coastal" ~ "Côte",
    component == "offshore" ~ "Large",
    component == "port" ~ "Port",
    component == "estuary" ~ "Estuaire",
    component == "lagoon" ~ "Lagune",
    component == "freshwater_river" ~ "Rivière",
    component == "open_ocean" ~ "Haute mer",
    component == "seamount" ~ "Mont sous-marin",
    component == "canyon" ~ "Canyon",
    TRUE ~ component
  ))

# translate all methods in french
dt <- dt %>%
  mutate(method = case_when(
    method == "surface_transect" ~ "Transect de surface",
    method == "seabed_transect" ~ "Transect de fond",
    method == "dive_transect" ~ "Transect en plongée",
    method == "motionless_descended_from_surface" ~ "Pompe descendue depuis la surface",
    method == "dive_motionless" ~ "Point fixe en plongée",
    method == "AUV" ~ "Robot autonome",
    TRUE ~ method
  ))






#----- Background ---------------------
# Load regions
regions <- st_read("./data/processed_data/Data-viz/med_regions.geojson")

# Transform CRS to 3857
dt <- st_transform(dt, 3857)
regions <- st_transform(regions, 3857)

# Crop regions to bounding box of dt
bb <- st_bbox(dt)
world_crop <- st_crop(regions, bb)








#----- OUTPUT PATH -----
output_path <- "./figures/aquarium_2025/"
dir.create(output_path, showWarnings = FALSE)

#--------------- Map year black --------------------
years <- 2017:2024
for (y in years) {
  dt_cum <- dt %>% filter(year <= y)
  
  # dynamic title text
  if (y == 2017) {
    title_txt <- "2017"
  } else if (y == 2018) {
    title_txt <- "2018"
  } else {
    title_txt <- paste0("2018-", y)
  }
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = dt_cum, color = "white", size = 1.5) +
    ggtitle(title_txt) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"])) +
    theme_dark() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text = element_text(color = "white", size = 10),
      axis.title = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title = element_text(color = "white", size = 26, face = "bold"),
      panel.border = element_blank()
    )
  
  print(p)
  
  # ---- SAVE FIGURE ----
  ggsave(
    filename = file.path(output_path, paste0("year_", title_txt, ".png")),
    plot = p, width = 10, height = 8, dpi = 300, bg = "black"
  )
}



#--------------- Map methods (black) ------------------------
methods <- sort(na.omit(unique(dt$method)))

for (i in seq_along(methods)) {
  m    <- methods[i]
  dt_m <- dt %>% filter(method == m)
  
  method_color <- "#F333FF"
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = dt,    color = "grey40", size = 1.3) +
    geom_sf(data = dt_m,  color = method_color, size = 1.8) +
    ggtitle(m) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"])) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = method_color, size = 26, face = "bold"),
      panel.border     = element_blank()
    )
  
  print(p)
  
  # ---- SAVE FIGURE ----
  ggsave(
    filename = file.path(output_path, paste0("method_", m, ".png")),
    plot = p, width = 10, height = 8, dpi = 300, bg = "black"
  )
}


#--------------- Map components (black) ------------------------
components <- sort(na.omit(unique(dt$component)))

for (i in seq_along(components)) {
  m    <- components[i]
  dt_m <- dt %>% filter(component == m)
  
  comp_color <- "gold"
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = dt,    color = "grey40", size = 1.3) +
    geom_sf(data = dt_m,  color = comp_color, size = 1.8) +
    ggtitle(m) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"])) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = comp_color, size = 26, face = "bold"),
      panel.border     = element_blank()
    )
  
  print(p)
  
  # ---- SAVE FIGURE ----
  ggsave(
    filename = file.path(output_path, paste0("component_", m, ".png")),
    plot = p, width = 10, height = 8, dpi = 300, bg = "black"
  )
}







#------------------------------------------------- MAP COASTAL POINTS ------------------------------


#--------------- Data & Background ---------------------

mtdt_7 <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")
# predictors_sel_v.1.4 ----
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.4.gpkg")
pred <- pred %>% st_drop_geometry()

# Make dt an sf object with POINT geometry from latitude and longitude start
dt2 <- pred %>% st_as_sf(coords = c("x", "y"),
    crs = 4326,
    remove = FALSE
  )

dt2 <- st_transform(dt2, 3857)

# Load regions ----
regions <- st_read("./data/processed_data/Data-viz/med_regions.geojson")

# Transform CRS to 3857
regions <- st_transform(regions, 3857)

# Crop regions to bounding box of dt
bb <- st_bbox(dt)
world_crop <- st_crop(regions, bb)



# Map with same style as before ----

p2 <- ggplot() +
  geom_sf(data = world_crop, fill = "black", color = "grey70") +
  geom_sf(data = dt2, color = "white", size = 1.5) +
  ggtitle("Sites à 0-50m de profondeur") +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"])) +
  theme_dark() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    axis.text = element_text(color = "white", size = 10),
    axis.title = element_text(color = "white"),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_line(color = "grey25"),
    plot.title = element_text(color = "white", size = 26, face = "bold"),
    panel.border = element_blank()
  )

p2

# Save p2
ggsave(
  filename = file.path(output_path, paste0("coastal_points_2018-2024.png")),
  plot = p2, width = 10, height = 8, dpi = 300, bg = "black"
)











#------------------------------------------------ Extract 50m bathy -----------------------
bathy_corse <- terra::rast("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")
plot(bathy)

# Transform to 3857
bathy <- project(bathy, "EPSG:3857")
plot(bathy)

# Crop on dt bbox
bathy_crop <- terra::crop(bathy, c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
plot(bathy_crop)

# Extract 50m contour line
bathy_50m <- terra::as.contour(bathy_crop, levels = -50)
plot(bathy_50m)

# Convert to sf
bathy_50m_sf <- st_as_sf(bathy_50m)

# Convert multilines to single lines
bathy_50m_sf <- bathy_50m_sf %>%
  st_cast("LINESTRING")

# Export bathy_50m_sf
st_write(bathy_50m_sf, "./data/processed_data/Data-viz/bathy_50m_sf.gpkg", delete_dsn = TRUE)


# Load cleaned bathy (qgis manual cleaning)
b <- st_read("./data/processed_data/Data-viz/bathy_50m_cleaned.gpkg")
plot(b)

#------------------------------------------------ MAP ALL POINTS + 50M bathy --------------------



# Cumulative data up to 2024
dt_cum <- dt %>% filter(year <= 2024)
title_txt <- "2018-2024"

p <- ggplot() +
  geom_sf(data = world_crop, fill = "black", color = "grey70") +
  
   # points/tracks
  geom_sf(data = dt_cum, color = "grey40", size = 1.3) +
  
  # 50 m bathymetry contour
  geom_sf(data = b, color = "deepskyblue", size = 0.6) +
  
  ggtitle(title_txt) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"])
  ) +
  theme_dark() +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    axis.text        = element_text(color = "white", size = 10),
    axis.title       = element_text(color = "white"),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_line(color = "grey25"),
    plot.title       = element_text(color = "white", size = 26, face = "bold"),
    panel.border     = element_blank()
  )

p

# Save
ggsave(
  filename = file.path(output_path, "year_2018-2024_with_50m_bathy.png"),
  plot = p, width = 10, height = 8, dpi = 300, bg = "black"
)













bathy_col <- "deepskyblue"

p <- ggplot() +
  geom_sf(data = world_crop, fill = "black", color = "grey70") +
  
  # points/tracks
  geom_sf(data = dt_cum, color = "grey40", size = 1.3) +
  
  # 50 m bathymetry contour
  geom_sf(data = b, color = bathy_col, size = 0.6) +
  
  ggtitle(title_txt) +
  
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"])
  ) +
  
  # Annotation (top-left)
  annotate(
    "text",
    x = bb["xmin"] + 0.02 * (bb["xmax"] - bb["xmin"]),    # slight margin from left
    y = bb["ymax"] - 0.8 * (bb["ymax"] - bb["ymin"]),    # slight margin from top
    label = "Bathymétrie = 50 mètres",
    color = bathy_col,
    size = 5,
    hjust = 0,
    vjust = 1
  ) +
  
  theme_dark() +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    axis.text        = element_text(color = "white", size = 10),
    axis.title       = element_text(color = "white"),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_line(color = "grey25"),
    plot.title       = element_text(color = "white", size = 26, face = "bold"),
    panel.border     = element_blank()
  )

p











#------------------------------------------------- EMPTY GRID --------------------
grid <- st_read("./data/processed_data/prediction_extent/grid_v1.1_3857.gpkg")
plot(grid)



p_grid <- ggplot(grid) +
  geom_sf(fill = NA, color = "deepskyblue", linewidth = 0.1) +
  # ggtitle("Sites à 0-50m de profondeur") +
  # coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
  #          ylim = c(bb["ymin"], bb["ymax"])) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "white"),
    
    # Black background everywhere
    panel.background = element_rect(fill = "black", color = NA),
    plot.background  = element_rect(fill = "black",  color = NA),
    
    # Grid lines faded to dark grey
    panel.grid.major = element_line(color = "grey30", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    
    # Axis text in white
    axis.text       = element_text(color = "white"),
    axis.title      = element_text(color = "white")
  )

p_grid


ggsave(
  filename = file.path(output_path, "empty_grid_0-50m_depth_area.png"),
  plot = p_grid, width = 10, height = 8, dpi = 300, bg = "black"
)







































#------------------------------------------------- PREDICTION MAPS ----------------------
# --- Richness surface ----
# Load & prep predictions ----
pd <- readRDS("./output/predictions/rf/T1.33/T1.33_R_predictions_surface_2023-av-sept.rds")

# transform to 3857
pd_3857 <- lapply(pd, function(x) st_transform(x, 3857))




# Map ----
library(RColorBrewer)

Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces de poissons en", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_R), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "PuBu"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 20, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.33_prediction_surface_R_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )

}






















# --- Richness bottom ----
# Load & prep predictions ----
pd_b <- readRDS("./output/predictions/rf/T1.33/T1.33_R_predictions_bottom_2023-av-sept.rds")

# transform to 3857
pd_b_3857 <- lapply(pd_b, function(x) st_transform(x, 3857))


# Map ----
Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces de poissons en", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_R), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "PuBu"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 20, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.33_prediction_bottom_R_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )
  
}











# --- Richness bottom ----
# Load & prep predictions ----
pd_b <- readRDS("./output/predictions/rf/T1.33/T1.33_R_predictions_bottom_2023-av-sept.rds")

# transform to 3857
pd_b_3857 <- lapply(pd_b, function(x) st_transform(x, 3857))


# Map ----
Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces de poissons en", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_R), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "PuBu"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 20, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.33_prediction_bottom_R_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )
  
}













# --- Richness surface + bottom panel ----
library(patchwork)
library(lubridate)
library(RColorBrewer)

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

# Common function to build one map (no legend, title = real month in French)
make_map <- function(df, date) {
  d_obj     <- ymd(date) %m-% months(1)          # real month = date - 1 month
  month_num <- format(d_obj, "%m")
  year_txt  <- format(d_obj, "%Y")
  date_fmt  <- paste(mois_fr[[month_num]], year_txt)
  
  ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_R), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "PuBu"),
      name = ""
    ) +
    ggtitle(date_fmt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 8),
      axis.title       = element_blank(),http://127.0.0.1:37451/graphics/plot_zoom_png?width=1850&height=1016
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 14, face = "bold"),
      legend.position  = "none",          # <- remove legend completely
      panel.border     = element_blank()
    )
}

# Build surface & bottom plot lists in the same order of dates
dates <- names(pd_3857)

plots_surface <- lapply(dates, function(d) make_map(pd_3857[[d]],  d))
plots_bottom  <- lapply(dates, function(d) make_map(pd_b_3857[[d]], d))

# Combine: top row = surface, bottom row = bottom (2 x 6)
panel_2x6 <- wrap_plots(
  c(plots_surface, plots_bottom),
  nrow = 2
)

panel_2x6

# Optional: save panel
ggsave(
  file.path(output_path, "T1.33_R_surface_bottom_panel_2x6.png"),
  plot = panel_2x6,
  width = 18, height = 6, dpi = 300, bg = "black"
)












#panel2 ---------------------------------
panel_2x6_labeled <- panel_2x6 +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA)
    )
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  # Add global labels
  patchwork::plot_annotation() &
  theme(
    plot.margin = margin(40, 30, 40, 50),   # space for labels
    axis.text  = element_text(color = "white")
  )

# Add labels manually using patchwork's guide_area trick
library(grid)

final_panel <- (
  panel_2x6 /
    grid::textGrob("Longitude", gp = grid::gpar(col = "white", fontsize = 16))
) +
  grid::textGrob("Latitude", rot = 90, gp = grid::gpar(col = "white", fontsize = 16))

final_panel



plots_surface[[1]] <- plots_surface[[1]] + theme(axis.title.y = element_text(color="white"))
plots_bottom[[1]]  <- plots_bottom[[1]]  + theme(axis.title.y = element_text(color="white"))

for(i in seq_along(plots_bottom)) {
  plots_bottom[[i]] <- plots_bottom[[i]] + theme(axis.title.x = element_text(color="white"))
}



#-------------------------- RedList ------------------------------
# --- RedList surface ----
# Load & prep predictions ----
rm(pd)
pd <- readRDS("./output/predictions/rf/T1.37/T1.37_RedList_predictions_surface_2023-av-sept.rds")

# transform to 3857
pd_3857 <- lapply(pd, function(x) st_transform(x, 3857))




# Map ----
library(RColorBrewer)

Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  # df$prediction_RedList <- round(df$prediction_RedList, 0)
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces menacées en", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_RedList), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "Reds"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 18, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.37_prediction_surface_RedList_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )

}


















# --- RedList bottom ----
# Load & prep predictions ----
rm(pd_b)
pd_b <- readRDS("./output/predictions/rf/T1.37/T1.37_RedList_predictions_bottom_2023-av-sept.rds")

# transform to 3857
pd_b_3857 <- lapply(pd_b, function(x) st_transform(x, 3857))


# Map ----
Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)[1]) {
  
  df <- pd_3857[[date]]
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces de poissons menacés", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_RedList), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "OrRd"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 18, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.33_prediction_bottom_R_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )
  
}














#-------------------------- Crypto ------------------------------
# --- Crypto surface ----
# Load & prep predictions ----
rm(pd)
pd <- readRDS("./output/predictions/rf/T1.37/T1.37_Crypto_predictions_surface_2023-av-sept.rds")

# transform to 3857
pd_3857 <- lapply(pd, function(x) st_transform(x, 3857))




# Map ----
library(RColorBrewer)

Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  # df$prediction_Crypto <- round(df$prediction_Crypto, 0)
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces cachées en", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_Crypto), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "BuPu"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 18, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.37_prediction_surface_Crypto_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )
  
}


















# --- Crypto bottom ----
# Load & prep predictions ----
rm(pd_b)
pd_b <- readRDS("./output/predictions/rf/T1.37/T1.37_Crypto_predictions_bottom_2023-av-sept.rds")

# transform to 3857
pd_b_3857 <- lapply(pd_b, function(x) st_transform(x, 3857))


# Map ----
Sys.setlocale("LC_TIME", "fr_FR.UTF-8")

mois_fr <- c(
  "01" = "janvier", "02" = "février", "03" = "mars",
  "04" = "avril",   "05" = "mai",      "06" = "juin",
  "07" = "juillet", "08" = "août",     "09" = "septembre",
  "10" = "octobre", "11" = "novembre", "12" = "décembre"
)

date_obj <- as.Date(date)
month_num <- format(date_obj, "%m")
year_txt  <- format(date_obj, "%Y")

date_fmt <- paste(mois_fr[[month_num]], year_txt)

for (date in names(pd_3857)) {
  
  df <- pd_3857[[date]]
  
  # date in object is 1 month ahead -> subtract 1 month
  d_obj      <- ymd(date) %m-% months(1)
  month_num  <- format(d_obj, "%m")
  year_txt   <- format(d_obj, "%Y")
  date_fmt   <- paste(mois_fr[[month_num]], year_txt)
  
  title_txt <- paste("Estimation du nombre d'espèces de poissons cachés", date_fmt) 
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "black", color = "grey70") +
    geom_sf(data = df, aes(fill = prediction_Crypto), color = NA) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "BuPu"),
      name = ""
    ) +
    ggtitle(title_txt) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"])
    ) +
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      axis.text        = element_text(color = "white", size = 10),
      axis.title       = element_text(color = "white"),
      panel.grid.major = element_line(color = "grey40"),
      panel.grid.minor = element_line(color = "grey25"),
      plot.title       = element_text(color = "white", size = 18, face = "bold"),
      legend.title     = element_text(color = "white"),
      legend.text      = element_text(color = "white", size = 16),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      panel.border     = element_blank()
    )
  
  print(p)
  
  #Save
  ggsave(
    filename = file.path(
      output_path,
      paste0("T1.33_prediction_bottom_R_", date, ".png")
    ),
    plot = p,
    width = 10, height = 8, dpi = 300, bg = "black"
  )
  
}












#===================================================================
#=================== ARCHIVED CODE ------------------------------------

# #--------------- Map year black --------------------
# years <- 2017:2024
# for (y in years) {
#   dt_cum <- dt %>% filter(year <= y)
#   
#   # dynamic title text
#   if (y == 2017) {
#     title_txt <- "2017"
#   } else if (y == 2018) {
#     title_txt <- "2018"
#   } else {
#     title_txt <- paste0("2018-", y)
#   }
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt_cum, color = "white", size = 1.5) +
#     ggtitle(title_txt) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       
#       # axes + text + ticks in white
#       axis.text = element_text(color = "white", size = 10),
#       axis.title = element_text(color = "white"),
#       
#       # gridlines white or grey
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       
#       # white title
#       plot.title = element_text(color = "white", size = 26, face = "bold"),
#       
#       # remove panel border
#       panel.border = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# #--------------- Map methods (black) ------------------------
# methods <- sort(na.omit(unique(dt$method)))
# cols    <- grDevices::rainbow(length(methods))
# 
# # always same color
# for (i in seq_along(methods)) {
#   m    <- methods[i]
#   dt_m <- dt %>% filter(method == m)
#   
#   method_color <- "#F333FF"   # or cols[i]
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt,    color = "grey40", size = 1.3) +
#     geom_sf(data = dt_m,  color = method_color, size = 1.8) +
#     ggtitle(m) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background  = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       axis.text        = element_text(color = "white", size = 10),
#       axis.title       = element_text(color = "white"),
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       # TITLE COLOR → apply here
#       plot.title       = element_text(color = method_color,
#                                       size = 26, face = "bold"),
#       panel.border     = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# 
# 
# 
# #--------------- Map components (black) ------------------------
# components <- sort(na.omit(unique(dt$component)))
# cols_comp  <- grDevices::rainbow(length(components))
# 
# for (i in seq_along(components)) {
#   m    <- components[i]
#   dt_m <- dt %>% filter(component == m)  # <- use column name `component`
#   
#   #comp_color <- "#33FFF5"   # or cols_comp[i]
#   comp_color <- "gold"
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt,    color = "grey40", size = 1.3) +
#     geom_sf(data = dt_m,  color = comp_color, size = 1.8) +
#     ggtitle(m) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background  = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       axis.text        = element_text(color = "white", size = 10),
#       axis.title       = element_text(color = "white"),
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       plot.title       = element_text(color = comp_color,
#                                       size = 26, face = "bold"),
#       panel.border     = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

# #--------------- Map year black --------------------
# years <- 2017:2024
# for (y in years) {
#   dt_cum <- dt %>% filter(year <= y)
#   
#   # dynamic title text
#   if (y == 2017) {
#     title_txt <- "2017"
#   } else if (y == 2018) {
#     title_txt <- "2018"
#   } else {
#     title_txt <- paste0("2018-", y)
#   }
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt_cum, color = "white", size = 1.5) +
#     ggtitle(title_txt) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       
#       # axes + text + ticks in white
#       axis.text = element_text(color = "white", size = 10),
#       axis.title = element_text(color = "white"),
#       
#       # gridlines white or grey
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       
#       # white title
#       plot.title = element_text(color = "white", size = 26, face = "bold"),
#       
#       # remove panel border
#       panel.border = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# #--------------- Map years ------------------------
# 
# 
# years <- 2018:2024
# 
# for (y in years) {
#   dt_cum <- dt %>% filter(year <= y)
#   
#   # dynamic title
#   if (y == 2018) {
#     title_txt <- "2018"
#   } else {
#     title_txt <- paste0("2018-", y)
#   }
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "grey90", color = "grey60") +
#     geom_sf(data = dt_cum, color = "grey70", size = 1) +
#     ggtitle(title_txt) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_minimal()
#   
#   print(p)
# }
# 
# 
# 
# #--------------- Map methods ------------------------
# library(dplyr)
# library(ggplot2)
# library(sf)
# 
# # # world sf
# # world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# #
# # # bounding box of your points
# # bb <- st_bbox(dt)
# #
# # # crop world to bounding box
# # world_crop <- st_crop(world, bb)
# 
# # vector of methods (drop NA if needed)
# methods <- sort(na.omit(unique(dt$method)))
# 
# # pick one color per method
# cols <- grDevices::rainbow(length(methods))
# 
# for (i in seq_along(methods)) {
#   m <- methods[i]
#   dt_m <- dt %>% filter(method == m)
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "grey90", color = "grey60") +
#     # all points in grey
#     geom_sf(data = dt, color = "grey70", size = 1, alpha = 0.5) +
#     # current method highlighted
#     geom_sf(data = dt_m, color = cols[i], size = 1.5) +
#     ggtitle(paste("Méthode:", m)) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_minimal()
#   
#   print(p)
# }
# 
# 
# 
# 
# 
# #--------------- Map components ------------------------
# library(dplyr)
# library(ggplot2)
# library(sf)
# 
# 
# # vector of components (drop NA if needed)
# components <- sort(na.omit(unique(dt$component)))
# 
# # pick one color per component
# cols <- grDevices::rainbow(length(components))
# 
# for (i in seq_along(components)) {
#   m <- components[i]
#   dt_m <- dt %>% filter(component == m)
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "grey90", color = "grey60") +
#     # all points in grey
#     geom_sf(data = dt, color = "grey70", size = 1, alpha = 0.5) +
#     # current component highlighted
#     geom_sf(data = dt_m, color = cols[i], size = 1.5) +
#     ggtitle(paste("Milieu:", m)) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_minimal()
#   
#   print(p)
# }
# 
# 
# 
# 
# 
# 


# #--------------- Map year black --------------------
# years <- 2017:2024
# for (y in years) {
#   dt_cum <- dt %>% filter(year <= y)
#   
#   # dynamic title text
#   if (y == 2017) {
#     title_txt <- "2017"
#   } else if (y == 2018) {
#     title_txt <- "2018"
#   } else {
#     title_txt <- paste0("2018-", y)
#   }
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt_cum, color = "white", size = 1.5) +
#     ggtitle(title_txt) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       
#       # axes + text + ticks in white
#       axis.text = element_text(color = "white", size = 10),
#       axis.title = element_text(color = "white"),
#       
#       # gridlines white or grey
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       
#       # white title
#       plot.title = element_text(color = "white", size = 26, face = "bold"),
#       
#       # remove panel border
#       panel.border = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# #--------------- Map methods (black) ------------------------
# library(dplyr)
# library(ggplot2)
# library(sf)
# 
# methods <- sort(na.omit(unique(dt$method)))
# cols    <- grDevices::rainbow(length(methods))
# 
# # try other palettes (eg viridis)
# cols <- viridis::viridis(length(methods))
# cols <- viridis::magma(length(methods))
# 
# # manual colors
# cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
# 
# # bright colors 
# cols <- c("#FF5733", "#33FF57", "#3357FF", "#F333FF", "#33FFF5", "#F5FF33")
# 
# 
# 
# # for (i in seq_along(methods)) {
# #   m    <- methods[i]
# #   dt_m <- dt %>% filter(method == m)
# #   
# #   p <- ggplot() +
# #     geom_sf(data = world_crop, fill = "black", color = "grey70") +
# #     # all points in white, faint
# #     geom_sf(data = dt,    color = "grey40", size = 1.3) +
# #     # current method highlighted in color
# #     geom_sf(data = dt_m,  color = cols[i], size = 1.8) +
# #     ggtitle(paste("Méthode:", m)) +
# #     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
# #              ylim = c(bb["ymin"], bb["ymax"])) +
# #     theme_dark() +
# #     theme(
# #       plot.background  = element_rect(fill = "black", color = NA),
# #       panel.background = element_rect(fill = "black", color = NA),
# #       axis.text        = element_text(color = "white", size = 10),
# #       axis.title       = element_text(color = "white"),
# #       panel.grid.major = element_line(color = "grey40"),
# #       panel.grid.minor = element_line(color = "grey25"),
# #       plot.title       = element_text(color = "white", size = 26, face = "bold"),
# #       panel.border     = element_blank()
# #     )
# #   
# #   print(p)
# # }
# 
# 
# 
# # always same color
# for (i in seq_along(methods)) {
#   m    <- methods[i]
#   dt_m <- dt %>% filter(method == m)
#   
#   method_color <- "#F333FF"   # or cols[i]
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt,    color = "grey40", size = 1.3) +
#     geom_sf(data = dt_m,  color = method_color, size = 1.8) +
#     ggtitle(m) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background  = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       axis.text        = element_text(color = "white", size = 10),
#       axis.title       = element_text(color = "white"),
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       # TITLE COLOR → apply here
#       plot.title       = element_text(color = method_color,
#                                       size = 26, face = "bold"),
#       panel.border     = element_blank()
#     )
#   
#   print(p)
# }
# 
# 
# 
# 
# 
# 
# #--------------- Map components (black) ------------------------
# library(dplyr)
# library(ggplot2)
# library(sf)
# 
# components <- sort(na.omit(unique(dt$component)))
# cols_comp  <- grDevices::rainbow(length(components))
# 
# # for (i in seq_along(components)) {
# #   m    <- components[i]
# #   dt_m <- dt %>% filter(component == m)
# #   
# #   p <- ggplot() +
# #     geom_sf(data = world_crop, fill = "black", color = "grey70") +
# #     # all points in white, faint
# #     geom_sf(data = dt,    color = "white", size = 1,   alpha = 0.3) +
# #     # current component highlighted in color
# #     geom_sf(data = dt_m,  color = cols_comp[i], size = 1.8, alpha = 0.9) +
# #     ggtitle(paste("Milieu:", m)) +
# #     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
# #              ylim = c(bb["ymin"], bb["ymax"])) +
# #     theme_dark() +
# #     theme(
# #       plot.background  = element_rect(fill = "black", color = NA),
# #       panel.background = element_rect(fill = "black", color = NA),
# #       axis.text        = element_text(color = "white", size = 10),
# #       axis.title       = element_text(color = "white"),
# #       panel.grid.major = element_line(color = "grey40"),
# #       panel.grid.minor = element_line(color = "grey25"),
# #       plot.title       = element_text(color = "white", size = 26, face = "bold"),
# #       panel.border     = element_blank()
# #     )
# #   
# #   print(p)
# # }
# 
# 
# 
# 
# 
# 
# components <- sort(na.omit(unique(dt$component)))
# cols_comp  <- grDevices::rainbow(length(components))
# 
# for (i in seq_along(components)) {
#   m    <- components[i]
#   dt_m <- dt %>% filter(component == m)  # <- use column name `component`
#   
#   #comp_color <- "#33FFF5"   # or cols_comp[i]
#   comp_color <- "gold"
#   
#   p <- ggplot() +
#     geom_sf(data = world_crop, fill = "black", color = "grey70") +
#     geom_sf(data = dt,    color = "grey40", size = 1.3) +
#     geom_sf(data = dt_m,  color = comp_color, size = 1.8) +
#     ggtitle(m) +
#     coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#              ylim = c(bb["ymin"], bb["ymax"])) +
#     theme_dark() +
#     theme(
#       plot.background  = element_rect(fill = "black", color = NA),
#       panel.background = element_rect(fill = "black", color = NA),
#       axis.text        = element_text(color = "white", size = 10),
#       axis.title       = element_text(color = "white"),
#       panel.grid.major = element_line(color = "grey40"),
#       panel.grid.minor = element_line(color = "grey25"),
#       plot.title       = element_text(color = comp_color,
#                                       size = 26, face = "bold"),
#       panel.border     = element_blank()
#     )
#   
#   print(p)
# }
# 
# 

#----------------------------------------------------- MAP PREDICTIONS





#------------------------------------------------------------- TO DO ------------------------------------

#------------------------------------------------- NUMBER OF SPECIES ------------------------------



#--------------- Load data ------------------------
#--------------- Total number of species  ------------------------
#--------------- Number of species per site = HOTSPOTS ----------------------

















