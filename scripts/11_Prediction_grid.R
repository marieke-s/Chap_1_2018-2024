#------------- Description ----
# The aim of this script is to create an empty grid to use for predictions. 





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

















#------------- Load data --------------
bathy50m <- terra::rast("./data/processed_data/prediction_extent/extent_bathy50m_ZEE-FR-MC-Med.tif")

dist_shore <- terra::rast("./data/processed_data/prediction_extent/debug_shore_fraction.tif")





#--------------------------------------------------- PREDICTION GRID V.1.0 ----------------
#------------ Generate Med empty grid EPSG:2154 -------------

# 1. Original bbox in WGS84
xmin <- 1.7126959079564994
ymin <- 41.14251524871988
xmax <- 10.809848598828735
ymax <- 44.53082479002648

bb_4326 <- st_bbox(
  c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  crs = st_crs(4326)
)

# 2. Transform bbox to EPSG:2154
bb_2154_poly <- st_as_sfc(bb_4326) |>
  st_transform(2154)

bb_2154 <- st_bbox(bb_2154_poly)

# Drop attributes to get plain numerics
xmin_2154 <- as.numeric(bb_2154["xmin"])
ymin_2154 <- as.numeric(bb_2154["ymin"])
xmax_2154 <- as.numeric(bb_2154["xmax"])
ymax_2154 <- as.numeric(bb_2154["ymax"])

# 3. Cell size = 1.2 km (in meters)
cellsize <- 1200

# 4. Compute number of cells (rounding UP)
x_range <- xmax_2154 - xmin_2154
y_range <- ymax_2154 - ymin_2154

nx <- ceiling(x_range / cellsize)
ny <- ceiling(y_range / cellsize)

# 5. Expanded bbox that fits an integer number of cells
xmax_expanded <- xmin_2154 + nx * cellsize
ymax_expanded <- ymin_2154 + ny * cellsize

bb_2154_expanded <- st_bbox(
  c(
    xmin = xmin_2154,
    ymin = ymin_2154,
    xmax = xmax_expanded,
    ymax = ymax_expanded
  ),
  crs = st_crs(2154)
)

# 6. Create grid (all cells exactly 1200m x 1200m)
grid_2154 <- st_make_grid(
  bb_2154_expanded,
  cellsize = cellsize,
  square   = TRUE
) |>
  st_sf() |>
  mutate(
    id    = dplyr::row_number(),
    value = 0
  )


plot(st_geometry(grid_2154))
beepr::beep()

# 8. Export whatever you prefer: full expanded grid or clipped
st_write(grid_2154, "./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154.gpkg",
  delete_dsn = TRUE)


#------------ Crop empty grid on coastal extent -------------
# 1. Extract with area of intersecting pieces
crs(bathy50m)  # should be EPSG:2154

# reproject bathy50m to EPSG:2154 if needed
bathy50m <- terra::project(bathy50m, "EPSG:2154") # Lambert-93

# Extract with area of intersecting pieces
extraction <- exact_extract(
  bathy50m,
  grid_2154,
  include_area = TRUE   # column "area" = area of each cell-polygon piece
)
beepr::beep()

# Area of each grid cell (in m^2)
cell_area <- as.numeric(st_area(grid_2154))

# Compute % coverage using ONLY pixels where bathy50m == 1 (and non-NA)
coverage_fraction <- vapply(
  seq_along(extraction),
  function(i) {
    df <- extraction[[i]]
    if (is.null(df) || nrow(df) == 0) return(0)
    
    # keep only pixels where value == 1 and not NA
    df1 <- df[!is.na(df$value) & df$value == 1, , drop = FALSE]
    
    if (nrow(df1) == 0) return(0)
    
    covered_area <- sum(df1$area, na.rm = TRUE)
    covered_area / cell_area[i]
  },
  numeric(1)
)
beepr::beep()

# Add coverage to grid and keep only polygons with >= x% coverage
grid_2154 <- grid_2154 %>%
  mutate(coverage = coverage_fraction)

x = 0.20

grid_2154_xpc <- grid_2154 %>%
  filter(coverage >=x)

# Optional: check
summary(grid_2154_xpc$coverage)

# Optional: visualize
plot(st_geometry(grid_2154_xpc))

# 5. Export
st_write(grid_2154_10pc, "./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154_20pc-bathy50m.gpkg",
  delete_dsn = TRUE)


#------------ Manual cell selection QGIS ------------
# On QGIS we manually select cells to include in prediction extent based on med_empty-grid_1.2km_2154_10pc-bathy50m.gpkg and simplifying the shape for visualisation purpose.

# man_grid <- st_read("./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154_manual-sel_50m.gpkg")


#------------ Export grid_v1.0_ april to september 2023---------------
# Path to the base grid
grid_path <- "./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154_manual-sel_50m.gpkg"

# Months to generate (always take next month on 1st day so starts 1rst of june for may and ends 1rst of october for september)
dates_to_do <- c("2023-05-01", "2023-08-01", "2023-09-01", "2023-06-01", "2023-07-01", "2023-10-01")

for (d in dates_to_do) {
  print(paste("----------------", "Processing date:", d, "----------------"))

  # Read the base grid
  man_grid <- st_read(grid_path, quiet = TRUE)
  
  # Set date and depth sampling columns
  man_grid$date <- as.Date(d)
  man_grid$depth_sampling_surface <- 1
  man_grid$depth_sampling_40m <- 40
  
  # Reproject in 4326
  man_grid <- sf::st_transform(man_grid, crs = 4326)

  # file name
  fname <- paste0("./data/processed_data/prediction_extent/grid_v1.0_", d, ".gpkg")
  
  # Write file
  st_write(man_grid, fname, delete_dsn = TRUE)
}


#------------ Convert grids to .geojson ---------------
# install.packages("sf")   # run once if you don't have sf installed
library(sf)

# 1. Set this to your repo path (or use getwd() if you're already in it)
repo_dir <- "./data/processed_data/prediction_extent/"   # <-- change this

# 2. Find all .gpkg files recursively
gpkg_files <- list.files(
  path       = repo_dir,
  pattern = "^grid_v1\\.0_.*\\.gpkg$",
  recursive  = TRUE,
  full.names = TRUE
)

# Safety check
if (length(gpkg_files) == 0) {
  message("No .gpkg files found in: ", repo_dir)
} else {
  message("Found ", length(gpkg_files), " .gpkg file(s).")
}

# 3. Loop over each .gpkg and convert to .geojson
for (gpkg_path in gpkg_files) {
  message("Processing: ", gpkg_path)
  
  # Read the GeoPackage
  g <- st_read(gpkg_path, quiet = TRUE)
  
  # Build output path: same directory, same basename, .geojson extension
  geojson_path <- sub("\\.gpkg$", ".geojson", gpkg_path)
  
  # Write GeoJSON
  st_write(
    g,
    geojson_path,
    driver   = "GeoJSON",
    delete_dsn = TRUE  # overwrite if file already exists
  )
  
  message("  -> Saved: ", geojson_path)
}

message("Done.")

#------------ grid_v1.1 : april to september 2023 with cols : regions, dates, surface, 40m. ---------------
# Path to the base grid
grid_path <- "./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154_manual-sel_50m.gpkg"
grid <- st_read(grid_path, quiet = TRUE)

grid_corse <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01/grille_corse.geojson")
grid_medest <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01/grille_med_est.geojson")
grid_medouest <- st_read("./data/processed_data/predictors/Prediction_grid_v1.0/2023-07-01/grille_med_ouest.geojson")

plot(grid_medest[1])
plot(grid_medouest[1])
plot(grid_corse[1])

# Create region column
grid <- grid %>%
  mutate(
    region = case_when(
      id %in% grid_corse$id     ~ "corse",
      id %in% grid_medest$id    ~ "med-est",
      id %in% grid_medouest$id  ~ "med-ouest",
      TRUE ~ NA_character_      # fallback (should not happen)
    )
  )

plot(grid[4])

# Create date columns 
grid <- grid %>%
  mutate(
    '2023-05-01' = "2023-05-01",
    '2023-06-01' = "2023-06-01",
    '2023-07-01' = "2023-07-01",
    '2023-08-01' = "2023-08-01",
    '2023-09-01' = "2023-09-01",
    '2023-10-01' = "2023-10-01"
  )

# Create depth_sampling columns
grid <- grid %>%
  mutate(
    depth_sampling_surface = 1,
    depth_sampling_40m     = 40
  )

# Reproject in 4326
grid <- sf::st_transform(grid, crs = 4326)

# Export grid 
st_write(grid, "./data/processed_data/prediction_extent/grid_v1.1.gpkg",
  delete_dsn = TRUE)

str(grid)





#------------ grid_v1.1_3857 -----------------
grid <- st_read("./data/processed_data/prediction_extent/grid_v1.1.gpkg", quiet = TRUE)
names(grid)

grid_3857 <- sf::st_transform(grid, crs = 3857)
st_crs(grid_3857)

st_write(grid_3857, "./data/processed_data/prediction_extent/grid_v1.1_3857.gpkg",
  delete_dsn = TRUE)



#------------ Check extraction ----------------
# Read all .geojson in the folder 
# grid_files <- list.files("./data/processed_data/predictors/Prediction_grid_v1.0/CUR-WIND/", pattern = "*\\.geojson$", full.names = TRUE)
grid_files <- list.files("/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_MARS3D_2018-2024/grid_v1.0/CUR-WIND", pattern = "*\\.geojson$", full.names = TRUE)

all_grids <- lapply(grid_files, st_read, quiet = TRUE)

t <- all_grids[[1]]
colnames(all_grids[[1]])
unique(t$date)

# Print a summary of each grid
for (i in seq_along(all_grids)) {
  cat("Summary of grid file:", grid_files[i], "\n")
  # make all cols numeric 
  all_grids[[i]] <- all_grids[[i]] %>%
    mutate(across(everything(), as.numeric)) %>%
    st_drop_geometry() %>%
    dplyr::select(-c("id",            "value",        "date",           "depth_sampling_surface", "depth_sampling_40m" ))
  print(summary(all_grids[[i]]))
  cat("\n")
}

# Check for NAs in each grid and print counts per column
for (i in seq_along(all_grids)) {
  if (anyNA(all_grids[[i]])) {
    print(paste("------------------- NA in :",  grid_files[i], "---------------"))
    cat("NA counts per column:\n")
    print(colSums(is.na(all_grids[[i]])))
  }
}

#----------------- FOR CSV (CHL/SST)
# Lister les fichiers CSV
csv_files <- list.files(
  "/run/user/1000/gvfs/sftp:host=marbec-data.ird.fr/BiodivMed/output/Extraction_SST_2018-2024/grid_v1.0",
  pattern = "\\.csv$", 
  full.names = TRUE
)

# Lire tous les CSV
all_tables <- lapply(csv_files, read.csv)

# Check first
t <- all_tables[[1]]
colnames(t)
if ("date" %in% colnames(t)) unique(t$date)

# Résumé de chaque CSV
for (i in seq_along(all_tables)) {
  cat("Summary of CSV file:", csv_files[i], "\n")
  
  # Colonnes à retirer si elles existent
  cols_to_drop <- c("id", "value", "date", 
                    "depth_sampling_surface", "depth_sampling_40m")
  cols_to_drop <- intersect(cols_to_drop, colnames(all_tables[[i]]))
  
  # Rendre toutes les colonnes numériques puis enlever les colonnes indésirables
  all_tables[[i]] <- all_tables[[i]] %>%
    mutate(across(everything(), as.numeric)) %>%
    dplyr::select(-all_of(cols_to_drop))
  
  print(summary(all_tables[[i]]))
  cat("\n")
}

# Vérifier les NA dans chaque CSV
for (i in seq_along(all_tables)) {
  if (anyNA(all_tables[[i]])) {
    cat("------------------- NA in :",  csv_files[i], "---------------\n")
    cat("NA counts per column:\n")
    print(colSums(is.na(all_tables[[i]])))
  }
}


# bind to geom and export to check 





#----- CHeck boats extraction -----
corse2 <- readr::read_csv("./data/processed_data/predictors/Prediction_grid_v1.1/mtdt_7_boats_month_corse_3857.csv")
corse1 <- readr::read_csv("./data/processed_data/predictors/Prediction_grid_v1.1/mtdt_7_boats_month_corse.csv")

# correlation
corse1 <- corse1 %>%
  # keep replicates that are in corse2
  filter(replicates %in% corse2$replicates)

cor(corse1$Boat_density_month, corse2$Boat_density_month)
names(corse1)
