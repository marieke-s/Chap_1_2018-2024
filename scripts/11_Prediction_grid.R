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

  
  # file name
  fname <- paste0("./data/processed_data/prediction_extent/grid_v1.0_", d, ".gpkg")
  
  # Write file
  st_write(man_grid, fname, delete_dsn = TRUE)
}
