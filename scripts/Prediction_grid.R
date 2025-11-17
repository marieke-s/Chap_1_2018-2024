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


# Export grid_v1.0_20230701 ---------------
man_grid <- st_read("./data/processed_data/prediction_extent/med_empty-grid_1.2km_2154_manual-sel_50m.gpkg")

# We need to set a date for climatic data extraction. We first want to predict on the month of June 2023 : we set the date to 1/07/2023 for all the polygons :
man_grid$date <- as.Date("2023-07-01")

# We also need to set a depth_samping for marsd extraction. We create 2 depth_sampling cols : 1 for surface and 1 for 40m depth. 
man_grid$depth_sampling_surface <- 1
man_grid$depth_sampling_40m <- 40

st_write(man_grid, "./data/processed_data/prediction_extent/grid_v1.0_20230701.gpkg",
  delete_dsn = TRUE)


#------------------------------------- Prediction extent ------------
#--- Distance to shore ----
distshore <- terra::rast("./data/processed_data/prediction_extent/debug_shore_fraction.tif")
plot(distshore)



#--- Bathymetry 50m extent -----
# load bathy raster 
bathy <- raster::raster("./data/raw_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged.tif")


plot(bathy)

# remove all values below -50 and above 0
bathy50m <- bathy %>%
  raster::calc(fun = function(x) {
    x[x < -50 | x > -0.99] <- NA
    return(x)
  })

plot(bathy50m)

# Export bathy50m
raster::writeRaster(bathy50m, filename = "./data/processed_data/predictors/Bathymetry/MNT_MED_CORSE_SHOM100m_merged-50-0.tif", format = "GTiff", overwrite=TRUE)

# Extent 50m set bathy50m NA to 0 and others to 1
extent50m <- bathy50m %>%
  raster::calc(fun = function(x) {
    x[!is.na(x)] <- 1
    return(x)
  })

plot(extent50m)

raster::writeRaster(extent50m, filename = "./data/processed_data/prediction_extent/extent_bathy50m.tif", format = "GTiff", overwrite=TRUE)





#--- Bathymetry 50 x EEZ Fr + MC ----
# QGIS : we crop the bathy 50m extent on EEZ France + Monaco on Med only
bath50mfr <- terra::rast("./data/processed_data/prediction_extent/extent_bathy50m_ZEE-FR-MC-Med.tif")
plot(bath50mfr)

# Reproject 
bath50mfr <- terra::project(bath50mfr, "EPSG:2154") # Lambert-93





#--- Comput 1.2km empty grid ----
# Based on bathy50mfr extent, we compute a 1.2km grid for prediction extent. 
# 1.2km resolution is chosen based on the mean and median buffer area of our site (see "./figure/Mtdt/Map_Hits/mtdt_7_sel_v1.0_buffer_area").

# 1.2km grid : raster 
r <- rast(ext(bath50mfr), res = 1200, crs = crs(bath50mfr))
values(r) <- 1:ncell(r)           # unique value per cell
plot(r)

# 1.2 grid : vector  
# Each pixel = 1 polygon with 1 ID
grid_vec <- as.polygons(r, aggregate = FALSE, values = TRUE, dissolve = FALSE)

# count nb of polygons in grid_vec
nrow(as.data.frame(grid_vec))
colnames(as.data.frame(grid_vec))

# plot with color by lyr.1
plot(grid_vec, "lyr.1")

# export 
st_write(as_sf(grid_vec), "./test.gpkg")

grid_vec
plot(grid_vec, col = NA, border = "black")

# Plot with color by ID
plot(grid_vec, "ID")





# # Presence/absence in fine bathy
# bath_pres <- !is.na(bath50mfr)

# Vectorise bath50mfr
vect_mer <- terra::as.polygons(bath50mfr, dissolve = TRUE)
plot(vect_mer)
class(vect_mer)

# export vect_mer
terra::writeVector(vect_mer, filename = "./data/processed_data/prediction_extent/prediction_extent_12km_bathy50m_ZEE_FR-MC_Med.gpkg", overwrite=TRUE)


# set vect_mer$frition to 1 inside polygons and to 0 outside


dim(vect_mer)

# rasterisation avec somme de couverture (fraction d'occupation)
r_frac <- terra::rasterize(vect_mer, r, field = "friction",
                           fun = "sum", background = 0, touches = TRUE)

plot(r_frac)

seuil_couverture <- 0.8

r_frac <- r_frac / max(r_frac[], na.rm = TRUE)  # application du seuil
values(r_frac) <- ifelse(values(r_frac) >= seuil_couverture, 1, NA)


plot(r_frac)













# Resample fine → coarse grid using max (any presence = 1)
r_masked <- resample(bath_pres, r, method = "max")

plot(r_masked)

# Export
terra::writeRaster(r_masked, filename = "./data/processed_data/prediction_extent/prediction_extent_12km_bathy50m_ZEE_FR-MC_Med.tif", overwrite=TRUE)

#--- Vectorize the grid ---
# Vectorize to convert to a gpkg file 
r_vect <- terra::as.polygons(r_masked, dissolve = TRUE)

plot(r_vect)

# Export
terra::writeVector(r_vect, filename = "./data/processed_data/prediction_extent/prediction_extent_12km_bathy50m_ZEE_FR-MC.gpkg", overwrite=TRUE)



d