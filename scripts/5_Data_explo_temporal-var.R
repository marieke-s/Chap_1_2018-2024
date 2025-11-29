#------------- Description ----
# The aim of this script is to check temporal variations in the predictors





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
library(rnaturalearth)



# Load functions
source("./utils/Fct_Data-Prep.R")





#------------- Data prepration-----


# Load data 
pred <- st_read("./data/processed_data/predictors/predictors_raw_v3.0.gpkg")
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv")
occ <- readr::read_csv("./data/processed_data/eDNA/occ_pooled_v1.1.csv")


# Merge data
full <- mtdt %>%
  left_join(
    div,
    by = "replicates"
  ) %>%
  left_join(
    pred %>% st_drop_geometry(),
    by = "replicates"
  ) %>%
  left_join(
    occ, 
    by = "replicates"
  )

# Add season, month, year, yearmonth, yearseason columns
full$date <- lubridate::as_date	(full$date)
full <- full %>%
  mutate(
    year  = year(date),
    month = month(date),
    season = case_when(
      month %in% c(1, 2, 3)  ~ "Winter",
      month %in% c(4, 5, 6)  ~ "Spring",
      month %in% c(7, 8, 9)  ~ "Summer",
      month %in% c(10,11,12) ~ "Autumn",
      TRUE ~ NA_character_
    ),
    yearmonth  = format(date, "%Y%m"),
    yearseason = paste0(year, "_", season)   # <-- HERE
  ) %>%
  mutate(
    month = factor(month.abb[month], levels = month.abb)
  )


full$season <- factor(full$season,
                      levels = c("Winter", "Spring", "Summer", "Autumn"))










#------------- Samples SST and CHL temporal variations -----
# CHL (samples)
chl_plot <- plot_temporal_var(full, date, cop_chl_month_mean,
                                output_file = "./figures/chl_interannual.png")

chl_plot

# Temperature at sampling depth
temp_plot <- plot_temporal_var(full, date, temp_mean_1m,
                                 output_file = "./figures/temp_interannual.png")


temp_plot

# Salinity at sampling depth
sal_plot <- plot_temporal_var(full, date, sal_mean_1m,
                                output_file = "./figures/sal_interannual.png")

sal_plot 

# SST
sst_plot <- plot_temporal_var(full, date, cop_analysed_sst_month_mean,
                                output_file = "./figures/sst_interannual.png")

sst_plot




#------------- Check CHL outliers -------------


# Samples in June 2024 ------------------
f <- full %>%
  filter(year==2024 & month=="Jun")

# Make a map of chl in June 2024
ggplot() +
  geom_sf(data = f, aes(color = cop_chl_month_mean), size = 3) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Chlorophyll-a in June 2024",
    color = "Chl-a (mg/m³)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )


f <- f %>% filter(cop_chl_month_mean > 3)
f %>% pull(date) # 2024-06-06


# Map chl ncdf for June 6, 2024 ----------

# 1. Open NetCDF and read coordinates + time

nc_file <- "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc"
nc <- nc_open(nc_file)

lon  <- ncvar_get(nc, "longitude")   # length 517
lat  <- ncvar_get(nc, "latitude")    # length 277
time <- ncvar_get(nc, "time")        # seconds since 1981-01-01

origin <- as.POSIXct("1981-01-01", tz = "UTC")
dates  <- as.Date(origin + time)     # one per time slice


# 2. Select July 2023 indices

idx_june23 <- which(format(dates, "%Y-%m-%d") == "2023-06-06")
length(idx_june23)    # should be 31


# 3. Read CHL for those days and compute monthly mean
#    CHL[longitude, latitude, time]

chl_sub <- ncvar_get(
  nc,
  varid = "CHL",
  start = c(1, 1, min(idx_june23)),
  count = c(-1, -1, length(idx_june23))
)
nc_close(nc)

# Replace fill value with NA
chl_sub[chl_sub == -999] <- NA

# Mean over time dimension (3rd dim) → matrix [lon, lat]
chl_june_mean <- apply(chl_sub, c(1, 2), mean, na.rm = TRUE)


# 4. Convert to data frame for ggplot

grid <- expand.grid(
  lon = lon,
  lat = lat
)

grid$chl_mean <- as.vector(chl_june_mean)

# Optionally drop NAs
grid <- grid %>% filter(!is.na(chl_mean))


# 5. Build map

world <- ne_countries(scale = "medium", returnclass = "sf")

gg_chl_july23 <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "white") +
  geom_raster(
    data = grid,
    aes(x = lon, y = lat, fill = chl_mean),
    interpolate = FALSE
  ) +
  scale_fill_viridis_c(
    name = "CHL (mg/m³)",
    option = "viridis"
  ) +
  coord_sf(
    xlim = range(grid$lon),
    ylim = range(grid$lat),
    expand = FALSE
  ) +
  labs(
    title = "Mean Chlorophyll-a – June 2023",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

gg_chl_july23



# --> We observe a chl bloom in that region on that date, which explains the high value in the sample.





# Temporal variations of samples CHL -----
# Interannual seasonal variation of chl-a 
chl1 <- ggplot(full, aes(x = date, y = cop_chl_month_mean, group = year, color = as.factor(year))) +
  geom_line(alpha = 0.5, size = 1) +  # Transparent lines for trends
  geom_point(size = 2, alpha = 0.7) +  # Smaller points with transparency
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 1.2, linetype = "dashed") +  # Smooth trend line
  labs(
    title = "Interannual Variation of Chlorophyll-a",
    x = "Date",
    y = "Chlorophyll-a (mg/l³)",
    color = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Better x-axis readability
  scale_color_viridis_d()  # Use colorblind-friendly colors

chl1

# Save as high-resolution PNG
# ggsave("interannual_variation_chl.png", chl, width = 12, height = 10, dpi = 300)


# Temporal variations of CHL ncdf -----
library(ncdf4)
library(dplyr)
library(ggplot2)
library(lubridate)
library(viridis)


# 1. Open NetCDF + read time axis

nc_file <- "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc"
nc <- nc_open(nc_file)

time_raw <- ncvar_get(nc, "time")   # seconds since 1981-01-01
origin   <- as.POSIXct("1981-01-01", tz = "UTC")
dates    <- as.Date(origin + time_raw)

# Check range
range(dates)


# 2. Restrict to 2018–2023

start_date <- as.Date("2018-01-01")
end_date   <- as.Date("2024-12-31")

idx <- which(dates >= start_date & dates <= end_date)

length(idx)   # number of days used


# 3. Compute spatial mean CHL for each day (over full grid)
#    CHL[lon, lat, time]

n_days <- length(idx)
chl_mean_daily <- numeric(n_days)

for (j in seq_along(idx)) {
  t_idx <- idx[j]
  
  # read one time slice: full lon/lat for that day
  chl_slice <- ncvar_get(
    nc,
    varid  = "CHL",
    start  = c(1, 1, t_idx),
    count  = c(-1, -1, 1)
  )
  
  # handle fill values
  chl_slice[chl_slice == -999] <- NA
  
  # spatial mean for that day
  chl_mean_daily[j] <- mean(chl_slice, na.rm = TRUE)
}
    
  nc_close(nc)

# Map check for one day (optional)
#image(chl_slice, col = viridis::viridis(100), main = paste("CHL on", dates[idx][1]))



# 4. Build data frame in the same spirit as 'full'

full <- tibble(
  date                = dates[idx],
  year                = year(date),
  cop_chl_month_mean  = chl_mean_daily   # name it like your existing code
)

# If you really want monthly means (instead of daily), you can do:
# full <- full %>%
#   mutate(month = floor_date(date, "month")) %>%
#   group_by(year, month) %>%
#   summarise(cop_chl_month_mean = mean(cop_chl_month_mean, na.rm = TRUE),
#             .groups = "drop") %>%
#   rename(date = month)


# 5. Plot: same style as your 'chl' example

chl <- ggplot(full, aes(x = date,
                        y = cop_chl_month_mean,
                        group = year,
                        color = as.factor(year))) +
  geom_line(alpha = 0.5, size = 1) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess",
              se = FALSE,
              color = "black",
              size = 1.2,
              linetype = "dashed") +
  labs(
    title = "Interannual Variation of Chlorophyll-a (Domain Mean)",
    x = "Date",
    y = "Chlorophyll-a (mg/m³)",
    color = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_color_viridis_d()


chl
  
  
# chl prediction grid 06/2023-----
grid_chl <- readr::read_csv("./data/processed_data/predictors/Prediction_grid_v1.1/Cop_CHL_stats_2023-07-01_FULL.csv")
grid <- st_read("./data/processed_data/prediction_extent/grid_v1.1.gpkg")

# Merge both to have chl and geometry
grid_chl_sf <- grid %>%
left_join(
  grid_chl,
  by = "id")

names(grid_chl_sf)

# Plot chl on map
ggplot() +
geom_sf(data = grid_chl_sf, aes(fill = Cop_CHL_month_mean), color = NA) +
scale_fill_viridis_c(option = "plasma", na.value = "lightgrey") +
labs(
  title = "Chlorophyll-a on 2023-06-01",
  fill  = "Chl-a (mg/m³)"
) +
theme_minimal() +
theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12)
)


# Make a summary 
# summary(grid_chl_sf$Cop_CHL_month_mean)


# nc chl 06/2023 -----



# 1. Open NetCDF and read coordinates + time

nc_file <- "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc"
nc <- nc_open(nc_file)

lon  <- ncvar_get(nc, "longitude")   # length 517
lat  <- ncvar_get(nc, "latitude")    # length 277
time <- ncvar_get(nc, "time")        # seconds since 1981-01-01

origin <- as.POSIXct("1981-01-01", tz = "UTC")
dates  <- as.Date(origin + time)     # one per time slice


# 2. Select July 2023 indices

idx_june23 <- which(format(dates, "%Y-%m") == "2023-06")
length(idx_june23)    # should be 31


# 3. Read CHL for those days and compute monthly mean
#    CHL[longitude, latitude, time]

chl_sub <- ncvar_get(
nc,
varid = "CHL",
start = c(1, 1, min(idx_june23)),
count = c(-1, -1, length(idx_june23))
)
nc_close(nc)

# Replace fill value with NA
chl_sub[chl_sub == -999] <- NA

# Mean over time dimension (3rd dim) → matrix [lon, lat]
chl_june_mean <- apply(chl_sub, c(1, 2), mean, na.rm = TRUE)


# 4. Convert to data frame for ggplot

grid <- expand.grid(
lon = lon,
lat = lat
)

grid$chl_mean <- as.vector(chl_june_mean)

# Optionally drop NAs
grid <- grid %>% filter(!is.na(chl_mean))


# 5. Build map

world <- ne_countries(scale = "medium", returnclass = "sf")

gg_chl_july23 <- ggplot() +
geom_sf(data = world, fill = "grey90", color = "white") +
geom_raster(
  data = grid,
  aes(x = lon, y = lat, fill = chl_mean),
  interpolate = FALSE
) +
scale_fill_viridis_c(
  name = "CHL (mg/m³)",
  option = "viridis"
) +
coord_sf(
  xlim = range(grid$lon),
  ylim = range(grid$lat),
  expand = FALSE
) +
labs(
  title = "Mean Chlorophyll-a – June 2023",
  x = "Longitude",
  y = "Latitude"
) +
theme_minimal() +
theme(
  panel.background = element_rect(fill = "lightblue"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
)

gg_chl_july23
  
  
  
# Samples CHL with ncdf CHL min and max -----

library(ggplot2)
library(viridis)

chl_envelope_plot <- ggplot() +
  # NCDF envelope: daily min–max across all grid cells
  geom_ribbon(
    data = chl_nc_range,
    aes(x = date, ymin = chl_min, ymax = chl_max),
    fill  = "grey50",
    alpha = 0.4
  ) +
  # Sample-based CHL: same aesthetics as before
  geom_line(
    data  = full,
    aes(x = date, y = cop_chl_month_mean, group = year, color = as.factor(year)),
    alpha = 0.5, size = 1
  ) +
  geom_point(
    data  = full,
    aes(x = date, y = cop_chl_month_mean, color = as.factor(year)),
    size  = 2, alpha = 0.7
  ) +
  geom_smooth(
    data  = full,
    aes(x = date, y = cop_chl_month_mean),
    method = "loess", se = FALSE,
    color  = "black", size = 1.2, linetype = "dashed"
  ) +
  labs(
    title = "Samples vs. NCDF Daily Chlorophyll Range",
    subtitle = "Grey band = daily min–max across all Med (Copernicus CHL)",
    x = "Date",
    y = "Chlorophyll-a (mg/m³)",
    color = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle= element_text(size = 12, hjust = 0.5),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_color_viridis_d()

chl_envelope_plot















# summary stats chl-----
  
# summary of chl_month_mean per year 
library(dplyr)
library(sf)

summary_chl <- full %>%
  st_drop_geometry() %>%        # <--- important
  group_by(year) %>%
  summarise(
    mean_chl = mean(cop_chl_month_mean, na.rm = TRUE),
    sd_chl   = sd(cop_chl_month_mean, na.rm = TRUE),
    min_chl  = min(cop_chl_month_mean, na.rm = TRUE),
    max_chl  = max(cop_chl_month_mean, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

summary_chl

  
  
  
  
  
# max value from chl ncdf-----
nc <- nc_open("./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc")
time <- ncvar_get(nc, "time")
origin <- as.POSIXct("1981-01-01", tz = "UTC")
dates <- as.Date(origin + time)

sel <- dates >= as.Date("2018-01-01") & dates <= as.Date("2024-12-31")
idx  <- which(sel)

# max value over all pixels & all days in that period
max_pixel <- sapply(idx, function(ti) {
  x <- ncvar_get(nc, "CHL", start = c(1,1,ti), count = c(-1,-1,1))
  x[x == -999] <- NA
  max(x, na.rm = TRUE)
  })

nc_close(nc)

max(max_pixel, na.rm = TRUE)



# Summary chl values from ncdf-----
library(ncdf4)
library(dplyr)
library(lubridate)

# ---- open file 
nc <- nc_open("./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc")

time <- ncvar_get(nc, "time")
origin <- as.POSIXct("1981-01-01", tz = "UTC")
dates <- as.Date(origin + time)

# restrict to 2018–2023
sel <- dates >= as.Date("2018-01-01") & dates <= as.Date("2024-12-31")
idx  <- which(sel)

# ---- compute daily spatial means
daily_means <- sapply(idx, function(ti) {
  x <- ncvar_get(nc, "CHL", start = c(1,1,ti), count = c(-1,-1,1))
  x[x == -999] <- NA
  mean(x, na.rm = TRUE)

})

nc_close(nc)

# ---- build df -
chl_df <- tibble(
  date = dates[idx],
  year = year(date),
  chl  = daily_means
)

# yearly CHL summary  ----
summary_chl <- chl_df %>%
  group_by(year) %>%
  summarise(
    mean_chl = mean(chl, na.rm = TRUE),
    sd_chl   = sd(chl, na.rm = TRUE),
    min_chl  = min(chl, na.rm = TRUE),
    max_chl  = max(chl, na.rm = TRUE),
    n        = n()
  )

summary_chl







  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  













#------------- Temporal variations of sst ncdf -----
library(ncdf4)
library(dplyr)
library(ggplot2)
library(lubridate)
library(viridis)


# 1. Open NetCDF + read time axis ----

nc_file <- "./data/raw_data/predictors/SST/SST_MED_SST_L4_NRT_OBSERVATIONS_010_004_c_V2_SST_20130101-20250101.nc"
nc <- nc_open(nc_file)

time_raw <- ncvar_get(nc, "time")   # seconds since 1981-01-01
origin   <- as.POSIXct("1981-01-01", tz = "UTC")
dates    <- as.Date(origin + time_raw)

# Check range
range(dates)


# 2. Restrict to 2018–2023 ----

start_date <- as.Date("2018-01-01")
end_date   <- as.Date("2024-12-31")

idx <- which(dates >= start_date & dates <= end_date)

length(idx)   # number of days used


# 3. Compute spatial mean sst for each day (over full grid) -----
#    sst[lon, lat, time]

n_days <- length(idx)
sst_mean_daily <- numeric(n_days)

for (j in seq_along(idx)) {
  t_idx <- idx[j]
  
  # read one time slice: full lon/lat for that day
  sst_slice <- ncvar_get(
    nc,
    varid  = "analysed_sst",
    start  = c(1, 1, t_idx),
    count  = c(-1, -1, 1)
  )
  
  # handle fill values
  sst_slice[sst_slice == -999] <- NA
  
  # spatial mean for that day
  sst_mean_daily[j] <- mean(sst_slice, na.rm = TRUE)

}
nc_close(nc)

# Map check for one day (optional)
#image(sst_slice, col = viridis::viridis(100), main = paste("sst on", dates[idx][1]))



# 4. Build data frame in the same spirit as 'full' ----

fullNC <- tibble(
  date                = dates[idx],
  year                = year(date),
  cop_sst_month_mean  = sst_mean_daily   # name it like your existing code
)

# If you really want monthly means (instead of daily), you can do:
# full <- full %>%
#   mutate(month = floor_date(date, "month")) %>%
#   group_by(year, month) %>%
#   summarise(cop_sst_month_mean = mean(cop_sst_month_mean, na.rm = TRUE),
#             .groups = "drop") %>%
#   rename(date = month)


# 5. Plot: same style as your 'sst' example ----

sst <- ggplot(fullNC, aes(x = date,
                          y = cop_sst_month_mean,
                          group = year,
                          color = as.factor(year))) +
  geom_line(alpha = 0.5, size = 1) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess",
              se = FALSE,
              color = "black",
              size = 1.2,
              linetype = "dashed") +
  labs(
    title = "Interannual Variation of sst (Domain Mean)",
    x = "Date",
    y = "sst (kelvin)",
    color = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_color_viridis_d()


sst



#------------- Temporal variations of sst ncdf -----
fullNC %>%
  filter(year == 2024) %>%
  pull(cop_sst_month_mean) %>%
  summary()






full %>%
  filter(year == 2024) %>%
  pull(cop_analysed_sst_month_mean) %>%
  summary()





#------------- Temporal variations of samples Temp -----
# Interannual seasonal variation of temp 
temp <- ggplot(full, aes(x = date, y = temp_mean_1m, group = year, color = as.factor(year))) +
  geom_line(alpha = 0.5, size = 1) +  # Transparent lines for trends
  geom_point(size = 2, alpha = 0.7) +  # Smaller points with transparency
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 1.2, linetype = "dashed") +  # Smooth trend line
  labs(
    title = "Interannual Variation of Temperature at sampling depth",
    x = "Date",
    y = "Temperature (°C)",
    color = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Better x-axis readability
  scale_color_viridis_d()  # Use colorblind-friendly colors

temp

# Save as high-resolution PNG
# ggsave("interannual_variation_chl.png", chl, width = 12, height = 10, dpi = 300)



library(ncdf4)
library(dplyr)
library(lubridate)
library(tibble)

nc_file <- "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc"
nc <- nc_open(nc_file)

time <- ncvar_get(nc, "time")
origin <- as.POSIXct("1981-01-01", tz = "UTC")
dates <- as.Date(origin + time)

# Restrict to same period as samples
sel <- dates >= as.Date("2018-01-01") & dates <= as.Date("2024-12-31")
idx <- which(sel)

daily_stats <- t(sapply(idx, function(ti) {
  x <- ncvar_get(nc, "CHL", start = c(1,1,ti), count = c(-1,-1,1))
  x[x == -999] <- NA
  c(
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}))

nc_close(nc)

chl_nc_range <- tibble(
  date    = dates[idx],
  year    = year(date),
  chl_min = daily_stats[, "min"],
  chl_max = daily_stats[, "max"]
)















#------------------------------------- Data prep ------------------

# Load 4 ncdf (sal-temp ; sst ; chl ; wind-vel)
sst <- nc_open("./data/raw_data/predictors/SST/SST_MED_SST_L4_NRT_OBSERVATIONS_010_004_c_V2_SST_20130101-20250101.nc")
chl <- nc_open("./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc")
grid <- st_read("./data/processed_data/prediction_extent/grid_v1.1.gpkg")

grid
sst
chl

plot(grid)
#------------------------------------- INTERANNUAL VARIATIONS OF PREDICTORS -----

library(ncdf4)
library(dplyr)
library(lubridate)
library(sf)
library(tibble)

extract_nc_ts_over_grid <- function(nc_file, var_name, grid, 
                                    start_date = "2018-01-01",
                                    end_date   = "2024-12-31") {
  nc <- nc_open(nc_file)
  
  # coords + time
  lon   <- ncvar_get(nc, "longitude")
  lat   <- ncvar_get(nc, "latitude")
  time  <- ncvar_get(nc, "time")
  origin <- as.POSIXct("1981-01-01", tz = "UTC")
  dates  <- as.Date(origin + time)
  
  # restrict time
  sel_time <- which(dates >= as.Date(start_date) & dates <= as.Date(end_date))
  dates_sel <- dates[sel_time]
  
  # subset lon/lat to grid bbox
  bb <- st_bbox(grid)
  lon_idx <- which(lon >= bb["xmin"] & lon <= bb["xmax"])
  lat_idx <- which(lat >= bb["ymin"] & lat <= bb["ymax"])
  
  # fill value
  fv <- ncatt_get(nc, var_name, "_FillValue")$value
  
  vals <- numeric(length(sel_time))
  
  for (i in seq_along(sel_time)) {
    ti <- sel_time[i]
    x <- ncvar_get(
      nc, var_name,
      start = c(min(lon_idx), min(lat_idx), ti),
      count = c(length(lon_idx), length(lat_idx), 1)
    )
    x[x == fv] <- NA
    vals[i] <- mean(x, na.rm = TRUE)
  }
  
  nc_close(nc)
  
  tibble(
    date  = dates_sel,
    year  = year(dates_sel),
    value = vals
  )
}



# CHL over your grid
chl_ts <- extract_nc_ts_over_grid(
  nc_file  = "./data/raw_data/predictors/Chl/cmems_obs-oc_med_bgc-plankton_my_l4-gapfree-multi-1km_P1D_20130101-20250101.nc",
  var_name = "CHL",
  grid     = grid
)

# SST (in Kelvin; convert to °C if you want: value - 273.15)
sst_ts <- extract_nc_ts_over_grid(
  nc_file  = "./data/raw_data/predictors/SST/SST_MED_SST_L4_NRT_OBSERVATIONS_010_004_c_V2_SST_20130101-20250101.nc",
  var_name = "analysed_sst",
  grid     = grid
) %>%
  mutate(value = value - 273.15)  # now in °C





library(ggplot2)
library(viridis)

library(dplyr)
library(lubridate)

mean_interannual <- function(df) {
  df %>%
    mutate(doy = yday(date)) %>%
    group_by(doy) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value   = sd(value, na.rm = TRUE),
      .groups    = "drop"
    )
}

chl_clim <- mean_interannual(chl_ts)
sst_clim <- mean_interannual(sst_ts)

library(dplyr)
library(lubridate)

mean_interannual <- function(df) {
  df %>%
    mutate(doy = yday(date)) %>%
    group_by(doy) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value   = sd(value, na.rm = TRUE),
      .groups    = "drop"
    )
}

chl_clim <- mean_interannual(chl_ts)
sst_clim <- mean_interannual(sst_ts)


ggplot(chl_clim, aes(x = doy, y = mean_value)) +
  geom_line(color = "darkgreen") +
  geom_ribbon(aes(ymin = mean_value - sd_value,
                  ymax = mean_value + sd_value),
              alpha = 0.2, fill = "darkgreen") +
  labs(
    title = "Mean interannual CHL cycle over grid (2018–2024)",
    x = "Day of year",
    y = "CHL (mg/m³)"
  ) +
  theme_minimal()





ggplot(sst_clim, aes(x = doy, y = mean_value)) +
  geom_line(color = "darkgreen") +
  geom_ribbon(aes(ymin = mean_value - sd_value,
                  ymax = mean_value + sd_value),
              alpha = 0.2, fill = "darkgreen") +
  labs(
    title = "Mean interannual CHL cycle over grid (2018–2024)",
    x = "Day of year",
    y = "CHL (mg/m³)"
  ) +
  theme_minimal()




library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)

plot_one_year_multi_year_lines <- function(df, y_label, title, 
                                           years = NULL,
                                           output_file = NULL) {
  df2 <- df %>%
    # optionally filter years
    { if (!is.null(years)) filter(., year %in% years) else . } %>%
    # remove Feb 29 to make DOY comparable across years
    filter(!(month(date) == 2 & mday(date) == 29)) %>%
    mutate(
      doy       = yday(date),
      clim_date = as.Date("2001-01-01") + doy - 1,  # dummy year
      year      = as.factor(year)
    )
  
  p <- ggplot(df2, aes(x = clim_date, y = value, color = year, group = year)) +
    geom_line(alpha = 0.7, size = 1) +
    geom_point(alpha = 0.4, size = 1) +
    labs(
      title = title,
      x     = "Month",
      y     = y_label,
      color = "Year"
    ) +
    scale_x_date(
      date_labels = "%b",
      date_breaks = "1 month"
    ) +
    scale_color_viridis_d() +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title   = element_text(size = 14, face = "bold"),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  }
  
  p
}

# from previous step:
# sst_ts <- extract_nc_ts_over_grid(...)

# convert to °C if needed:
sst_ts <- sst_ts %>%
  mutate(value = value - 273.15)


sst_oneyear <- plot_one_year_multi_year_lines(
  df       = sst_ts,
  y_label  = "SST (°C)",
  title    = "SST – annual cycle over grid (one line per year)",
  years    = 2018:2023,   # 6 lines; remove this arg to include all years
  output_file = "./figures/sst_oneyear_2018_2023.png"
)

sst_oneyear




chl_oneyear <- plot_one_year_multi_year_lines(
  df       = chl_ts,
  y_label  = "Chlorophyll-a (mg/m³)",
  title    = "CHL – annual cycle over grid (one line per year)",
  years    = 2018:2023,   # again, 6 lines
  output_file = "./figures/chl_oneyear_2018_2023.png"
)

chl_oneyear


#------------------------------------- INTERSEASONAL VARIATIONS OF PREDICTORS -----
