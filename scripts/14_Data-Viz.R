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




#------------- Load functions ------------------------------------------------
source("./utils/Fct_Data-Prep.R")

#--------------------------------- MAP GRID CELLS x sampling x time ------------------------------
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









