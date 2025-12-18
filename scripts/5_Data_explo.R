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

#------------- Load and prep data ----
# Mtdt_7
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")

# Mtdt_7_sel_v1.0
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.0.gpkg")

# predictors_raw_v2.0
pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v2.0.gpkg")

# predictors_raw_v2.2
pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v2.2.gpkg")

# div_indices_v1.0
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))


# predictors_raw_sel_v1.3 ----
pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_sel_v1.3.gpkg")

# remove space and majuscule in habitat names
pred_raw$grouped_main_habitat <- gsub(" ", "_", pred_raw$grouped_main_habitat)
pred_raw$grouped_main_habitat <- tolower(pred_raw$grouped_main_habitat)

# set character as factor
pred_raw$grouped_main_habitat <- as.factor(pred_raw$grouped_main_habitat)

# predictors_raw_v3.1 ----
pred_raw <- sf::st_read("./data/processed_data/predictors/predictors_raw_v3.1.gpkg")
# mtdt_7_sel_v1.1 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")
# mtdt_7 ---
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")

# div_indices_sel_v1.1.gpkg ----
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# Prep ----
# Add coords to mtdt
mtdt <- pred_raw %>%
  dplyr::select(c(x,y,replicates)) %>%
  st_drop_geometry() %>%
  left_join(mtdt, by = "replicates") %>%
  filter(replicates %in% mtdt$replicates)

# Add coordinates to div
div<- pred_raw %>%
  dplyr::select(c(x,y,replicates)) %>%
  st_drop_geometry() %>%
  left_join(div, by = "replicates") %>%
  filter(replicates %in% div$replicates)

# Compute year, month and season from date 
mtdt <- mtdt %>%
  mutate(
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% c(1, 2, 3) ~ "Winter",
      month %in% c(4, 5, 6) ~ "Spring",
      month %in% c(7, 8, 9) ~ "Summer",
      month %in% c(10, 11, 12) ~ "Autumn",
      TRUE ~ NA_character_
    ))



# Replace month number by month name
mtdt <- mtdt %>%
  mutate(
    month = factor(month.abb[month], levels = month.abb)
  )

mtdt$season <- factor(mtdt$season,
                      levels = c("Winter", "Spring", "Summer", "Autumn"))


# Add area_km2 and dist_seabed_depthsampling, bathy_mean, bathy_max to mtdt
mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, area_km2, dist_seabed_depthsampling)
    , by = "replicates"
  )



# Add bathy_mean and bathy_max to mtdt
mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, bathy_mean, bathy_max)
    , by = "replicates"
  )

# Merge div indices to mtdt and predictors
full <- mtdt %>%
  left_join(
    div %>% st_drop_geometry(),
    by = "replicates"
  ) %>%
  left_join(
    pred_raw %>% st_drop_geometry(),
    by = "replicates"
  )



unique(mtdt$month)


# Export 2023 june -----


st_write(
  mtdt %>% filter(year == 2023 & month == "6"),
  "./data/processed_data/Mtdt/mtdt_7_sel_1.1_2023_june.gpkg",
  delete_dsn = TRUE
)
#------------------------------------- Explore Mtd_7 ------------------------------------------------------------------------------------------------------------------------

#--- Check NAs -----
sapply(mtdt, function(x) sum(is.na(x)))






######## Temporal variability #####

#--- Map + hist season | month | year -----

# Cols to plot 
cat_cols <- c("year", "month", "season")

map_categorical_plots(mtdt,
                      cols_to_plot = cat_cols,
                      version = "mtdt_7_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)

map_categorical_plots(mtdt,
                      cols_to_plot = "year",
                      version = "mtdt_7_sel_v1.1_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)

map_categorical_plots(mtdt,
                      cols_to_plot = "season",
                      version = "mtdt_7_sel_v1.1_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)


#--- Map + hist 2023 months ----
map_categorical_plots(
  mtdt %>% filter(year == 2023),
  cols_to_plot = "month",
  version = "mtdt_7_2023_months",
  output_directory = "./figures/Mtdt/Map_Hist", 
  separate_maps = TRUE
)



#--- Count table season | month | year -----

# Cols to count 
cat_cols <- c("year", "month", "season")

export_count_tables(
  mtdt,
  cat_cols = cat_cols,
  output_directory = "./figures/Mtdt/Tables",
  filename_stub = "counts_per_year_month_season",
  version = "mtdt_7_temporal_var",
  layout = "horizontal"  
)

# Count points in spring + summer
mtdt %>%
  filter(season %in% c("Spring", "Summer")) %>%
  nrow() # 698 : 94%






######## Sampling effort  ####






#--- Map + hist of sampling effort -----
cols_to_plot <- c("area_km2", "dist_seabed_depthsampling", "depth_sampling", "bathy_mean", "bathy_max", "PCR_replicates", "field_replicates", "estimated_volume_total")

cols_to_plot <- c("area_km2")

map_index_plots(df = mtdt, 
                version = "mtdt_7_sel_v1.0_buffer_area", 
                output_directory = "./figures/Mtdt/Map_Hist", 
                cols_to_plot = cols_to_plot)



colnames(mtdt)














#--- Map + hist of categorical sampling & seafloor depth -----
# Make bathy categorical 
mtdt <- mtdt %>%
  mutate(
    bathy_max_cat = case_when(
      bathy_max > 0 & bathy_max <= 10 ~ "0-10m",
      bathy_max > 10 & bathy_max <= 20 ~ "10-20m",
      bathy_max > 20 & bathy_max <= 30 ~ "20-30m",
      bathy_max > 30 & bathy_max <= 40 ~ "30-40m",
      bathy_max > 40 & bathy_max <= 50 ~ "40-50m",
      bathy_max > 50 & bathy_max <= 60 ~ "50-60m",
      bathy_max > 60 & bathy_max <= 70 ~ "60-70m",
      bathy_max > 70 & bathy_max <= 80 ~ "70-80m",
      bathy_max > 80 & bathy_max <= 90 ~ "80-90m",
      bathy_max > 80 & bathy_max <= 90 ~ "80-90m",
      bathy_max > 90 & bathy_max <= 100 ~ "90-100m",
      bathy_max > 100 ~ "> 100m",
      
      TRUE ~ NA_character_
    ),
    depth_sampling_cat = case_when(
      depth_sampling > 0 & depth_sampling <= 10 ~ "0-10m",
      depth_sampling > 10 & depth_sampling <= 20 ~ "10-20m",
      depth_sampling > 20 & depth_sampling <= 30 ~ "20-30m",
      depth_sampling > 30 & depth_sampling <= 40 ~ "30-40m",
      depth_sampling > 40 & depth_sampling <= 50 ~ "40-50m",
      depth_sampling > 50 & depth_sampling <= 60 ~ "50-60m",
      depth_sampling > 60 & depth_sampling <= 70 ~ "60-70m",
      depth_sampling > 70 & depth_sampling <= 80 ~ "70-80m",
      depth_sampling > 80 & depth_sampling <= 90 ~ "80-90m",
      depth_sampling > 80 & depth_sampling <= 90 ~ "80-90m",
      depth_sampling > 90 & depth_sampling <= 100 ~ "90-100m",
      depth_sampling > 100 ~ "> 100m",
      TRUE ~ NA_character_
    )
  ) 

# Make factor to ensure correct order in plots
mtdt$bathy_max_cat <- factor(mtdt$bathy_max_cat, levels = c("0-10m", "10-20m", "20-30m", "30-40m", "40-50m", "50-60m", "60-70m", "70-80m", "80-90m", "90-100m", "> 100m"))
mtdt$depth_sampling_cat <- factor(mtdt$depth_sampling_cat, levels = c("0-10m", "10-20m", "20-30m", "30-40m", "40-50m", "50-60m", "60-70m", "70-80m", "80-90m", "90-100m", "> 100m"))



# Map
map_categorical_plots(mtdt,
                      cols_to_plot = c("bathy_max_cat", "depth_sampling_cat"),
                      version = "mtdt_7_bathy_depth_cat",
                      output_directory = "./figures/Predictors/Map_Hist", 
                      separate_maps = TRUE)

# Results : strong spatial bias of depth_sampling : in Occitanie almost only surface sampling and 'deep' sampling > 10m mostly in PACA and Corsica. 






# With sel.v1.1

# Make bathy categorical 
full <- full %>%
  mutate(
    bathy_mean_cat = case_when(
      bathy_mean > 0 & bathy_mean <= 10 ~ "0-10m",
      bathy_mean > 10 & bathy_mean <= 20 ~ "10-20m",
      bathy_mean > 20 & bathy_mean <= 30 ~ "20-30m",
      bathy_mean > 30 & bathy_mean <= 40 ~ "30-40m",
      bathy_mean > 40 & bathy_mean <= 50 ~ "40-50m",
      bathy_mean > 50 & bathy_mean <= 60 ~ "50-60m",
      bathy_mean > 60 & bathy_mean <= 70 ~ "60-70m",
      bathy_mean > 70 & bathy_mean <= 80 ~ "70-80m",
      bathy_mean > 80 & bathy_mean <= 90 ~ "80-90m",
      bathy_mean > 80 & bathy_mean <= 90 ~ "80-90m",
      bathy_mean > 90 & bathy_mean <= 100 ~ "90-100m",
      bathy_mean > 100 ~ "> 100m",
      
      TRUE ~ NA_character_
    )
    )


# Make factor to ensure correct order in plots
full$bathy_mean_cat <- factor(full$bathy_mean_cat, levels = c("0-10m", "10-20m", "20-30m", "30-40m", "40-50m", "50-60m", "60-70m", "70-80m", "80-90m", "90-100m", "> 100m"))



# Map
map_categorical_plots(full,
                      cols_to_plot = "bathy_mean_cat",
                      version = "predictors_raw_sel_v1.3_bathy_depth_MEAN_cat",
                      output_directory = "./figures/Predictors/Map_Hist", 
                      separate_maps = TRUE)

#--- Map + hist of region -----
# Map
map_categorical_plots(mtdt,
                      cols_to_plot = c("region"),
                      version = "mtdt_7_bathy_",
                      output_directory = "./figures/Predictors/Map_Hist", 
                      separate_maps = TRUE)

table(mtdt$region)
# Corse        Occitanie     PACA 
# 299 39%      138 18%       325 43%

#--- Map + hist year -----
map_categorical_plots(mtdt, 
                      cols_to_plot = c("year"),
                      version = "mtdt_7_bathy",
                      output_directory = "./figures/Predictors/Map_Hist", 
                      separate_maps = TRUE)

######## Lockdown ####
#--- Map + hist lockdown -----
map_categorical_plots(mtdt,
                      cols_to_plot = "lockdown",
                      version = "mtdt_7_lockdown",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)









#--- Count season | month | year -----



export_count_tables(
  mtdt,
  cat_cols = "lockdown",
  output_directory = "./figures/Mtdt/Tables",
  filename_stub = "counts_lockdown",
  version = "mtdt_7",
  layout = "horizontal"  
)

# Count points in spring + summer
mtdt %>%
  filter(season %in% c("Spring", "Summer")) %>%
  nrow() # 698 : 94%





######## Methods ####
#--- Map + hist method -----
map_categorical_plots(mtdt,
                      cols_to_plot = "method",
                      version = "mtdt_7_method",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)


map_categorical_plots(mtdt,
                      cols_to_plot = "method",
                      version = "mtdt_7_sel_v1.1_method",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)






#--- Count methods -----



export_count_tables(
  mtdt,
  cat_cols = "method",
  output_directory = "./figures/Mtdt/Tables",
  filename_stub = "counts_method",
  version = "mtdt_7",
  layout = "horizontal"  
)

# Count points in spring + summer
mtdt %>%
  filter(season %in% c("Spring", "Summer")) %>%
  nrow() # 698 : 94%





######## lockdown x season #####

# no lowckdown points + no winter autumn point count
mtdt %>%
  filter(lockdown == "0" & season %in% c("Summer", "Spring")) %>%
  nrow() # 627 : 84%


######## method x season #####

# Heatmap 
mtdt$season <- factor(
  mtdt$season,
  levels = c("Winter", "Spring", "Summer", "Autumn")
)

# Compute counts and percentages
heatmap_data <- mtdt %>%
  count(method, season) %>%
  mutate(percentage = 100 * n / sum(n))

# Plot with black buffer using two text layers
ggplot(heatmap_data, aes(x = method, y = season, fill = n)) +
  geom_tile(color = "black") +
   geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    color = "black", size = 3.5
  ) +
  scale_y_discrete(limits = rev(levels(mtdt$season))) +
  scale_fill_distiller(palette = "RdPu") +
  labs(
    title = "Heatmap of Method vs Season",
    x = "Method",
    y = "Season",
    fill = "Count"
  ) +
  theme(
    plot.background = element_rect(color = "black", fill = "black"),
    panel.background = element_rect(color = "purple", fill = "black"),
    panel.grid = element_blank(),
    axis.title.x = element_text(color = "white"),
    axis.title.y = element_text(color = "white"),
    axis.text.x = element_text(color = "white", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "white"),
    plot.title = element_text(color = "white"),
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(fill = "black", color = NA),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white")
  )

ggsave(
  filename = "./figures/Mtdt/Method_vs_Season_heatmap_mtdt7.jpg",
  width = 8,
  height = 6,
  dpi = 300
)




######## lockdown x season x samplig effort #####

# no lowckdown points + no winter autumn point count
mtdt %>%
  # filter(lockdown == "0" & season %in% c("Summer", "Spring")) %>%
  # filter(field_replicates == 2, PCR_replicates == 24) %>%
  # filter(estimated_volume_total >= 55 & estimated_volume_total <= 65) %>%
  filter(bathy_max <= 50) %>%
  # filter(method == "surface_transect" | method == "seabed_transect") %>%
  nrow() # 121

mtdt %>%
  filter(method == "seabed_transect") %>%
  mutate(estimated_volume_total = as.numeric(estimated_volume_total)) %>%
  with(hist(estimated_volume_total)) # almost all 120 L 

mtdt %>%
  filter(method == "seabed_transect") %>%
  mutate(estimated_volume_total = as.numeric(PCR_replicates)) %>%
  with(hist(PCR_replicates)) # almost all 24
















######## season/month x depth #####
# Scatterplot ----

# bathy_max_cat vs season colored by method
ggplot(mtdt, aes(x = season, y = bathy_max, color = method)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Scatterplot of Bathymetry Max vs Season Colored by Method",
    x = "Season",
    y = "Bathymetry Max (m)",
    color = "Method"
  )

ggsave(
  filename = "./figures/Mtdt/scatterplot_bathy_max_vs_season_method_mtdt7.jpg",
  width = 8,
  height = 6,
  dpi = 300
)

# Heatmap ----

# season x depth_sampling

mtdt$season <- factor(
  mtdt$season,
  levels = c("Winter", "Spring", "Summer", "Autumn")
)

# Compute counts and percentages
heatmap_data <- mtdt %>%
  count(depth_sampling_cat, season) %>%
  mutate(percentage = 100 * n / sum(n))

# Plot with black buffer using two text layers
ggplot(heatmap_data, aes(x = depth_sampling_cat, y = season, fill = n)) +
  geom_tile(color = "black") +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    color = "black", size = 3.5
  ) +
  scale_y_discrete(limits = rev(levels(mtdt$season))) +
  scale_fill_distiller(palette = "RdPu") +
  labs(
    title = "Heatmap of Method vs Season",
    x = "depth_sampling_cat",
    y = "Season",
    fill = "Count"
  ) +
  theme(
    plot.background = element_rect(color = "black", fill = "black"),
    panel.background = element_rect(color = "purple", fill = "black"),
    panel.grid = element_blank(),
    axis.title.x = element_text(color = "white"),
    axis.title.y = element_text(color = "white"),
    axis.text.x = element_text(color = "white", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "white"),
    plot.title = element_text(color = "white"),
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(fill = "black", color = NA),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white")
  )

ggsave(
  filename = "./figures/Mtdt/depth_sampling_cat_vs_Season_heatmap_mtdt7.jpg",
  width = 8,
  height = 6,
  dpi = 300
)








#------------------------------------- Explore raw predictors ---------------------------------------------------------------------------------------------------------------
#--- Data prep ----
# numerical cols
cols <- pred_raw %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-c("y"))%>%
  colnames()


#--- Check NAs, Inf, Negative -----
# numerical cols
cols <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

# NAs
sapply(pred_raw, function(x) sum(is.na(x))) 
# in predictors_raw_v1.2 : 5 NAs for habitat predictors

# Inf
sapply(pred_raw[, cols] %>% st_drop_geometry, function(x) sum(is.infinite(x)))
# No Inf values found

# Negative 
sapply(pred_raw[, cols] %>% st_drop_geometry, function(x) sum(x < 0, na.rm = TRUE))
# 5 in dist_seabed_depthsampling 



#--- Map + hist + summary : numerical data -----

cols_to_plot <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

map_index_plots(df = pred_raw, version = "predictors_raw_num_v1.2", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = cols_to_plot)
map_index_plots(df = pred_raw, version = "predictors_raw_num_v1.3", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = "Boat_density_month")


rm(cols_to_plot)


#--- Map + hist + summary : categorical data -----

cols_to_plot <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.character), where(is.factor)) %>%
  colnames()


map_categorical_plots(df = pred_raw, version = "predictors_raw_cat_v1.2", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = cols_to_plot)



# Map
http://127.0.0.1:18049/graphics/plot_zoom_png?width=1707&height=935
#--- Map + hist env predictors by season -----
# Add season to pred_raw
pred_raw <- pred_raw %>%
  left_join(
    mtdt %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, season) %>%
      distinct(),
    by = "replicates"
  )


# CHL TEMP SAL WS CUR
# Env cols to plot : all column containing "chl", "temp", "sst", "sal", "ws", "vel", "wind"
env_cols <- colnames(pred_raw)[grepl("chl|temp|sst|sal|ws|vel|wind", colnames(pred_raw), ignore.case = TRUE)]
env_cols <- "Boat_density_month"

map_index_plots(df = pred_raw, 
                version = "predictors_raw_num_v1.3_SEASON", 
                output_directory = "./figures/Predictors/Explo_raw/Map_Hist/by_season", 
                cols_to_plot = env_cols,
                factor = "season", 
                panel_text_scale = "s")


rm(env_cols, cols_to_plot)


#--- Summary env predictors by season -----
# 1) Pick your environmental columns
env_cols <- colnames(pred_raw)[grepl("chl|temp|sst|sal|ws|vel|wind", colnames(pred_raw), ignore.case = TRUE)]

# 2) Helper to format numbers (tweak digits if you like)
fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
}

# 3) Build & save a 2x2 grid of tables (one table per season) for a single variable
make_var_figure <- function(df, var, out_dir = "season_stats_figs", digits = 3) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  # Ensure we have the pieces we need
  stopifnot("season" %in% names(df))
  stopifnot(var %in% names(df))
  
  # Keep only the variable + season; coerce to numeric (in case of factors/characters)
  df_use <- df %>%
    dplyr::select(season, !!sym(var)) %>%
    mutate(season = as.factor(season),
           value  = suppressWarnings(as.numeric(!!sym(var))))
  
  # Compute summaries by season
  stats <- df_use %>%
    group_by(season) %>%
    summarize(
      min    = min(value, na.rm = TRUE),
      max    = max(value, na.rm = TRUE),
      mean   = mean(value, na.rm = TRUE),
      sd     = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # If all values are NA for a season, the stats above become Inf/-Inf/NaN; clean that up
  stats <- stats %>%
    mutate(across(-season, ~ ifelse(is.finite(.x), .x, NA_real_)))
  
  # Long format per season for table display
  stats_long <- stats %>%
    pivot_longer(cols = -season, names_to = "Statistic", values_to = "Value") %>%
    mutate(Value = fmt_num(Value, digits = digits))
  
  # Keep seasons in factor order; if you have a specific order, set levels(pred_raw$season) beforehand
  seasons <- levels(df_use$season)
  seasons <- seasons[seasons %in% unique(as.character(df_use$season))]
  
  # Make one small table per season
  grobs <- lapply(seasons, function(s) {
    dat <- stats_long %>%
      filter(season == s) %>%
      dplyr::select(Statistic, Value)
    
    # Nice compact table
    gridExtra::tableGrob(
      dat,
      rows = NULL,
      theme = gridExtra::ttheme_minimal(
        base_size = 10,
        core = list(padding = unit(c(4, 6), "pt")),
        colhead = list(fg_params = list(fontface = "bold"))
      )
    ) %>%
      # Add a subtitle strip with the season name
      gridExtra::arrangeGrob(
        top = grid::textGrob(
          label = as.character(s),
          gp = grid::gpar(fontface = "bold", fontsize = 11)
        )
      )
  })
  
  # If fewer than 4 seasons present, still arrange them; 2x2 grid is default target
  n <- length(grobs)
  ncol <- 2
  nrow <- ceiling(n / ncol)
  
  title <- grid::textGrob(
    label = paste0(var, " — seasonal statistics"),
    gp = grid::gpar(fontface = "bold", fontsize = 14)
  )
  
  lay <- do.call(gridExtra::arrangeGrob, c(grobs, list(ncol = ncol, nrow = nrow, top = title)))
  
  # Save as PNG
  safe_name <- paste0("season_stats_", make.names(var), ".png")
  outfile <- file.path(out_dir, safe_name)
  png(filename = outfile, width = 2000, height = 1500, res = 200)
  grid::grid.draw(lay)
  dev.off()
  
  invisible(outfile)
}

# 4) Run for each environmental variable and export one figure per variable
outputs <- vapply(env_cols, function(v) make_var_figure(pred_raw, v), FUN.VALUE = character(1))

# Optional: print paths to the saved figures
cat("Saved figures:\n", paste0(" - ", outputs), sep = "\n")



#--- Habitat surface proportion -----

habitat_cols <- c("habitats_artificiels_mean",
                  "matte_morte_p_oceanica_mean",
                  "algues_infralittorales_mean", 
                  "soft_bottom_mean",
                  "meadow_mean", 
                  "rock_mean",
                  "coralligenous_mean")

# Sum the habitat surface proportions and check if they equal 1
pred_raw <- pred_raw %>%
  mutate(habitat_sum = rowSums(across(all_of(habitat_cols)), na.rm = TRUE))

summary(pred_raw$habitat_sum)
hist(pred_raw$habitat_sum, breaks = 50)
boxplot(pred_raw$habitat_sum)

# Result : almost all equal to 1. Some are not --> can be due to the removed habitats. 

# Map geom where habitat_sum != 1
t <- pred_raw %>%
  filter(habitat_sum != 1)

ggplot() +
  geom_sf(data = t, aes(color = habitat_sum), size = 2) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Habitat Surface Proportion Sum across Sampling Locations",
       color = "Habitat Sum")

pred_raw <- pred_raw %>%
  dplyr::select(-habitat_sum)
rm(t, habitat_cols)

#--- variance-to-mean ratio ------

sapply(cols, function(h) {
  var_to_mean <- var(pred_raw[[h]]) / mean(pred_raw[[h]])
  var_to_mean
})

# Print cols with variance-to-mean ratio > 1 if NA, ignore it
ratios <- sapply(cols, \(h) var(pred_raw[[h]], na.rm = TRUE) / mean(pred_raw[[h]], na.rm = TRUE))
ratios <- ratios[!is.na(ratios) & ratios > 1]
sort(ratios)

r <- sort(sapply(cols, \(h) var(pred_raw[[h]], na.rm=1)/mean(pred_raw[[h]], na.rm=1)), decreasing=TRUE)
barplot(r[r > 1], las=2)

# make barplot without canyon_shape_area
r_no_canyon <- r[names(r) != "canyon_shape_area"]
par(mar=c(10,4,2,1))
barplot(r_no_canyon[r_no_canyon > 1], las=2,  cex.names = 0.7)

# make barplot without canyon_shape_area nor canyon_shape_length
r_no_canyon_length <- r[names(r) != "canyon_shape_area" & names(r) != "canyon_shape_length"]
par(mar=c(10,4,2,1))
barplot(r_no_canyon_length[r_no_canyon_length > 1], las=2, cex.names = 0.7)



#--- canyon_shape_area has very high variance-to-mean ratio ---
hist(pred_raw$canyon_shape_area, breaks=50)

# map pred_raw with geom fill color by canyon_shape_area
# Add canyon to the map
canyon <- st_read("./data/raw_data/predictors/Canyon/canyon_med.geojson")
plot(canyon)

# Project crs canyon to pred_raw crs
canyon <- st_transform(canyon, st_crs(pred_raw))

# Crop canyon to pred_raw extent
canyon <- st_crop(canyon, st_bbox(pred_raw))


ggplot() +
  geom_sf(data = canyon, fill = "lightblue", color = "blue", alpha = 0.5) +
  geom_sf(data = pred_raw, aes(color = canyon_shape_area), size = 2) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Canyon Shape Area across Sampling Locations",
       color = "Canyon Shape Area")


# map canyon with fill color by arear_km2
ggplot() +
  geom_sf(data = canyon, aes(fill = area_km2), color = "blue", alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Canyon Area (km²)",
       fill = "Area (km²)")

# Conclusion : huge canyon off the cast of Occitanie. 



# Interpreting the variance-to-mean ratio : if the ratio is > 1, the variance is greater than the mean, indicating a high level of dispersion in the data.

#--- boxplot ----
plot_list <- lapply(cols, function(cn) {
  ggplot(pred_raw, aes(y = .data[[cn]])) +
    geom_boxplot(outlier.alpha = 0.3, fill = "steelblue") +
    theme_minimal(base_size = 9) +
    labs(title = cn, y = NULL, x = NULL) +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
})

# Split into groups of 10
plot_groups <- split(plot_list, ceiling(seq_along(plot_list) / 10))

# Save each as a 2x5 panel
for (i in seq_along(plot_groups)) {
  panel <- wrap_plots(plot_groups[[i]], ncol = 5, nrow = 2)
  
  ggsave(
    filename = sprintf("boxplot_panel_%02d.png", i),
    plot = panel,
    width = 12, height = 6, dpi = 300
  )
}

# Results, data with strong outliers : 
# dist_shore
# canyon data
# roughness_min
# slope_min
# tri_min
# sal_max_1y
# chl_week_mean
# chl_week_max




#--- 1.5xIQR rule -----
# Identify variables with outliers using 1.5xIQR rule
# see here for method : https://www.khanacademy.org/math/statistics-probability/summarizing-quantitative-data/box-whisker-plots/a/identifying-outliers-iqr-rule
outlier_vars <- names(Filter(find_outliers, pred_raw[cols] %>% st_drop_geometry()))  # Apply only to numeric columns

# Print variable names that contain outliers
if (length(outlier_vars) > 0) {
  print("Variables with outliers:")
  print(outlier_vars)
  print("Percentage of numerical variables with outliers:")
  print(paste(round(length(outlier_vars)/length(cols)*100,2), "%"))
  print("Number of numerical variables with outliers:")
  print(length(outlier_vars))
} else {
  print("No variables with outliers found.")
}

#--- Max >= 10*median rule --------
# Identify high variability predictors : if max >= to 10*median 

for (col in names(pred_raw)) {
  if (is.numeric(pred_raw[[col]])) {
    max_val <- max(pred_raw[[col]], na.rm = TRUE)
    med_val <- median(pred_raw[[col]], na.rm = TRUE)
    threshold <- 10 * med_val

    if (max_val >= threshold) {
      message(col)
    } 
  }
}



#--- Near zero var -----
# Check near 0 variance variable 
print(caret::nearZeroVar(pred)) # output designate columns 28 and 31
summary(pred[, c(28, 31)]) # they correspond to WIND_w_min and WIND_m_min

caret::nearZeroVar(pred_raw %>% st_drop_geometry(), saveMetrics = TRUE)
# freqRatio 
# the ratio of frequencies for the most common value over the second most common value --> the closest to 1, the more balanced the values are
# percentUnique
# the percentage of unique data points out of the total number of data points
# zeroVar
# a vector of logicals for whether the predictor has only one distinct value
# nzv
# a vector of logicals for whether the predictor is a near zero variance predictor

# Results : 
# nvz : canyon_type, habitats_artificiels_mean, rock_mean + win_min 7day 1 month and 1year
# zeroVar :  win_min 7day 1 month and 1year

# Env predictors : more variability in mean data compared to min and max. 

#--- Normality ----
#--- Histograms + Normality test ---------------------
# Function to create histogram with Shapiro-Wilk normality test
create_hist_plot <- function(var_name, values) {
  
  # Remove NAs to avoid errors
  values <- na.omit(values)
  
  # Perform Shapiro-Wilk Normality Test
  shapiro_test <- shapiro.test(values)
  p_value <- shapiro_test$p.value
  interpretation <- ifelse(p_value >= 0.05, "Normal", "Not Normal")
  
  # Create histogram
  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_histogram(bins = 50, fill = "gray", color = "black", alpha = 0.6) +
    labs(title = var_name, x = "", y = "Frequency") +
    theme_minimal(base_size = 12) +
    annotate("text", x = max(values, na.rm = TRUE) * 0.7, 
             y = max(hist(values, plot = FALSE)$counts) * 0.9,
             label = paste0("SW p = ", round(p_value, 3), "\n(", interpretation, ")"),
             size = 4, hjust = 0, vjust = 1, color = "black", fontface = "bold",
             fill = "white", label.size = 0.5)
  
  return(p)
}


# Generate histograms for each variable 
ind <- pred_raw[, cols] %>%
  st_drop_geometry()

plots <- lapply(names(ind), function(var) create_hist_plot(var, ind[[var]]))

# Save all plots
for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("./figures/Predictors/Hist_Normality-test/predictors_raw_v1.2_hist_norm_", names(ind)[i], ".png"),
    plot = plots[[i]],
    width = 8, height = 6, dpi = 300
  )
}


#--- Table summary -----
# Var-to-mean ratio + IQG + median rule + NZV + Normality 
#  Packages 
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(sf)
  library(caret)
  library(readr)   # for write_csv (optional)
})

#  Helpers 
has_outliers_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA)
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr <- qs[2] - qs[1]
  if (iqr == 0) return(FALSE)
  lower <- qs[1] - 1.5 * iqr
  upper <- qs[2] + 1.5 * iqr
  any(x < lower | x > upper)
}

safe_shapiro_p <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  if (n > 5000) {
    set.seed(123)
    idx <- base::sample.int(n, 5000, replace = FALSE)
    x <- x[idx]
  }
  out <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
  as.numeric(out)
}

fmt_flag <- function(x) ifelse(is.na(x), NA, ifelse(x, "Yes", "No"))

#  Data selection (sf-safe) 
num_df <- pred_raw %>%
  sf::st_drop_geometry() %>%
  dplyr::select(where(is.numeric))

vars <- names(num_df)

#  NZV metrics from caret 
nzv_metrics <- caret::nearZeroVar(num_df, saveMetrics = TRUE) %>%
  tibble::rownames_to_column("Variable") %>%
  tibble::as_tibble() %>%
  dplyr::select(Variable, freqRatio, percentUnique, zeroVar, nzv)

#  Core diagnostics per variable 
diag_tbl <- purrr::map_dfr(vars, function(v) {
  x <- num_df[[v]]
  m   <- mean(x, na.rm = TRUE)
  vr  <- stats::var(x,  na.rm = TRUE)
  vmr <- ifelse(is.finite(m) && m != 0, vr / m, NA_real_)
  
  out_iqr <- has_outliers_iqr(x)
  
  mx   <- suppressWarnings(max(x, na.rm = TRUE))
  med  <- stats::median(x, na.rm = TRUE)
  max_ge_10x_med <- ifelse(is.finite(mx) && is.finite(med), mx >= 10 * med, NA)
  
  p_sw <- safe_shapiro_p(x)
  normality <- dplyr::case_when(
    is.na(p_sw)   ~ NA_character_,
    p_sw >= 0.05  ~ "Normal",
    TRUE          ~ "Not Normal"
  )
  
  tibble::tibble(
    Variable = v,
    Variance_to_Mean = vmr,
    Outliers_IQR = out_iqr,
    Max_ge_10x_Median = max_ge_10x_med,
    Shapiro_p = p_sw,
    Normality = normality
  )
})

#  Join NZV + finalize summary 
summary_tbl <- diag_tbl %>%
  dplyr::left_join(nzv_metrics, by = "Variable") %>%
  dplyr::mutate(
    NearZeroVar = dplyr::coalesce(nzv, FALSE),
    ZeroVar     = dplyr::coalesce(zeroVar, FALSE),
    `Variance / Mean Ratio`        = Variance_to_Mean,
    `Outlier (1.5×IQR)`            = fmt_flag(Outliers_IQR),
    `Max ≥ 10×Median`              = fmt_flag(Max_ge_10x_Median),
    `Near-Zero Variance`           = fmt_flag(NearZeroVar),
    `Zero Variance`                = fmt_flag(ZeroVar),
    `Normality (Shapiro-Wilk p)`   = Shapiro_p,
    `Normality Interpretation`     = Normality,
    `Key Conclusion` = purrr::pmap_chr(
      list(Variance_to_Mean, Outliers_IQR, Max_ge_10x_Median, NearZeroVar, ZeroVar, Normality),
      function(vmr, oi, max10, nzv, zv, normi) {
        con <- character(0)
        if (isTRUE(zv))  con <- c(con, "Constant (drop)")
        else if (isTRUE(nzv)) con <- c(con, "Near-zero variance (consider drop)")
        if (!is.na(vmr) && vmr > 1) con <- c(con, "High dispersion")
        if (isTRUE(oi))             con <- c(con, "IQR outliers")
        if (isTRUE(max10))          con <- c(con, "Heavy right tail")
        if (!is.na(normi) && normi == "Not Normal") con <- c(con, "Non-normal (transform?)")
        if (!length(con)) "No red flags" else paste(unique(con), collapse = "; ")
      }
    )
  ) %>%
  dplyr::select(
    Variable,
    `Variance / Mean Ratio`,
    `Outlier (1.5×IQR)`,
    `Max ≥ 10×Median`,
    `Near-Zero Variance`,
    `Zero Variance`,
    `Normality (Shapiro-Wilk p)`,
    `Normality Interpretation`,
    freqRatio,
    percentUnique,
    `Key Conclusion`
  ) %>%
  dplyr::arrange(dplyr::desc(`Variance / Mean Ratio`))

# Inspect / export
print(summary_tbl, n = Inf)
 readr::write_csv(summary_tbl, "./output/predictors_raw_v1.2_variance_outliers_norm.csv")

#--- BOATS (before extraction) ----

# Explore ncdf ----
library(ncdf4)
library(lubridate)

nc <- nc_open("./data/raw_data/predictors/Boats/boats_med_2018_2024_bymonth/boats_metropole_100m.nc")

lon   <- ncvar_get(nc, "lon")
lat   <- ncvar_get(nc, "lat")
boats <- ncvar_get(nc, "boats")  # [lon, lat, yearmonth]

d <- dim(boats)
n_months <- d[3]

# Months from 2018-01 to 2024-12 (84 months)
ym <- format(seq(as.Date("2018-01-01"), by = "1 month", length.out = n_months),
             "%Y%m")

length(ym)   # should be 84
head(ym)
tail(ym)

# boats[lon, lat, month]
summary_ym <- t(apply(boats, 3, function(x) {
  x <- as.numeric(x)
  c(
    n    = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    sd   = sd(x,  na.rm = TRUE)
  )
}))

summary_ym <- as.data.frame(summary_ym)
summary_ym$yearmonth <- ym

summary_ym

# Plot mean boats per month ----
library(dplyr)
library(ggplot2)

summary_ym <- summary_ym %>%
  mutate(
    yearmonth = as.character(yearmonth),
    date = as.Date(paste0(yearmonth, "01"), format = "%Y%m%d")
  )

ggplot(summary_ym, aes(x = date, y = mean)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Mean Boats Detection Over Time (Corse)",
    x = "Year-Month",
    y = "Mean Boats Detected"
  ) +
  theme_minimal()



library(lubridate)

n_months <- dim(boats)[3]
ym <- format(seq(ymd("2018-01-01"), by = "1 month", length.out = n_months), "%Y%m")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

grid_idx <- expand.grid(
  lon = lon,
  lat = lat
)
nrow(grid_idx)   # 1400 * 2747 = 3,845,800 cells




# Map for a specific month ----
library(ggplot2)
library(sf)

month_i <- which(ym == "201807")   # July 2018

df <- grid_idx
df$boats <- as.vector(boats[,,month_i])
df <- df[df$boats > 0, ]   # optional: remove zero/no-detection

df <- df %>% tidyr::drop_na(lon, lat)

sf_month <- st_as_sf(df, coords = c("lon", "lat"), crs = 3857)

ggplot(sf_month) +
  geom_sf(aes(color = boats), size = 0.1) +
  scale_color_viridis_c() +
  labs(
    title = paste("Boats – Corse –", ym[month_i]),
    color = "boats"
  ) +
  theme_minimal()



# Loop through all months and save maps ----
out_dir <- "./figures/Predictors/boats/boats_fr_monthly/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(ym)) {
  message("Mapping month: ", ym[i])
  
  # Base grid + values for this month
  df <- grid_idx
  df$boats <- as.vector(boats[,,i])
  
  # Remove rows with NA coords or NA boats
  df <- df %>%
    filter(!is.na(lon), !is.na(lat), !is.na(boats))
  
  # Optional: remove zero detections for performance
  df <- df %>% filter(boats > 0)
  
  if (nrow(df) == 0) {
    message("  No non-zero boats for ", ym[i], " — skipping.")
    next
  }
  
  sf_month <- st_as_sf(df, coords = c("lon", "lat"), crs = 3857)
  
  g <- ggplot(sf_month) +
    geom_sf(aes(color = boats), size = 0.15) +
    scale_color_viridis_c() +
    labs(
      title = paste("Boats – Corse –", ym[i]),
      color = "boats"
    ) +
    theme_minimal()
  
  ggsave(
    file.path(out_dir, paste0("boats_", ym[i], ".png")),
    g,
    width = 7, height = 5, dpi = 300
  )
}






#--- BOATS (after extraction) ----
boats <- st_read("./data/processed_data/predictors/Boats/mtdt_7_boats_month_ALL.gpkg")
plot(boats)
sf::sf_use_s2(FALSE)
coords <- st_coordinates(st_centroid(boats))

boats <- boats %>%
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  )


# Map + hist + summary : boats -----


# all together
map_index_plots(df = boats, version = "with_mtdt7", output_directory = "./figures/Predictors/Explo_raw/Map_Hist", cols_to_plot = c("Boat_density_month"))


# by season | month | year | yearmonth
boats$date <- as.Date(boats$date)
boats <- boats %>%
  mutate(
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% c(1, 2, 3) ~ "Winter",
      month %in% c(4, 5, 6) ~ "Spring",
      month %in% c(7, 8, 9) ~ "Summer",
      month %in% c(10, 11, 12) ~ "Autumn",
      TRUE ~ NA_character_),
    yearmonth = format(date, "%Y%m"))



map_index_plots(df = boats, 
                version = "with_mtdt7_by_month", 
                output_directory = "./figures/Predictors/Explo_raw/Map_Hist/by_season", 
                cols_to_plot = c("Boat_density_month"),
                factor = "month", 
                panel_text_scale = "s")



#--- Cleanup ----
rm(diag_tbl, fmt_flag, has_outliers_iqr, make_var_figure, num_df, safe_shapiro_p, summary_tbl, nzv_metrics, vars)
rm(ad_res, dag_res, ind, ks_res, normality_results, plots, shapiro_res, col, cols, i, max_val, med_val, median_val, normality_status, outlier_vars, threshold)

dev.off()
gc()
#------------------------------------- Explore div indices --------------
#--- Data prep ----
indicators <- div







#--- Check NAs ----
sapply(indicators, function(x) sum(is.na(x))) 

#--- Hist + summary ----

## Make a numeric copy of indicators (except 'replicates')
indicators_num <- indicators %>%
  mutate(across(-replicates, ~ as.numeric(.)))



## List of columns to plot
cols_to_plot <- names(indicators_num)[names(indicators_num) != "replicates"]
cols_to_plot <- c("R", "Crypto", "Elasmo")

## Compute histogram plots safely
plots_list <- lapply(cols_to_plot, function(col_name) {
  gg_hist_summary(indicators_num[[col_name]], col_name = col_name)
})

## Name the list for easier reference
names(plots_list) <- cols_to_plot

# Display plots
plots_list


























# Save plots 

# version string (set manually)
version <- "div_indices_v1.0"  # <-- change this as needed

# create output folder if needed
output_dir <- "./figures/Div_indices/Hist/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# loop through each plot in your list
for (col_name in names(plots_list)) {
  file_name <- paste0("Hist-summary_", col_name, "_", version, ".jpg")
  file_path <- file.path(output_dir, file_name)
  
  # save the plot
  ggsave(
    filename = file_path,
    plot = plots_list[[col_name]],
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message("Saved: ", file_path)
}




rm(plots_list, indicators_num, cols_to_plot, col_name, file_name, file_path, output_dir, version)


#--- Map + hist + summary -----

map_index_plots(df = indicators, version = "div_indices_sel_v1.1", output_directory = "./figures/Div_indices/Map_Hist", cols_to_plot = colnames(indicators[4:13]))


























# Div x sampling effort | region | season-month | depth | habitat | etc. [TODO] ---


#--- R x sampling effort [to do] -----
# Make a panel of 4 scatter plots : 
# R x area_km2
# R x estimated_volume_total
# R x PCR_replicates
# R x depth_sampling







#--- R x bathy [to do] -----


#-------------------------------------  Link R - variables explo -----------------------------------------
#--- Load data ----
# predictors_sel_v.1.5 
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.5.gpkg")

# remove space and majuscule in habitat names
pred$grouped_main_habitat <- gsub(" ", "_", pred$grouped_main_habitat)
pred$grouped_main_habitat <- tolower(pred$grouped_main_habitat)

# set character as factor
pred$grouped_main_habitat <- as.factor(pred$grouped_main_habitat)

# check levels
levels(pred$grouped_main_habitat)
table(pred$grouped_main_habitat)

unique(pred$grouped_main_habitat)

# mtdt_7_sel_v1.1 
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

# div_indices_sel_v1.1.gpkg 
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# tot 
tot <- pred %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(mtdt), by = "replicates") %>%
  left_join(div, by = "replicates")

# export tot 
 st_write(st_as_sf(tot), "./data/processed_data/tot_v1.0.gpkg", delete_dsn = TRUE)


#--- R ~ volume/area/depth_sampling  ----------------------------
predictor.variable.names = c("estimated_volume_total",
                             "depth_sampling", "area_km2")

spatialRF::plot_training_df(
  data = tot,
  dependent.variable.name = "R",
  predictor.variable.names = predictor.variable.names,
  ncol = 3,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

# Save plot
ggsave(
  filename = "./figures/Div_indices/R_vs_volume_area_depth-sampling.png",
  width = 12,
  height = 4,
  dpi = 300
)


# Correlation test
cor.test(tot$R, tot$estimated_volume_total, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  tot$R and tot$estimated_volume_total
# t = 11.261, df = 635, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3411226 0.4707572
# sample estimates:
#       cor 
# 0.4079942 


cor.test(tot$R, tot$depth_sampling, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  tot$R and tot$depth_sampling
# t = -1.8355, df = 635, p-value = 0.0669
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.149487695  0.005063804
# sample estimates:
#         cor 
# -0.07264805 


cor.test(tot$R, tot$area_km2, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  tot$R and tot$area_km2
# t = -0.7347, df = 635, p-value = 0.4628
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.10658545  0.04864995
# sample estimates:
#         cor 
# -0.02914347 




#--- R ~ method/project ----------------------------


library(ggplot2)
library(forcats)


# Remove AUV and glider methods
dt <- tot %>%
  filter(!method %in% c("AUV", "glider"))

# By method 
p <- ggplot(dt, aes(x = method, y = R)) +
  geom_violin(fill = "grey90", color = "grey60", width = 0.9, alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "black") +
  labs(
    title = "R by sampling method",
    x = "Sampling method",
    y = "R"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
p

# Save plot 
ggsave(
  filename = "./figures/Div_indices/R_by_method.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


# By project 
p <- ggplot(dt, aes(x = project, y = R)) +
  geom_violin(fill = "grey90", color = "grey60", width = 0.9, alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "black") +
  labs(
    title = "R by sampling method",
    x = "Sampling method",
    y = "R"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
p

# Save plot 
ggsave(
  filename = "./figures/Div_indices/R_by_project.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


#--- R ~ habitat ----------------------------

# Main habitat
p <- ggplot(tot, aes(x = grouped_main_habitat, y = R)) +
  geom_violin(fill = "grey90", color = "grey60", width = 0.9, alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "black") +
  labs(
    title = "R by sampling method",
    x = "Sampling method",
    y = "R"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
p

# Save plot 
ggsave(
  filename = "./figures/Div_indices/R_by_habitat.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)




# Nb habitat / km2
spatialRF::plot_training_df(
  data = tot,
  dependent.variable.name = "R",
  predictor.variable.names = "nb_habitat_per_km2",
  ncol = 2,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

# Save plot
ggsave(
  filename = "./figures/Div_indices/R_vs_nb-habitat_per_km2.png",
  width = 6,
  height = 4,
  dpi = 300
)

#--- R ~ regions ----------------------------
p <- ggplot(tot, aes(x = region, y = R)) +
  geom_violin(fill = "grey90", color = "grey60", width = 0.9, alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "black") +
  labs(
    title = "R by region",
    x = "Region",
    y = "R"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
p

# Save plot
ggsave(
  filename = "./figures/Div_indices/R_by_region.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


#--- R ~ year ----------------------------
tot$year <- year(tot$date) %>% as.factor()


p <- ggplot(tot, aes(x = year, y = R)) +
  geom_violin(fill = "grey90", color = "grey60", width = 0.9, alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "steelblue", alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "point", size = 2.2, color = "black") +
  labs(
    title = "R by year",
    x = "year",
    y = "R"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
p

# Save plot
ggsave(
  filename = "./figures/Div_indices/R_by_year.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)




#--- Characteristics of R hotspots ----------------------------
# Hist of R --
ggplot(tot, aes(x = R)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of R values",
    x = "R",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# Boxplot plot of R --
boxplot(tot$R)


# Map of hotspots --
# Add geometry
tot <- pred %>%
  dplyr::select(replicates, geom) %>%
  left_join(st_drop_geometry(tot), by = "replicates")

# Hotspots definition : R >  60
hot <- tot %>% filter(R > 60)

# Map hotspots

ggplot() +
  # non-hotspots: transparent fill, visible outline
  geom_sf(data = tot, 
          fill = "white",          # no fill
          color = "black",    # outline color
          size = 0.1,
          alpha = 0.5) +      # translucent outline
  
  # hotspots: colored by R
  geom_sf(data = hot, 
          aes(color = R), 
          size = 2.2) +
  
  scale_color_viridis_c() +
  
  theme_minimal() +
  labs(title = "R Hotspots (R ≥ 60)",
       color = "R")


# Save plot
ggsave(
  filename = "./figures/Div_indices/R_hotspots_map.png",
  width = 8,
  height = 6,
  dpi = 300
)


# Characteristics of hotspots ---
nrow(hot) # 45 sites

# Hist
ggplot(hot, aes(x = R)) +
  geom_histogram(bins = 20, fill = "orange", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of R values in Hotspots",
    x = "R",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
# Boxplot
boxplot(hot$R)


# Methods used in hotspots
table(hot$method)
# dive_motionless                     dive_transect motionless_descended_from_surface                   seabed_transect                  surface_transect 
# 9                                 3                                 2                                10                                21 

table(hot$estimated_volume_total) 
# 60  64  90 120 
# 24   1   1  19 

summary(hot$estimated_volume_total)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 60.00   60.00   60.00   86.09  120.00  120.00 

boxplot(hot$estimated_volume_total)


# !! Hotspot --> no volume < 60

hot$year <- year(hot$date)
table(hot$year)
# 2018 2019 2020 2021 2022 2023 2024 
# 4    3    1   13    2   13    9 

table(hot$region)
# Corse Occitanie      PACA 
# 14        14        17 

tot %>% filter(estimated_volume_total < 60) %>% pull(R) %>% summary

#--- Max R points ----------------------------
# Extract point with max R 
max_r_point <- tot %>% filter(R == max(R, na.rm = TRUE))


#------------------------------------ Occurences data explo ------------------------------
#--- Load data ----
occ <- readr::read_delim("./data/processed_data/eDNA/occ_sel-on-mtdt_3.csv")
traits <- read.csv("./data/processed_data/Traits/species_traits_NA-resolved_v1.0.csv", sep = ";", header=T)85

## Nb species ----
dim(occ) # 262 spp 
# 258 spp occ our sel 
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))


## Nb de familles ----
names(traits)
length(unique(traits$Family))
sort(unique(traits$Family))

## Nb of elasmo ----
traits %>% filter(elasmo == 1) %>% dim # 25 (à vérif avec Célia) --> 23 manuscript Celia
# 23 elasmobranch species (10 sharks and 11 rays and 2 skates) were identified. 
# Notable examples include the critically endangered bull ray (Aetomylaeus bovinus) and angelshark (Squatina squatina).


t <- traits %>% filter(elasmo == 1)
t$species
# [1] "Aetomylaeus bovinus"                        "Alopias superciliosus"                     
# [3] "Alopias vulpinus"                           "Cetorhinus maximus"                        
# [5] "Dasyatis pastinaca"                         "Dasyatis tortonesei"                       
# [7] "Dipturus oxyrinchus"                        "Etmopterus spinax"                         
# [9] "Galeus melastomus"                          "Isurus oxyrinchus"                         
# [11] "Mobula mobular"                             "Mustelus mustelus"                         
# [13] "Myliobatis aquila"                          "Prionace glauca"                           
# [15] "Pteroplatytrygon violacea"                  "Raja brachyura"                            
# [17] "Raja miraletus"                             "Raja undulata"                             
# [19] "Rostroraja alba"                            "Scyliorhinus canicula"                     
# [21] "Scyliorhinus stellaris"                     "Squalus blainville"                        
# [23] "Squatina squatina"                          "Torpedo marmorata"                         
# [25] "Raja asterias.Raja clavata.Raja polystigma"

tt <- traits %>% filter(elasmo == 1) %>% filter(IUCN_combined == "Threatened" | IUCN_combined == "CR" | IUCN_combined == "EN" | IUCN_combined == "VU")
nrow(tt) # 15
tt$species 
# [1] "Aetomylaeus bovinus"    "Alopias superciliosus"  "Alopias vulpinus"      
# [4] "Cetorhinus maximus"     "Dasyatis pastinaca"     "Etmopterus spinax"     
# [7] "Isurus oxyrinchus"      "Mobula mobular"         "Mustelus mustelus"     
# [10] "Myliobatis aquila"      "Raja undulata"          "Rostroraja alba"       
# [13] "Scyliorhinus stellaris" "Squatina squatina"      "Torpedo marmorata"  






##  Nb of RedList + examples ----
traits %>% filter(IUCN_combined == "Threatened" | IUCN_combined == "CR" | IUCN_combined == "EN" | IUCN_combined == "VU") %>% dim # 23 
t <- traits %>% filter(IUCN_combined == "Threatened" | IUCN_combined == "CR" | IUCN_combined == "EN" | IUCN_combined == "VU")
t$species
# [1] "Aetomylaeus bovinus"     "Alopias superciliosus"   "Alopias vulpinus"       
# [4] "Anguilla anguilla"       "Balistes capriscus"      "Cetorhinus maximus"     
# [7] "Dasyatis pastinaca"      "Epinephelus costae"      "Epinephelus marginatus" 
# [10] "Etmopterus spinax"       "Istiophorus platypterus" "Isurus oxyrinchus"      
# [13] "Mobula mobular"          "Mola mola"               "Mustelus mustelus"      
# [16] "Myliobatis aquila"       "Pomatomus saltatrix"     "Raja undulata"          
# [19] "Rostroraja alba"         "Schedophilus ovalis"     "Scyliorhinus stellaris" 
# [22] "Squatina squatina"       "Torpedo marmorata"    





# Nb of cryptiques species + examples  ----



## Richesse ----
summary(div$R)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   22.00   36.00   35.76   49.00   85.00 

# 35 espèces en moyenne par site +/- 17
# standard deviation 
sd(div$R)


# "The overall diversity map (Figure 2) highlights three major biodiversity hotspots with comparable species richness (90 to 108 species). The northwestern region of Corsica near Calvi, the Cerbère-Banyuls Marine Nature Reserve, and the Calanques National Park and its surrounding waters near Marseille all support up to 100 species," --> correct with up to 85 species



## Elasmo ----

summary(div$Elasmo)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   1.262   2.000   8.000 
hist(div$Elasmo)

# Count how many sites have at least one elasmo
sum(div$Elasmo >= 1) # 389 sites with at least one elasmo
389/637 * 100 # 61.067

# Map of elsamo > 0
div_sf <- pred %>%
  dplyr::select(replicates, geom) %>%
  left_join(st_drop_geometry(div), by = "replicates")
ggplot() + 
  geom_sf(data = div_sf, aes(color = Elasmo), size = 0.5) +
  scale_color_viridis_c() +
  labs(title = "Number of Elasmobranch Species Detected per Site",
       color = "Elasmobranchs") +
  theme_minimal()



# " These detections represent approximately 60% of the Teleostei and 38% of the Elasmobranchii listed in our regional checklist for the northwestern Mediterranean. "











#-------------------------------------------- INFOGRAPHICS --------------------------------
library(tidyverse)

#--- Load data ----
occ    <- readr::read_delim("./data/processed_data/eDNA/occ_sel-on-mtdt_3.csv")
traits <- read.csv("./data/processed_data/Traits/species_traits_used_for_div_indices_v1.0.csv",
                   sep = ";", header = TRUE)
div    <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))
mtdt_7 <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")
names(mtdt_7)
sum(mtdt_7$field_replicates)


#--- Key numbers ----
n_spp            <- occ %>% ncol()              # total species
n_elasmo         <- traits %>% filter(elasmo == 1) %>% nrow() 
n_threat         <- traits %>% filter(IUCN_combined %in% c("Threatened","CR","EN","VU")) %>%
  nrow()
n_threat_elasmo  <- traits %>% filter(elasmo == 1,
                                      IUCN_combined %in% c("Threatened","CR","EN","VU")) %>%
  nrow()

rich_mean <- mean(div$R, na.rm = TRUE)
rich_sd   <- sd(div$R,   na.rm = TRUE)
rich_min  <- min(div$R,  na.rm = TRUE)
rich_max  <- max(div$R,  na.rm = TRUE)

# optional: fixed numbers you already checked in the manuscript
n_sharks <- 10
n_rays   <- 11
n_skates <- 2

#--- Data frame for the infographic ----
info_df <- tibble(
  block = factor(c("Total species", "Elasmobranchs",
                   "Threatened species", "Richness per site"),
                 levels = c("Total species", "Elasmobranchs",
                            "Threatened species", "Richness per site")),
  big   = c(
    n_spp,
    n_elasmo,
    n_threat,
    sprintf("%.1f ± %.1f", rich_mean, rich_sd)
  ),
  small = c(
    "",
    sprintf("%d sharks, %d rays, %d skates", n_sharks, n_rays, n_skates),
    sprintf("%d threatened elasmobranchs", n_threat_elasmo),
    sprintf("min = %d, max = %d", rich_min, rich_max)
  )
)

#--- Simple infographic panel ----
ggplot(info_df, aes(x = 0, y = 0)) +
  geom_tile(width = 0.95, height = 0.95, alpha = 0.05) +
  geom_text(aes(label = block), y = 0.35, fontface = "bold", size = 5) +
  geom_text(aes(label = big),   y = 0.05, fontface = "bold", size = 7) +
  geom_text(aes(label = small), y = -0.25, size = 4) +
  facet_wrap(~ block, ncol = 2) +
  coord_fixed() +
  theme_void() +
  theme(strip.text = element_blank())

ggsave(
  filename = "./figures/infographics_v1.0.png",
  width = 8,
  height = 6,
  dpi = 300
)










#--------------------- waffle infographic -----------------------
llibrary(waffle)

# PARAMETERS
n_sq <- 10                 # number of squares (default = 1)
sq_col <- "deepskyblue"   # color of the square

# Build waffle grid (force a single tiny square)
df <- data.frame(values = n_sq)

ggplot(df, aes(values = values, fill = "sq")) +
  geom_waffle(
    n_rows = 1,
    size = 0,
    colour = NA,
    flip = FALSE
  ) +
  scale_fill_manual(values = c("sq" = sq_col)) +
  coord_equal() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.position  = "none",
    plot.margin      = margin(40, 40, 40, 40)   # ensures the square stays small inside a large black frame
  )



waffle(
  parts = c(85, 36, 1),
  rows = 3,
  size = 1,
  colors = c("deepskyblue", "deepskyblue", "deepskyblue")
) +
  theme_void() +
  theme(
    legend.position = "none")







library(waffle)
library(ggplot2)
library(patchwork)   # to stack plots

p1 <- waffle(
  parts  = c(one = 1),
  rows   = 1,
  size   = 1,
  colors = "deepskyblue"
) +
  theme_void() +
  theme(legend.position = "none")

p2 <- waffle(
  parts  = c(thirtysix = 36),
  rows   = 1,
  size   = 1,
  colors = "deepskyblue"
) +
  theme_void() +
  theme(legend.position = "none")

p3 <- waffle(
  parts  = c(eightyfive = 85),
  rows   = 1,
  size   = 1,
  colors = "deepskyblue"
) +
  theme_void() +
  theme(legend.position = "none")

(p1 / p2 / p3) +
  plot_annotation(theme = theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA)
  ))



library(waffle)
library(ggplot2)
library(patchwork)

make_line <- function(n_fill, total = 85, col = "deepskyblue") {
  waffle(
    parts  = c(fill = n_fill, empty = total - n_fill),
    rows   = 1,
    size   = 1,
    colors = c("fill" = col, "empty" = "white")
  ) +
    theme_void() +
    # theme(
    #   legend.position  = "none",
    #   plot.background  = element_rect(fill = "black", color = NA),
    #   panel.background = element_rect(fill = "black", color = NA)
    # )
}

p1 <- make_line(1)
p2 <- make_line(36)
p3 <- make_line(85)

p1 / p2 / p3   # three lines: 1, 36, 85 squares, all same size



# save as png
ggsave(
  filename = "./figures/infographic_waffle_v1.0.png",
  plot = p1 / p2 / p3 ,
  width = 4,
  height = 6,
  dpi = 300
)




















## 1) Min / mean / max richness as blue squares ----
rich_min  <- min(div$R, na.rm = TRUE)
rich_mean <- round(mean(div$R, na.rm = TRUE))
rich_max  <- max(div$R, na.rm = TRUE)

rich_stats <- tibble::tibble(
  site_stat = c("Min", "Mean", "Max"),
  richness  = c(rich_min, rich_mean, rich_max)
)

ggplot(rich_stats, aes(values = richness, fill = "Richness")) +
  geom_waffle(n_rows = 10, color = "white", size = 0.1) +
  scale_fill_manual(values = c(Richness = "deepskyblue")) +
  facet_wrap(~ site_stat, nrow = 1) +
  theme_void() +
  theme(legend.position = "none")




























