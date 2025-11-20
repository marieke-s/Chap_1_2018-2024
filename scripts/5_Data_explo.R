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


# mtdt_7_sel_v1.1 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

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

# Add area_km2 and dist_seabed_depthsampling, bathy_mean, bathy_max to mtdt
mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, area_km2, dist_seabed_depthsampling)
    , by = "replicates"
  )

# Replace month number by month name
mtdt <- mtdt %>%
  mutate(
    month = factor(month.abb[month], levels = month.abb)
  )

mtdt$season <- factor(mtdt$season,
                      levels = c("Winter", "Spring", "Summer", "Autumn"))

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



map_index_plots(df = pred_raw, 
                version = "ppredictors_raw_num_v1.2_SEASON", 
                output_directory = "./figures/Predictors/Map_Hist/by_season", 
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
#--- R x region [to do] -----


#--**************************************************************************************************** ******** AUTRES CODES *******  ---------------------------------------------------------------------------------------------------------------

#--***********************************************************************************************************************************------
#------------------------------------------------------------ Resp var distribution ---------------------

#-------------------------------------------------------------- TEST BlockCV ------------------
#------- cv_spatial_autocor ------------------
# /!\ cv_spatial_autocor from BlockCV works only with presence-absence response variable

# Make pred a spatial object
spat_pred <- st_as_sf(pred, coords = c("x", "y"), crs = 2154) # Projected coordinate system for France
str(spat_pred)

spat_pred <- spat_pred %>% dplyr::select("species_richness")
head(spat_pred)
spat_pred_sp <- as(spat_pred, "Spatial")
head(spat_pred_sp)

# Plot spat_pred 
plot(spat_pred, max.plot = 13)

# # Define cols 
# cols <- colnames(spat_pred)
# 
# # Remove species_richness and geometry
# cols <- cols[!cols %in% c("species_richness")]
# cols <- cols[!cols %in% c("geometry")]

# cv_spatial_autocor
testblocks <-  cv_spatial_autocor(x = spat_pred_sp, column = "species_richness")

# generate a "occ" column" with random 0 and 1 values
set.seed(123)
spat_pred$occ <- sample(0:1, nrow(spat_pred), replace = TRUE)
head(spat_pred)
spat_pred <- spat_pred %>% dplyr::select(-species_richness)





#------- cv_block_size ---------------------
cv_block_size(x = spat_pred)

#------- block size manually ---------------------
# Compute extent of the spatial object
ext <- st_bbox(spat_pred)

# Compute extent in km
extent_x_km <- (ext["xmax"] - ext["xmin"]) / 1000
extent_y_km <- (ext["ymax"] - ext["ymin"]) / 1000

cat("Extent width (km):", extent_x_km, "\n")
cat("Extent height (km):", extent_y_km, "\n")

# Compute expected number of X km blocks
X <- 100
num_blocks_x <- extent_x_km / X
num_blocks_y <- extent_y_km / X

cat("Expected number of blocks:", floor(num_blocks_x * num_blocks_y), "\n")

rm(ext, num_blocks_x, num_blocks_y)

# Create grid with cell size of X km (X * 1000 meters)
grid <- st_make_grid(spat_pred, cellsize = X * 1000, what = "polygons")

# Convert to sf object
grid_sf <- st_sf(geometry = grid)

# Ensure CRS matches spat_pred
st_crs(grid_sf) <- st_crs(spat_pred)

# Print number of blocks
cat("Number of blocks created:", nrow(grid_sf), "\n")

# Plot the grid and spatial points
ggplot() +
  geom_sf(data = grid_sf, fill = NA, color = "blue", linewidth = 0.5) +  # Blocks in blue
  geom_sf(data = spat_pred, color = "red", size = 2) +  # Points in red
  theme_minimal() +
  labs(title = "10 km Spatial Blocks Over Study Area",
       subtitle = paste("Total Blocks:", nrow(grid_sf)),
       caption = "Blue = Blocks, Red = Spatial Points")




#------- cv_spatial ---------------------
cur_folds <- cv_spatial(
  x = pa_data, # Données présence-absence (objet créé dans ce script au dessus)
  column = "Observed", # Colonne contenant l'info présence/absence
  r = baseline[["bio5"]], # Facultatif, quel fond de carte utiliser pour le graphe
  size = blocksize, # Taille des blocs
  k = nb_of_folds, # Nombre de folds
  selection = "random", # Méthode d'attribution des blocs en folds.Cf. ?spatialBlock
  iteration = 50, # Nombre de runs pour répartir les blocs en folds. Nombre d'essais pour trouver des plis équilibrés --> prend le meilleur, le + équilibré
  biomod2 = TRUE) # Créer une table au format biomod2 ici


# Test on our data
cur_folds <- cv_spatial(
  x = spat_pred, # Données présence-absence (objet créé dans ce script au dessus)
  column = "species_richness", # Colonne contenant l'info présence/absence
  size = 10000, # Taille des blocs
  k = 5, # Nombre de folds
  selection = "random", # Méthode d'attribution des blocs en folds.Cf. ?spatialBlock
  iteration = 50, # Nombre de runs pour répartir les blocs en folds. Nombre d'essais pour trouver des plis équilibrés --> prend le meilleur, le + équilibré
  biomod2 = TRUE) # Créer une table au format biomod2 ici

str(cur_folds)

#------ Conclusion Block CV ------
# BlockCV doesnt appear to function on our data. Manually created blocks are uneven. So we will use the Cyril Hautecoeur's function to create the folds.




#-------------------------- 1. Get spatial autocorrelation range ------------------







# #------- manually - species richness alone ------------------
# # Compute empirical variogram
# vgm_emp <- variogram(species_richness ~ 1, data = spat_pred_sp)
# 
# # Plot the variogram
# plot(vgm_emp)
# 
# # Fit a variogram model
# vgm_fit <- fit.variogram(vgm_emp, model = vgm("Sph"))
# 
# # Plot fitted variogram
# plot(vgm_emp, model = vgm_fit)
# 
# # Create a spatial weights matrix
# coords <- st_coordinates(spat_pred)
# nb <- dnearneigh(coords, d1 = 0, d2 = max(dist(coords)) / 3)
# lw <- nb2listw(nb, style = "W")
# 
# # Compute Moran's I
# moran_test <- moran.test(spat_pred$species_richness, lw)
# 
# #-- Results --
# print(vgm_emp)
# # np       dist     gamma dir.hor dir.ver   id
# # 1  112   8430.903  88.66964       0       0 var1
# # 2  131  19989.667 125.19847       0       0 var1
# # 3  143  33185.973 135.33217       0       0 var1
# # 4  137  46347.953 120.41241       0       0 var1
# # 5  138  59453.639 182.67029       0       0 var1
# # 6  142  73066.971 143.14789       0       0 var1
# # 7  143  86158.530 169.80420       0       0 var1
# # 8  104  99104.567 188.44712       0       0 var1
# # 9  146 112245.658 203.54452       0       0 var1
# # 10 111 126161.573 113.92342       0       0 var1
# # 11 131 138893.768 129.83969       0       0 var1
# # 12  84 151563.871 184.10119       0       0 var1
# # 13  82 165376.495 138.86585       0       0 var1
# # 14  94 178142.020 103.36702       0       0 var1
# # 15 123 191668.327 178.76423       0       0 var1
# 
# 
# # Interpretation
# # Your empirical variogram contains three key columns:
# #   
# # - dist: The average distance between point pairs in each bin (in meters).
# # - gamma: The semivariance, representing the dissimilarity between points at that distance.
# # - np: The number of point pairs contributing to each semivariance estimate.
# # 
# # 🔹 Observations from Your Empirical Variogram:
# #   
# #   At short distances (~8 400 m), the semivariance (gamma) is 88.67, meaning points close to each other have similar species richness.
# # As distance increases (e.g., ~199 000 m), semivariance increases to 178.76, meaning points become more dissimilar.
# # The semivariance fluctuates, indicating some spatial correlation but not a perfect trend.
# 
# 
# 
# 
# print(vgm_fit)
# # model    psill    range
# # 1   Nug 61.26653     0.00
# # 2   Sph 86.83085 39096.16
# 
# 
# # Interpretation
# # This fitted spherical variogram model describes spatial autocorrelation.
# # 
# # 🔹 Key Parameters
# # 1️⃣ Nugget (psill = 61.27)
# # 
# # Represents spatial variation at very short distances (measurement error or micro-scale variation).
# # Since the nugget is large, it suggests a high proportion of variation is unexplained by spatial structure (species richness varies a lot at short distances).
# # 
# # 2️⃣ Partial Sill (psill = 86.83)
# # 
# # Represents the structured spatial variation (the amount of variation explained by spatial autocorrelation).
# # 86.83 indicates spatial structure in richness distribution.
# #
# # 3️⃣ Range (range = 39,096.16 m ≈ 39 km)
# # 
# # The range is the distance at which spatial autocorrelation becomes negligible.
# # After ~39 km, points are no longer correlated in terms of species richness.
# # 
# # 🔹 Interpretation of Variogram Fit
# # ✅ There is spatial autocorrelation up to ~39 km, beyond which species richness varies randomly.
# # ✅ The high nugget effect suggests considerable local variability or measurement noise.
# 
# 
# 
# 
# print(moran_test)
# # Moran I test under randomisation
# # 
# # data:  spat_pred$species_richness  
# # weights: lw    
# # 
# # Moran I statistic standard deviate = -0.37761, p-value = 0.6471
# # alternative hypothesis: greater
# # sample estimates:
# #   Moran I statistic       Expectation          Variance 
# # -0.018238979      -0.011111111       0.000356305 
# 
# 
# # Interpretation
# # 🔹 What This Means:
# #   
# #   Moran’s I statistic = -0.0182 → Very close to 0, suggesting a random spatial pattern (no strong clustering or dispersion).
# # p-value = 0.6471 → Not significant (p > 0.05), meaning species richness does not exhibit significant spatial autocorrelation.
# # Z-score = -0.3776 → Also near 0, confirming the absence of strong spatial structure.
# 
# 
# 
#  
# # Comparison Morand's I vs. Variogram Analysis
# # 
# # The variogram suggested some spatial autocorrelation up to ~39 km, but Moran’s I shows that the overall global spatial pattern is weak.
# # This could mean local spatial trends exist, but the overall spatial pattern is mostly random.
# 
# 
# 
# 
# 
# 
# 

#------- manually for all variables - METHOD n°1 ------------------
# Extract numeric data from SpatialPointsDataFrame
spat_pred_data <- spat_pred_sp@data

# Identify numeric variables
var_names <- names(spat_pred_data)[sapply(spat_pred_data, is.numeric)]

# Initialize list to store variogram and Moran's I results
results_list <- list()

# Create spatial weights matrix for Moran's I
coords <- st_coordinates(spat_pred)
nb <- dnearneigh(coords, d1 = 0, d2 = max(dist(coords)) / 3)
lw <- nb2listw(nb, style = "W")

# Loop through each variable
for (var in var_names) {
  
  # Compute empirical variogram
  vgm_emp <- variogram(as.formula(paste(var, "~ 1")), data = spat_pred_sp)
  
  # Fit a spherical variogram model
  vgm_fit <- tryCatch(
    fit.variogram(vgm_emp, model = vgm("Sph")),
    error = function(e) NULL  # Handle cases where fitting fails
  )
  
  # Compute Moran’s I for spatial autocorrelation
  moran_test <- tryCatch(
    moran.test(spat_pred_sp[[var]], lw),
    error = function(e) NULL  # Handle cases where Moran's I fails
  )
  
  # Store results if both computations succeed
  if (!is.null(vgm_fit) && !is.null(moran_test)) {
    results_list[[var]] <- data.frame(
      variable = var,
      model = vgm_fit$model[2],  
      nugget = vgm_fit$psill[1],  
      partial_sill = vgm_fit$psill[2],  
      total_sill = sum(vgm_fit$psill),  
      range = vgm_fit$range[2],  
      moran_i = moran_test$estimate[1],  # Moran's I statistic
      moran_p = moran_test$p.value       # p-value for significance
    )
  }
}

# Convert results list to dataframe
results_df_1 <- bind_rows(results_list)

# Print final table
print(results_df_1)



#------- manually for all variables - METHOD n°2 ------------------
# Extract numeric data from SpatialPointsDataFrame
spat_pred_data <- spat_pred_sp@data

# Identify numeric variables
var_names <- names(spat_pred_data)[sapply(spat_pred_data, is.numeric)]

# Initialize list to store variogram and Moran's I results
results_list <- list()

# Create spatial weights matrix for Moran's I
coords <- st_coordinates(spat_pred)
nb <- dnearneigh(coords, d1 = 0, d2 = max(dist(coords)) / 3)
lw <- nb2listw(nb, style = "W")

# Loop through each variable
for (var in var_names) {
  
  # Compute empirical variogram using autoKrige
  auto_vario <- tryCatch(
    autoKrige(as.formula(paste(var, "~ 1")), spat_pred_sp),
    error = function(e) NULL  # Handle cases where variogram fitting fails
  )
  
  # Compute Moran’s I for spatial autocorrelation
  moran_test <- tryCatch(
    moran.test(spat_pred_sp[[var]], lw),
    error = function(e) NULL  # Handle cases where Moran's I fails
  )
  
  # Store results if both computations succeed
  if (!is.null(auto_vario) && !is.null(moran_test)) {
    vgm_fit <- auto_vario$var_model  # Extract fitted variogram model
    
    results_list[[var]] <- data.frame(
      variable = var,
      model = vgm_fit$model[2],   # Automatically selected model type
      nugget = vgm_fit$psill[1],   # Nugget effect
      partial_sill = vgm_fit$psill[2],  # Partial sill
      total_sill = sum(vgm_fit$psill),  # Total sill = Nugget + Partial Sill
      range = vgm_fit$range[2],   # Spatial range
      moran_i = moran_test$estimate[1],  # Moran's I statistic
      moran_p = moran_test$p.value       # p-value for significance
    )
  } else {
    message(paste("❌ No valid variogram or Moran's I computation for", var))
  }
}

# Convert results list to dataframe
results_df_2 <- bind_rows(results_list)

# Print final table
print(results_df_2)



#------- Mean distance between points [ TO CHECK ]------------------

# Compute pairwise distances (in meters, assuming projected CRS)
dist_matrix <- geosphere::distm(st_coordinates(spat_pred_4326))

# Extract upper triangle (excluding self-distances)
dist_values <- dist_matrix[upper.tri(dist_matrix)]

# Compute mean distance
mean(dist_values) # 237 km --> High value because we have points in Corsica and France
boxplot(dist_values/1000)
max(dist_values/1000)
min(dist_values/1000) # weid... --> CHECK

#------- Visualization ------------------
#---- Map point with color scale for col ----
var <- c("gravity_mean", "number_habitat", "CHL_week_max")

for (i in var) {
  col <- i 
  p <- ggplot() +
    geom_sf(data = spat_pred, aes(color = .data[[col]])) +  # Fix here
    scale_color_viridis_c() +
    theme_minimal() +
    labs(title = "Spatial Distribution",
         subtitle = paste0("Color scale represents ", col),
         caption = "Data Source: eDNA Sampling Points")
  
  print(p)  # Explicitly print each plot inside the loop
}

rm(i, var)

hist(spat_pred$gravity_mean, breaks = 90)

#---- Fig Autocorrelation range ----
# Create data frames for plotting, ensuring Moran's I results are included
range_data_1 <- results_df_1 %>% 
  select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 1: Spherical")

range_data_2 <- results_df_2 %>% 
  select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 2: AutoKrige")

# Combine both datasets and process
range_data <- bind_rows(range_data_1, range_data_2) %>%
  mutate(
    range_km = range / 1000,  # Convert meters to km
    sig_stars = case_when(
      moran_p < 0.001 ~ "***",
      moran_p < 0.01  ~ "**",
      moran_p < 0.05  ~ "*",
      TRUE            ~ "NS"  # Prevent empty string errors
    )
  )

# ✅ Ensure variables are ordered by **descending** range
range_data <- range_data %>%
  arrange(desc(range_km)) %>%
  mutate(variable = factor(variable, levels = rev(unique(variable))))  # Reverse order

# Compute summary statistics for **each method**
summary_stats <- range_data %>%
  group_by(method) %>%
  summarise(
    min_range = min(range_km, na.rm = TRUE),
    max_range = max(range_km, na.rm = TRUE),
    median_range = median(range_km, na.rm = TRUE),
    mean_range = mean(range_km, na.rm = TRUE)
  ) %>%
  ungroup()

# ✅ Format summary text for annotation
summary_text <- summary_stats %>%
  mutate(
    text = paste0(method, ": ",
                  "\nMin: ", round(min_range, 2), " km",
                  " | Max: ", round(max_range, 2), " km",
                  "\nMedian: ", round(median_range, 2), " km",
                  " | Mean: ", round(mean_range, 2), " km\n")
  ) %>%
  pull(text) %>%
  paste(collapse = "\n")  # Combine into single text block

# Generate bar plot
range_plot <- ggplot(range_data, aes(x = variable, y = range_km, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Adjust width for spacing
  coord_flip() +  # Flip axes for better readability
  theme_light() +
  labs(title = "Comparison of Spatial Autocorrelation Ranges",
       x = "Variable",
       y = "Range (km)",
       fill = "Method") +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

# Custom filling colors
range_plot <- range_plot +
  scale_fill_manual(
    name = "Method",  # Legend title
    values = c("Method 1: Spherical" = "#003049", 
               "Method 2: AutoKrige" = "#669bbc")
  )


# ✅ Add Moran's I values and significance stars next to each bar
range_plot <- range_plot +
  geom_text(aes(label = paste0("I: ", round(moran_i, 2), " ", sig_stars), 
                color = sig_stars),  # Use sig_stars for color
            position = position_dodge(width = 0.7), 
            hjust = -0.3, size = 4) +
  
  # ✅ Fix: Use valid names in scale_color_manual
  scale_color_manual(
    name = "Moran's I Significance",
    values = c("***" = "red", "**" = "orange", "*" = "blue", "NS" = "black"),
    labels = c("*** p < 0.001", "** p < 0.01", "* p < 0.05", "NS")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 4)))  # Make legend symbols bigger

# ✅ Add summary statistics as a **text box** at the bottom
range_plot <- range_plot +
  annotate("text", x = length(unique(range_data$variable)) -12, 
           y = min(range_data$range_km) + 3000,  # Place below the lowest bar
           label = summary_text, size = 4, hjust = 0, color = "black",
           fontface = "bold", alpha = 0.8)

# Print the final plot
print(range_plot)





# To do 
# 
# mk : 
# change colors
# modify the legend


#------- Results and Interpretation ------------------
# Method

# We used different methods to compute the spatial autocorrelation range for each variable: 
# You applied two different variogram fitting methods:
#   
# Manual Method (gstat::fit.variogram() with Sph model)
# Automated Method (automap::autoKrige(), which selects the best model automatically, as in BlockCV)
# 
# The results show large differences between the two methods in terms of:
#   
#   Selected models (e.g., "Sph" vs. "Ste" or "Gau")
# Estimated range, nugget, and sill values








# Results

# Define the output file path
output_file <- "variogram_moran_interpretation.txt"

# Create a text string with the interpretation
interpretation <- "
### Interpretation and Comparison of Spatial Autocorrelation Analysis

You performed two methods to analyze spatial autocorrelation:

1. **Method 1:** Using `variogram()` and manually fitting a spherical variogram model (`vgm(\"Sph\")`).
2. **Method 2:** Using `autoKrige()`, which **automatically** selects and fits a variogram model.

Both methods include:
- **Empirical variogram fitting** (to estimate spatial autocorrelation structure)
- **Moran's I test** (to assess spatial dependence)

---

### Key Variogram Parameters:
- **Nugget (`nugget`)**: Variance at zero distance (measurement error or microscale variation).
- **Partial Sill (`partial_sill`)**: The spatially structured variance.
- **Total Sill (`total_sill = nugget + partial_sill`)**: The total variance.
- **Range (`range`)**: The distance where spatial correlation becomes negligible.
- **Moran’s I (`moran_i`)**: Measures the strength of spatial autocorrelation.
- **Moran's p-value (`moran_p`)**: Determines if spatial autocorrelation is significant.

---

### Comparing Variogram Results:
| Variable           | Method 1 - Spherical (Manually Fitted) | Method 2 - AutoKrige (Automatically Selected) | Interpretation |
|--------------------|--------------------------------------|--------------------------------------|------------------------------|
| R                 | Sph, Range = 39 km | Ste, Range = 38 km | Similar range, AutoKrige selects Stein model. |
| FD                | Sph, Range = 19.6 km | Sph, Range = 14.4 km | AutoKrige detects shorter range (local correlation). |
| LFI               | Sph, Range = 28.7 km | Ste, Range = 39.7 km | AutoKrige finds a slightly longer correlation range. |
| Crypto            | Sph, Range = 85.9 km | Ste, Range = 40.6 km | AutoKrige finds a much shorter correlation range. |
| DP_B_ratio        | Sph, Range = 63.8 km | Sph, Range = 6.1 km | AutoKrige detects much shorter correlation. |
| RedList           | Sph, Range = 29.8 km | Gau, Range = 24.3 km | AutoKrige detects slightly shorter correlation. |
| Chondri           | Sph, Range = 32.3 km | Gau, Range = 3.1 km  | AutoKrige finds a very localized correlation pattern. |
| Commercial        | Sph, Range = 31 km  | Ste, Range = 35.5 km | Similar results for both methods. |
| High_commerc      | Sph, Range = 3120 km | Ste, Range = 193 km  | AutoKrige finds a much shorter correlation range. |
| PD                | Sph, Range = 182.9 km | Ste, Range = 5 km  | AutoKrige suggests very localized dependence. |
| Vulner            | Sph, Range = 57.4 km  | Ste, Range = 5.6 km  | AutoKrige detects much shorter correlation. |
| dist_to_reserve   | Sph, Range = 86.6 km  | Gau, Range = 30.8 km | AutoKrige detects shorter correlation. |
| gravity_mean      | Sph, Range = 123.1 km | Ste, Range = 6.3 km  | AutoKrige suggests much faster decay of correlation. |
| bathy_max_log     | Sph, Range = 1265 km | Ste, Range = 2435 km | AutoKrige finds a much longer correlation range. |
| bathy_min_log     | Sph, Range = 132 km | Gau, Range = 6.5 km | AutoKrige detects short-range correlation. |
| dist_shore        | Sph, Range = 187 km | Gau, Range = 3.7 km | AutoKrige finds a very localized correlation pattern. |
| number_habitat    | Sph, Range = 1821 km | Ste, Range = 3754 km | Both methods find long-range correlation. |
| SAL_w_min         | Sph, Range = 193 km | Ste, Range = 41.6 km | AutoKrige detects shorter correlation. |
| SAL_w_max         | Sph, Range = 63.9 km  | Gau, Range = 5 km  | AutoKrige detects much shorter correlation. |
| TEMP_w_min        | Sph, Range = 163 km | Gau, Range = 2.4 km  | AutoKrige suggests very localized dependence. |
| TEMP_w_max        | Sph, Range = 141 km | Ste, Range = 144 km | Similar results for both methods. |
| CHL_week_max      | Sph, Range = 339 km | Sph, Range = 191 km | Both use spherical, but AutoKrige finds a shorter range. |
| CURRENT_w_max     | Sph, Range = 57 km  | Sph, Range = 53 km  | Similar results. |
| WIND_w_max        | Sph, Range = 19 km  | Ste, Range = 10 km  | Similar results. |

#### Key Observations:
1. **Manually fitted variograms (`Sph`) tend to estimate longer ranges.**
2. **AutoKrige often detects much shorter correlation ranges,** especially for environmental predictors.
3. **Variables with high variance (`bathy_max_log`, `number_habitat`) show large differences between methods.**
4. **For some variables (`CHL_week_max`, `WIND_w_max`), both methods give similar results.**

---

### Moran's I Interpretation:
| Variable           | Method 1 Moran’s I | p-value | Method 2 Moran’s I | p-value | Interpretation |
|--------------------|------------------|---------|------------------|---------|------------------------------|
| R                 | -0.0182 | 0.647 | -0.0182 | 0.647 | No significant spatial autocorrelation. |
| FD                | -0.0229 | 0.732 | -0.0229 | 0.732 | No significant spatial autocorrelation. |
| LFI               | -0.0274 | 0.804 | -0.0274 | 0.804 | No significant spatial autocorrelation. |
| DP_B_ratio        |  0.1566 | <0.001 |  0.1566 | <0.001 | Strong spatial autocorrelation. |
| RedList           | -0.0091 | 0.458 | -0.0091 | 0.458 | No spatial autocorrelation. |
| Chondri           | -0.0032 | 0.335 | -0.0032 | 0.335 | No spatial autocorrelation. |
| gravity_mean      |  0.1802 | <0.001 |  0.1802 | <0.001 | Strong positive spatial autocorrelation. |
| bathy_max_log     |  0.0565 | <0.001 |  0.0565 | <0.001 | Weak but significant spatial autocorrelation. |
| number_habitat    |  0.3021 | <0.001 |  0.3021 | <0.001 | Strong spatial autocorrelation. |
| CHL_week_max      |  0.2088 | <0.001 |  0.2088 | <0.001 | Strong spatial autocorrelation. |

✅ **Recommendation**: AutoKrige provides more localized correlations, which may be preferable for capturing finer-scale spatial patterns.

"

# Save interpretation to a text file
writeLines(interpretation, output_file)

# Confirmation message
cat("Interpretation saved to:", output_file, "\n")











#-------------------------- Choice of buffer size based on results ------------------
# Define the text content
buffer_cv_recommendation <- "
### Recommended Buffer Size for Buffer Cross-Validation

Based on the spatial autocorrelation analysis from the variogram and Moran’s I results, we need to determine an appropriate buffer size for spatial cross-validation. The buffer should be large enough to remove spatial dependence but not too large to avoid excessive data loss.

---

### 🔹 Selecting an Optimal Buffer Size
📌 **Guiding Principle:** The buffer size should be **at least half the range** of spatial autocorrelation to effectively reduce dependency in training and test sets.

#### **Step 1: Identify Key Spatial Ranges**
- The **range** values from the variograms indicate the distance where spatial correlation becomes negligible.
- AutoKrige tends to find **shorter spatial correlations**, whereas manual fitting (`Sph`) often detects **longer correlations**.

#### **Step 2: Choose a Representative Distance**
- **Short-range autocorrelation (local effects):**  
  - `gravity_mean` → **6.3 km**  
  - `bathy_min_log` → **6.5 km**  
  - `dist_shore` → **3.7 km**  
  - `SAL_w_max` → **5 km**  
  - `TEMP_w_min` → **2.4 km**  
- **Mid-range autocorrelation (regional effects):**  
  - `FD` → **14.4 km**  
  - `LFI` → **39.7 km**  
  - `Commercial` → **35.5 km**  
  - `CHL_week_max` → **191 km**  
- **Long-range autocorrelation (broad-scale patterns):**  
  - `number_habitat` → **3754 km**  
  - `bathy_max_log` → **2435 km**  
  - `High_commerc` → **193 km**  

#### **Step 3: Select a Balanced Buffer Size**
1. If we **focus on local effects**, a **buffer of ~10 km** would be reasonable, as most local variables have correlation ranges below 15 km.  
2. If we want **regional independence**, **30–50 km** could be a good choice.  
3. **For large-scale patterns**, a **100+ km buffer** may be necessary, but this would drastically reduce training data.

✅ **Final Recommendation:**  
- If the goal is to **eliminate short-range spatial dependence**, use **10–15 km**.  
- If we need **moderate spatial independence**, use **30 km**.  
- If large-scale independence is required, use **50–100 km** (but be cautious about data loss).

---"

writeLines(buffer_cv_recommendation, "buffer_cv_recommendation.txt")








#------- Clean env ------------------
rm(list = ls())
#--- *********************************************************************************************************************************** -----