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




#------------- Load functions ------------------------------------------------
source("./utils/Fct_Data-Prep.R")

#------------- Load and prep data ----
# Mtdt_7
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7.gpkg")

# predictors_raw_v1.2
pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v1.2.gpkg")

# div_indices_v1.0
div <- read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))


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


#------------- Explore Mtd_7 ------------------------------------------------------------------------------------------------------------------------

# Check NAs -----
sapply(mtdt, function(x) sum(is.na(x)))






######## Temporal variability #####
# Data prep ----
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




# Map + hist season | month | year -----

# Cols to plot 
cat_cols <- c("year", "month", "season")

map_categorical_plots(mtdt,
                      cols_to_plot = cat_cols,
                      version = "mtdt_7_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)



# Count table season | month | year -----

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
# Data prep ----
# Add area_km2 and dist_seabed_depthsampling, bathy_mean, bathy_max to mtdt
mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, area_km2, dist_seabed_depthsampling)
    , by = "replicates"
  )


mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, bathy_mean, bathy_max)
    , by = "replicates"
  )


# Map and hist of sampling effort -----
cols_to_plot <- c("area_km2", "dist_seabed_depthsampling", "depth_sampling", "bathy_mean", "bathy_max", "PCR_replicates", "field_replicates", "estimated_volume_total")

map_index_plots(df = mtdt, 
                version = "mtdt_7_sampling_effort", 
                output_directory = "./figures/Mtdt/Map_Hist", 
                cols_to_plot = cols_to_plot)



colnames(mtdt)














######## Lockdown ####
# Map + hist lockdown -----
map_categorical_plots(mtdt,
                      cols_to_plot = "lockdown",
                      version = "mtdt_7_lockdown",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)









# Count table season | month | year -----



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
# Map + hist method -----
map_categorical_plots(mtdt,
                      cols_to_plot = "method",
                      version = "mtdt_7_method",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)









# Count table season | month | year -----



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
  filter(lockdown == "0" & season %in% c("Summer", "Spring")) %>%
  filter(field_replicates == 2, PCR_replicates == 24) %>%
  filter(estimated_volume_total >= 55 & estimated_volume_total <= 65) %>%
  filter(bathy_max <= 35) %>%
  filter(method == "surface_transect" | method == "seabed_transect") %>%
  nrow() # 92

mtdt %>%
  filter(method == "seabed_transect") %>%
  mutate(estimated_volume_total = as.numeric(estimated_volume_total)) %>%
  with(hist(estimated_volume_total)) # almost all 120 L 

mtdt %>%
  filter(method == "seabed_transect") %>%
  mutate(estimated_volume_total = as.numeric(PCR_replicates)) %>%
  with(hist(PCR_replicates)) # almost all 24















#------------- Explore raw predictors ---------------------------------------------------------------------------------------------------------------
# Check NAs -----
sapply(pred_raw, function(x) sum(is.na(x))) 
# in predictors_raw_v1.2 : 5 NAs for habitat predictors

# Map + hist + summary : numerical data -----

cols_to_plot <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

map_index_plots(df = pred_raw, version = "predictors_raw_num_v1.2", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = cols_to_plot)


rm(cols_to_plot)


# Map + hist + summary : categorical data -----

cols_to_plot <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.character), where(is.factor)) %>%
  colnames()


map_categorical_plots(df = pred_raw, version = "predictors_raw_cat_v1.2", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = cols_to_plot)



# Map
http://127.0.0.1:18049/graphics/plot_zoom_png?width=1707&height=935
# Map + hist env predictors by season -----
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





# Summary env predictors by season -----
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



# Outliers ----
# numerical cols
cols <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()


# variance-to-mean ratio ------
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



# canyon_shape_area has very high variance-to-mean ratio -----
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

# 1.5xIQR rule -----


# Near zero var -----
# Check near 0 variance variable 
print(caret::nearZeroVar(pred)) # output designate columns 28 and 31
summary(pred[, c(28, 31)]) # they correspond to WIND_w_min and WIND_m_min

# Normality ----
# Spatial autocorrelation ----















#---------- Explore div indices --------------
# Data prep ----
indicators <- div







# Check NAs ----
sapply(indicators, function(x) sum(is.na(x))) 

# Hist + summary ----

## Make a numeric copy of indicators (except 'replicates')
indicators_num <- indicators %>%
  mutate(across(-replicates, ~ as.numeric(.)))



## List of columns to plot
cols_to_plot <- names(indicators_num)[names(indicators_num) != "replicates"]

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


# Map + hist + summary -----

map_index_plots(df = indicators, version = "div_indices_v1.0", output_directory = "./figures/Div_indices/Map_Hist", cols_to_plot = colnames(indicators[4:13]))

























