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
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0.csv") %>%
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


#------------------------------------- Explore Mtd_7 ------------------------------------------------------------------------------------------------------------------------

#--- Check NAs -----
sapply(mtdt, function(x) sum(is.na(x)))






######## Temporal variability #####
#--- Data prep ----
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




#--- Map + hist season | month | year -----

# Cols to plot 
cat_cols <- c("year", "month", "season")

map_categorical_plots(mtdt,
                      cols_to_plot = cat_cols,
                      version = "mtdt_7_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist", 
                      separate_maps = TRUE)



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
#--- Data prep ----
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


#--- Map and hist of sampling effort -----
cols_to_plot <- c("area_km2", "dist_seabed_depthsampling", "depth_sampling", "bathy_mean", "bathy_max", "PCR_replicates", "field_replicates", "estimated_volume_total")

map_index_plots(df = mtdt, 
                version = "mtdt_7_sampling_effort", 
                output_directory = "./figures/Mtdt/Map_Hist", 
                cols_to_plot = cols_to_plot)



colnames(mtdt)














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















#------------------------------------- Explore raw predictors ---------------------------------------------------------------------------------------------------------------
#--- Data prep ----
# numerical cols
cols <- pred_raw[2:150] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
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

#--- Normality [TODO] ----
#-------------------------- 0. CHECKS : normality --------------------------
# Explanation : We check if predictors are normally distributed in order to determine which correlation test to use.

# Initialize results dataframe
normality_results <- data.frame(Column = character(), Shapiro_p = numeric(), KS_p = numeric(), Normality = character(), stringsAsFactors = FALSE)

# Loop through each column
for (col in colnames(num)) {
  cat("\nColumn:", col, "\n")
  
  shapiro_res <- shapiro.test(num[[col]])
  ks_res <- ks.test(num[[col]], "pnorm", mean(num[[col]]), sd(num[[col]]))
  
  print(shapiro_res)
  print(ks_res)
  
  # Determine normality based on Shapiro-Wilk test
  normality_status <- ifelse(shapiro_res$p.value > 0.05, "Normal", "Not-Normal")
  
  # Store results in dataframe
  normality_results <- rbind(normality_results, 
                             data.frame(Column = col, 
                                        Shapiro_p = shapiro_res$p.value, 
                                        KS_p = ks_res$p.value, 
                                        Normality = normality_status))
}

# Check and save results
print(normality_results)
write.csv(normality_results, "./output/predictors_normality_test.csv", row.names = FALSE)

# Interpretation and Conclusion : 
# The majority of predictors are not normally distributed.
# We will thus need to use non parametric test to check correlation among predictors --> use spearman correlation instaed of pearson correlation.

#--- Spatial autocorrelation ----
















#------------------------------------- Explore div indices --------------
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


























# Div x sampling effort | region | season-month | depth | habitat | etc. [TODO] ---


#------------------------------------- Useful bits of code  ---------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------ Resp var distribution ---------------------
#----------- CHECKS : normality + poisson --------------------------
# Explanation : We check if response variables are normally distributed in order to determine which distribution to use in glm.
# Documentation for normality tests : https://www.geeksforgeeks.org/how-to-test-for-normality-in-r/

# Initialize results dataframe
normality_results <- data.frame(Column = character(), Shapiro_p = numeric(), KS_p = numeric(), Normality = character(), stringsAsFactors = FALSE)

for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Compute normality tests
  shapiro_res <- shapiro.test(ind[[col]])
  
  # Handle ties issue in KS test
  ks_res <- tryCatch({
    ks.test(ind[[col]], "pnorm", mean(ind[[col]]), sd(ind[[col]]))
  }, warning = function(w) {
    cat("⚠️ Warning: Ties detected in KS test. Using jittered values.\n")
    ks.test(jitter(ind[[col]]), "pnorm", mean(ind[[col]]), sd(ind[[col]]))
  })
  
  ad_res <- ad.test(ind[[col]])   # Anderson-Darling test
  dag_res <- agostino.test(ind[[col]])  # D’Agostino test
  
  # Print results
  print(shapiro_res)
  print(ks_res)
  print(ad_res)
  print(dag_res)
  
  # Determine normality based on multiple tests
  normality_status <- ifelse(shapiro_res$p.value > 0.05 & ks_res$p.value > 0.05 & 
                               ad_res$p.value > 0.05 & dag_res$p.value > 0.05, 
                             "Normal", "Not-Normal")
  
  # Store results in dataframe
  normality_results <- rbind(normality_results, 
                             data.frame(Column = col, 
                                        Shapiro_p = shapiro_res$p.value, 
                                        KS_p = ks_res$p.value, 
                                        AD_p = ad_res$p.value, 
                                        DAg_p = dag_res$p.value, 
                                        Normality = normality_status))
}

# Print results
print(normality_results)

# Save results
write.csv(normality_results, "./output/resp-var_normality_test.csv", row.names = FALSE)





#----------- FIGURE : Histograms + Normality test ---------------------
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
    geom_histogram(bins = 15, fill = "gray", color = "black", alpha = 0.6) +
    labs(title = var_name, x = "", y = "Frequency") +
    theme_minimal(base_size = 12) +
    annotate("text", x = max(values, na.rm = TRUE) * 0.7, 
             y = max(hist(values, plot = FALSE)$counts) * 0.9,
             label = paste0("SW p = ", round(p_value, 3), "\n(", interpretation, ")"),
             size = 4, hjust = 0, vjust = 1, color = "black", fontface = "bold",
             fill = "white", label.size = 0.5)
  
  return(p)
}


# Generate histograms for each variable in 'ind'
plots <- lapply(names(ind), function(var) create_hist_plot(var, ind[[var]]))

# Arrange all plots in a grid
grid.arrange(grobs = plots, ncol = 3, top = "Histograms with Shapiro-Wilk Normality Test Results")




#----------- FIGURE : QQ Plot + Poisson test --------------------------
# See doc here on how we check if data is Poisson distributed : https://www.geeksforgeeks.org/how-to-know-if-a-data-follows-a-poisson-distribution-in-r/

# Function to create QQ plot for Poisson distribution with chi-squared test
create_qq_plot <- function(var_name, values) {
  
  # Remove NA values
  values <- na.omit(values)
  
  # Convert non-integer values to nearest integers
  values <- round(values)
  
  # Compute mean & variance
  lambda <- mean(values)
  variance <- var(values)
  
  # Observed frequencies
  obs_freq <- table(values)
  
  # Expected frequencies based on Poisson distribution
  exp_freq <- dpois(as.numeric(names(obs_freq)), lambda) * length(values)
  
  # Perform Chi-Squared Test with Monte Carlo approximation
  chisq_test <- tryCatch(
    chisq.test(obs_freq, p = exp_freq, rescale.p = TRUE, simulate.p.value = TRUE, B = 10000),
    error = function(e) return(NULL)
  )
  
  # Interpretation
  if (!is.null(chisq_test)) {
    p_value <- chisq_test$p.value
    interpretation <- ifelse(p_value >= 0.05, "Poisson Distributed", "Not Poisson Distributed")
  } else {
    p_value <- NA
    interpretation <- "Chi-Squared Test Failed"
  }
  
  # Generate QQ Plot
  qqplot_data <- data.frame(
    x = qpois(ppoints(length(values)), lambda),
    y = sort(values)
  )
  
  # Fixed annotation position using `Inf` (top-right of plot)
  p <- ggplot(qqplot_data, aes(x = x, y = y)) +
    geom_point(color = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = var_name, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal(base_size = 12) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Mean: ", round(lambda, 2), 
                            "\nVar: ", round(variance, 2),
                            "\nChi-sq p = ", ifelse(is.na(p_value), "NA", round(p_value, 3)),
                            "\n(", interpretation, ")"),
             hjust = 1, vjust = 1, size = 4, fontface = "bold", color = "black")
  
  return(p)
}






# Generate QQ plots for each variable in 'ind'
qq_plots <- lapply(names(ind), function(var) create_qq_plot(var, ind[[var]]))

# Arrange all QQ plots in a grid using patchwork
final_plot <- patchwork::wrap_plots(qq_plots, ncol = 3) +
  patchwork::plot_annotation(title = "QQ Plots for Poisson Distribution with Chi-Squared Test")

# Display plot
print(final_plot)




#----------- CHECKS : 0 inflation --------------------------
# Load required package
if (!require("pscl")) install.packages("pscl", dependencies = TRUE)
library(pscl)

# Initialize results dataframe
zero_inflation_results <- data.frame(Column = character(), 
                                     Zero_Proportion = numeric(), 
                                     Vuong_p = numeric(), 
                                     Zero_Inflated = character(), 
                                     stringsAsFactors = FALSE)

# Loop over each column
for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Count zeros
  zero_count <- sum(ind[[col]] == 0, na.rm = TRUE)
  zero_proportion <- zero_count / length(ind[[col]])
  
  # Perform Vuong test for zero inflation in Poisson regression
  tryCatch({
    poisson_model <- glm(ind[[col]] ~ 1, family = poisson)
    zip_model <- zeroinfl(ind[[col]] ~ 1, dist = "poisson")
    vuong_test <- vuong(poisson_model, zip_model)  # Vuong test
    
    # Extract p-value
    vuong_p <- vuong_test$p
  }, error = function(e) {
    cat("⚠️ Error: Vuong test failed for", col, "\n")
    vuong_p <<- NA
  })
  
  # Print results
  cat("   - Zero Proportion:", round(zero_proportion, 3), "\n")
  cat("   - Vuong Test p-value:", ifelse(is.na(vuong_p), "Failed", round(vuong_p, 3)), "\n")
  
  # Determine zero inflation based on Vuong test (p < 0.05 suggests zero inflation)
  zero_inflation_status <- ifelse(!is.na(vuong_p) & vuong_p < 0.05, "Zero-Inflated", "Not Zero-Inflated")
  
  # Store results in dataframe
  zero_inflation_results <- rbind(zero_inflation_results, 
                                  data.frame(Column = col, 
                                             Zero_Proportion = zero_proportion, 
                                             Vuong_p = vuong_p, 
                                             Zero_Inflated = zero_inflation_status))
}

# Print final results
print(zero_inflation_results)

#----------- CHECKS : Overdispersion --------------------------
check_overdispersion <- function(dataset, response_vars, predictors) {
  
  # Store all plots
  plots <- list()
  
  # Loop through each response variable
  for (response_var in response_vars) {
    
    cat("\n🔹 Checking Overdispersion for:", response_var, "\n")
    
    # Define formula
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    
    # Fit Poisson Model
    poisson_model <- glm(formula, data = dataset, family = poisson(link = "log"))
    
    # Run Overdispersion Test
    aer_test <- AER::dispersiontest(poisson_model)
    dispersion_ratio <- aer_test$estimate["dispersion"]
    p_value <- aer_test$p.value
    overdisp_status <- ifelse(p_value < 0.05, "Significant Overdispersion", "No Strong Evidence of Overdispersion")
    
    # Generate Residual Plots
    residual_data <- data.frame(
      Fitted = poisson_model$fitted.values,
      Residuals = residuals(poisson_model, type = "pearson")
    )
    
    p <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "blue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Overdispersion Check -", response_var), 
           x = "Fitted Values", y = "Pearson Residuals") +
      theme_minimal(base_size = 12) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold",
               label = paste0(
                 "Dispersion: ", round(dispersion_ratio, 2), 
                 "\nP-value: ", round(p_value, 3),
                 "\n", overdisp_status
               ))
    
    plots[[response_var]] <- p
  }
  
  # Arrange all plots in a grid
  grid.arrange(grobs = plots, ncol = 2, top = "Overdispersion Analysis Across Variables")
}
check_overdispersion <- function(dataset, response_vars, predictors) {
  
  # Store all plots
  plots <- list()
  
  # Initialize a data frame to store results
  results_df <- data.frame(
    Response_Variable = character(),
    Dispersion_Index = numeric(),
    P_Value = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each response variable
  for (response_var in response_vars) {
    
    cat("\n🔹 Checking Overdispersion for:", response_var, "\n")
    
    # Define formula
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    
    # Fit Poisson Model
    poisson_model <- glm(formula, data = dataset, family = poisson(link = "log"))
    
    # Run Overdispersion Test
    aer_test <- AER::dispersiontest(poisson_model)
    dispersion_ratio <- aer_test$estimate["dispersion"]
    p_value <- aer_test$p.value
    overdisp_status <- ifelse(p_value < 0.05, "Significant Overdispersion", "No Strong Evidence of Overdispersion")
    
    # Save results in the data frame
    results_df <- rbind(results_df, data.frame(
      Response_Variable = response_var,
      Dispersion_Index = round(dispersion_ratio, 2),
      P_Value = round(p_value, 3),
      Interpretation = overdisp_status,
      stringsAsFactors = FALSE
    ))
    
    # Generate Residual Plots
    residual_data <- data.frame(
      Fitted = poisson_model$fitted.values,
      Residuals = residuals(poisson_model, type = "pearson")
    )
    
    p <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "blue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Overdispersion Check -", response_var), 
           x = "Fitted Values", y = "Pearson Residuals") +
      theme_minimal(base_size = 12) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold",
               label = paste0(
                 "Dispersion: ", round(dispersion_ratio, 2), 
                 "\nP-value: ", round(p_value, 3),
                 "\n", overdisp_status
               ))
    
    plots[[response_var]] <- p
  }
  
  # Print results table
  print(results_df)
  
  # Arrange all plots in a grid
  grid.arrange(grobs = plots, ncol = 2, top = "Overdispersion Analysis Across Variables")
  
  # Return the results dataframe
  return(results_df)
}

# Usage
check_overdispersion(dataset = tot, response_vars = colnames(ind), predictors = predictors)


#----------- CHECKS : mean and variance --------------------------
# Initialize results dataframe
mean_variance_results <- data.frame(Column = character(), 
                                    Mean = numeric(), 
                                    Variance = numeric(), 
                                    Overdispersed = character(), 
                                    stringsAsFactors = FALSE)

# Loop over each column
for (col in colnames(ind)) {
  cat("\n🔹 Column:", col, "\n")
  
  # Compute mean & variance
  mean_val <- mean(ind[[col]], na.rm = TRUE)
  var_val <- var(ind[[col]], na.rm = TRUE)
  
  # Print results
  cat("   - Mean:", round(mean_val, 3), "\n")
  cat("   - Variance:", round(var_val, 3), "\n")
  
  # Determine overdispersion (Variance > Mean)
  overdispersion_status <- ifelse(var_val > mean_val, "Yes", "No")
  
  # Store results in dataframe
  mean_variance_results <- rbind(mean_variance_results, 
                                 data.frame(Column = col, 
                                            Mean = mean_val, 
                                            Variance = var_val, 
                                            Overdispersed = overdispersion_status))
}

# Print final results
print(mean_variance_results)



#----------- DOC : quasipoisson vs nb distribution ---------------------
# how do I now if my overdispersion is "mild" or "severe" ?
#   ChatGPT said:
#   
#   ​Overdispersion in Poisson regression occurs when the observed variance in the count data exceeds the mean, violating the Poisson model's assumption that the mean equals the variance. Assessing the severity of overdispersion is crucial for selecting an appropriate statistical model. Here's how you can evaluate it:​
# en.wikipedia.org
# 
# 1. Calculate the Dispersion Statistic:
#   
#   The dispersion statistic (also known as the dispersion parameter or index of dispersion) is calculated as:​
# 
# Dispersion Statistic=∑(Yi−Y^i)2/Y^in−pDispersion Statistic=n−p∑(Yi​−Y^i​)2/Y^i​​​
# 
# where YiYi​ are the observed counts, Y^iY^i​ are the predicted counts from the model, nn is the number of observations, and pp is the number of parameters estimated in the model.​
# 
# 2. Interpret the Dispersion Statistic:
#   
#   Close to 1: Indicates that the variance is approximately equal to the mean, suggesting no overdispersion.​
# 
# Greater than 1: Indicates overdispersion.​
# 
# 3. Assess the Severity of Overdispersion:
#   
#   Mild Overdispersion: Dispersion statistic is slightly above 1 (e.g., between 1 and 2). In such cases, the standard errors may be slightly underestimated, but the Poisson model might still be adequate with robust standard errors.​
# 
# Moderate to Severe Overdispersion: Dispersion statistic significantly exceeds 2. This suggests that the Poisson model is likely inappropriate, and alternative modeling approaches should be considered.​
# 
# 4. Alternative Modeling Approaches:
#   
#   Quasi-Poisson Model: Adjusts the standard errors to account for overdispersion without altering the mean structure. Suitable for mild to moderate overdispersion.​
# 
# Negative Binomial Model: Introduces an additional parameter to model the overdispersion explicitly. Preferred for moderate to severe overdispersion.​
# 
# 5. Formal Testing:
#   
#   While calculating the dispersion statistic provides a straightforward assessment, formal statistical tests, such as Dean's test for overdispersion, can also be employed to evaluate the presence and extent of overdispersion. ​
# civil.colorado.edu
# 
# Note: It's essential to consider the context of your data and the implications of overdispersion on your model's inferences. Consulting subject-matter experts or a statistician can provide additional guidance tailored to your specific analysis.



#----------- DOC : over and underdispersion issues --------------------------
# is underdispersion an issue for normally distributed variables ? 
#   ChatGPT said:
#   🔍 Is Underdispersion an Issue for Normally Distributed Variables?
#   
#   No, underdispersion is not an issue for normally distributed variables. The concept of underdispersion primarily applies to count data models (Poisson, Negative Binomial, etc.), where variance is expected to follow a specific structure (e.g., equal to the mean in Poisson models).

# 🔹 Why Doesn't Underdispersion Matter for Normal Distributions?
# 
#     Normal Distribution Assumes Constant Variance
#         In a normal distribution:
#         Y∼N(μ,σ2)
#         Y∼N(μ,σ2) variance (σ2σ2) is a free parameter that is estimated from the data.
#         There's no requirement that variance follows the mean, unlike in Poisson models.
# Whether variance is "too small" or "too large" relative to expectations does not break any assumptions.
# 
# Overdispersion and Underdispersion are Concerns for Poisson-Based Models
# In Poisson models, variance must equal the mean (λλ).
# Underdispersion (Var<MeanVar<Mean) causes issues because standard Poisson models cannot accommodate it.
# Overdispersion (Var>MeanVar>Mean) causes similar problems.
# This is why alternative models (Quasi-Poisson, Negative Binomial, Generalized Poisson) exist.
# 
# In Normal Regression (GLM, GAM, etc.), Variance is Simply Estimated
# In a Gaussian GLM/GAM (family = gaussian()), we assume:
#   Y=Xβ+ε,ε∼N(0,σ2)
# Y=Xβ+ε,ε∼N(0,σ2)
# The model does not expect variance to follow a specific function of the mean, so underdispersion is not problematic.



# 🔍 Is Overdispersion an Issue for Normally Distributed Variables?
#   
#   No, overdispersion is not an issue for normally distributed variables because:
#   
#   The Normal (Gaussian) distribution does not assume a strict relationship between the mean and variance.
# Unlike Poisson or Negative Binomial models, where variance is linked to the mean, a Gaussian model explicitly estimates variance (σ2σ2) as a separate parameter.




# 
# is underdispersion or overdispersion an issue for negative binomial distribution ? 
#   ChatGPT said:
#   🔍 Is Underdispersion or Overdispersion an Issue for Negative Binomial (NB) Distribution?
#   🔹 Overdispersion: NO, NB is designed for it
# 
# Overdispersion is NOT an issue for the Negative Binomial distribution because it is explicitly designed to handle count data where variance exceeds the mean.
# 🔹 Underdispersion: YES, NB is NOT suited for it
# 
# Underdispersion is a problem because the Negative Binomial assumes variance is always greater than the mean:
#   Var(Y)=μ+κμ2
# Var(Y)=μ+κμ2
# 
# where κ>0κ>0 controls the dispersion. This means the variance must be greater than the mean. If the observed variance is less than expected, the NB model will likely perform poorly.
# 🔹 When is Overdispersion or Underdispersion a Problem?
#   Dispersion Type	Impact on NB Model	Solution
# Overdispersion (Var>MeanVar>Mean)	✅ NB models this well	Use glm.nb() or mgcv::nb()
# Underdispersion (Var<MeanVar<Mean)	❌ NB does not fit well	Use Generalized Poisson or COM-Poisson





#----------- DOC : random effect vs residuals ---------------------
#----------- RUN : Assign distribution to each resp var --------------------------

# Create a df to store the distribution of each response variable
distrib <- data.frame(response_var = colnames(ind))

# Assign distribution for each response variable : 
# poisson for count variables,
# gaussian for continuous variables,
# nb for overdispersed count variables,
# (using nb for overdispersion because quasipoisson is better suited for "mild" overdispersion (ie between 1 and 2) and nb for "severe" overdispersion (ie > 2))
distrib$Distribution <- c("nb", "poisson", "nb", "nb", "gaussian", "poisson", "poisson", "nb", "poisson", "gaussian", "gaussian")

distrib





#-------------------------------------------------------------- Check Sampling Bias ------------------
# Renames cols 
dt <- df %>% dplyr::rename("decimalLongitude" = "x", "decimalLatitude" = "y", "species" = "species_richness")

# Set decimalLatitude and decimalLongitude as numeric
dt$decimalLatitude <- as.numeric(dt$decimalLatitude)
dt$decimalLongitude <- as.numeric(dt$decimalLongitude)

# Add/remove sampling rate predictors
# Remove : airports, waterbodies
# Add : ports


# Compute sampling bias
bias <- sampbias::calculate_bias(x = dt, 
                                 gaz = NULL, # use default geographic gazetteers
                                 res = 0.01, 
                                 terrestrial = FALSE)

# Check results
str(bias)
summary(bias)
plot(bias)

proj <- project_bias(bias)
map_bias(proj, type = "log_sampling_rate")


# Interpretation
# Fig 1.
# Strong effect of cities and roads. 
# Negligable effect of airports and waterbodies (expected).


# Fig 2. (map)
# Estimation of sampling rate according to the predictors reflect well the reality. 

# Problem is that sampling rate is computer on all marine space of the extent, when our study actually focuses on the coast. 
# We will need to filter the data to keep only the coastal extent. 
# --> This can be done using restrict_sample parameter.


grid <- st_read("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1/data/processed_data/data_prep/eDNA/Pts_to_Grid/grid_med005_coastal.shp")
class(grid)

bias <- sampbias::calculate_bias(x = dt, 
                                 gaz = NULL, # use default geographic gazetteers
                                 res = 0.01, # 10km resolution
                                 terrestrial = FALSE, 
                                 restrict_sample = grid
)


# Results interpretation with coastal extent
# With coastal extent we see that bias are way lower. 


# Export results
saveRDS(bias, "./output/bias_sampling_coastal.rds")



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