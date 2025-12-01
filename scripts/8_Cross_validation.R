#------------- Description ---------------------
# Purpose: 
# This script aims to compute select predictors. 

# Author: Marieke Schultz

# Date script created: 16/11/2025

#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(ape)
library(automap)
library(blockCV)
library(dplyr)
library(fastDummies)
library(fishtree)
library(geiger)
library(geosphere)  # For distance calculations
library(ggplot2)
library(gstat)
library(gridExtra)
library(igraph)
library(picante)
library(Rarity)
library(readr)
library(reshape2)
library(rlang)
library(rnaturalearth)  # For land map background
library(sampbias)
library(sf)
library(spdep)
library(svglite)
library(tidyr)
library(units)
library(virtualspecies)



# Load functions
source("./utils/Fct_Data-Prep.R")
source("./utils/Fct_Modeling.R")

#------------- Load and prep data ------------------
# predictors_sel_v.1.1 ----
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.1.gpkg")
# mtdt_7_sel_v1.0 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.0.gpkg")

# div_indices_sel_v1.0.gpkg ----
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.0.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# full ----
full <- pred %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(mtdt), by = "replicates") %>%
  left_join(div, by = "replicates")
#-------------- Check data --------------
# print number of na in each pred column
sapply(pred, function(x) sum(is.na(x))) # no NA
#------------- Predictors spatial autocorrelation -----------------
#--- Prep -----
# For pred
# First make an sf object with 4326 projection
spat_pred <- st_as_sf(pred %>% st_drop_geometry(), coords = c("x", "y"), crs = 2154)
crs(pred) # 2154 = Projected coordinate system for France

# Convert sf object to Spatial object
spat_pred_sp <- as(spat_pred, "Spatial")

# Checks
class(spat_pred_sp)
plot(spat_pred_sp) # punctual data
crs(spat_pred_sp) # 2154, "RGF93 v1 / Lambert-93", LENGTHUNIT[\"metre\",1]]]
spat_pred_sp@coords # [500,]  9.370270  41.63922 --> projected coordinates in meters
spat_pred_sp@data 

# Extract numeric data from SpatialPointsDataFrame
spat_pred_data <- spat_pred_sp@data

# Identify numeric variables
var_names <- names(spat_pred_data)[sapply(spat_pred_data, is.numeric)] %>% 
  setdiff(c("x", "y")) # Exclude non-predictor numeric columns


#--- manually for all predictors - METHOD n°1 ------------------
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

#--- manually for all predictors - METHOD n°2 ------------------
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



#---  Interpretation of results ----
# Interpretation of results: https://chatgpt.com/share/e/691b318f-7c28-800e-ba16-49d7666f9eef
#--- Figure ----
# Create data frames for plotting, ensuring Moran's I results are included
range_data_1 <- results_df_1 %>% 
  dplyr::select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 1: Spherical")

range_data_2 <- results_df_2 %>% 
  dplyr::select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 2: AutoKrige")

# Combine both datasets and process
range_data <- bind_rows(range_data_1, range_data_2) %>%
  mutate(
    range_km = range, 
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
  coord_flip() +  # Flip axes 
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
            hjust = -0.1, size = 3.5) +
  
  # ✅ Fix: Use valid names in scale_color_manual
  scale_color_manual(
    name = "Moran's I Significance",
    values = c("***" = "red", "**" = "orange", "*" = "blue", "NS" = "black"),
    labels = c("*** p < 0.001", "** p < 0.01", "* p < 0.05", "NS")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 4)))  # Make legend symbols bigger

# ✅ Add summary statistics as a **text box** at the bottom
range_plot <- range_plot +
  annotate(
    "text",
    x = length(unique(range_data$variable)) * 0.5,     # middle vertically
    y = 65,               # push to the RIGHT side
    label = summary_text,
    hjust = 0,                                         # left-align the box
    size = 4,
    color = "black",
    fontface = "bold"
  )

# Print the final plot
print(range_plot)







#------------- Div_indices spatial autocorrelation -----------------
#--- Prep -----
# Add coordinates (pred$x and pred$y) to div
div <- div %>%
  left_join(pred %>% st_drop_geometry() %>% dplyr::select(replicates, x, y), by = "replicates")


# First make an sf object with 4326 projection
spat_pred <- st_as_sf(div %>% st_drop_geometry(), coords = c("x", "y"), crs = 2154)

# Convert sf object to Spatial object
spat_pred_sp <- as(spat_pred, "Spatial")

# Checks
class(spat_pred_sp)
plot(spat_pred_sp) # punctual data
crs(spat_pred_sp) # 2154, "RGF93 v1 / Lambert-93", LENGTHUNIT[\"metre\",1]]]
spat_pred_sp@coords # [500,]  9.370270  41.63922 --> projected coordinates in meters
spat_pred_sp@data 

# Extract numeric data from SpatialPointsDataFrame
spat_pred_data <- spat_pred_sp@data

# Identify numeric variables
var_names <- names(spat_pred_data)[sapply(spat_pred_data, is.numeric)] %>% 
  setdiff(c("x", "y")) # Exclude non-predictor numeric columns


#--- manually for all predictors - METHOD n°1 ------------------
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

#--- manually for all predictors - METHOD n°2 ------------------
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


#--- Figure ----
# Create data frames for plotting
range_data_1 <- results_df_1 %>% 
  dplyr::select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 1: Spherical")

range_data_2 <- results_df_2 %>% 
  dplyr::select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 2: AutoKrige")

# Combine both datasets and process
range_data <- bind_rows(range_data_1, range_data_2) %>%
  mutate(
    range_km = range, 
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
  coord_flip() +  # Flip axes 
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
            hjust = -0.1, size = 3.5) +
  
  # ✅ Fix: Use valid names in scale_color_manual
  scale_color_manual(
    name = "Moran's I Significance",
    values = c("***" = "red", "**" = "orange", "*" = "blue", "NS" = "black"),
    labels = c("*** p < 0.001", "** p < 0.01", "* p < 0.05", "NS")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 4)))  # Make legend symbols bigger

# ✅ Add summary statistics as a **text box** at the bottom
range_plot <- range_plot +
  annotate(
    "text",
    x = length(unique(range_data$variable)) * 0.5,     # middle vertically
    y = 65,               # push to the RIGHT side
    label = summary_text,
    hjust = 0,                                         # left-align the box
    size = 4,
    color = "black",
    fontface = "bold"
  )

# Print the final plot
print(range_plot)

beepr::beep()





#------------- Mean distance between points [ Data explo ]------------------

# Compute pairwise distances (in meters, assuming projected CRS)
dist_matrix <- geosphere::distm(st_coordinates(spat_pred))

# Extract upper triangle (excluding self-distances)
dist_values <- dist_matrix[upper.tri(dist_matrix)]

# Compute mean distance
mean(dist_values/1000) # 228 km --> High value because we have points in Corsica and France
boxplot(dist_values/1000)
max(dist_values/1000) # 547 km
min(dist_values/1000) # 0







#-------------  Environmental CV (spBlock_cv) ------------------
# --- Load data  (same as T.0.0) ----
# predictors_sel_v.1.3 ----
pred <- st_read("./data/processed_data/predictors/predictors_sel_v1.3.gpkg")

# remove space and majuscule in habitat names
pred$grouped_main_habitat <- gsub(" ", "_", pred$grouped_main_habitat)
pred$grouped_main_habitat <- tolower(pred$grouped_main_habitat)

# set character as factor
pred$grouped_main_habitat <- as.factor(pred$grouped_main_habitat)

# check levels
levels(pred$grouped_main_habitat)
table(pred$grouped_main_habitat)

unique(pred$grouped_main_habitat)

# mtdt_7_sel_v1.1 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

# div_indices_sel_v1.1.gpkg ----
div <- readr::read_csv2("./data/processed_data/Traits/div_indices_v1.0_sel_v1.1.csv") %>%
  mutate(DeBRa = as.numeric(DeBRa))

# tot ----
tot <- pred %>%
  st_drop_geometry() %>%
  left_join(st_drop_geometry(mtdt), by = "replicates") %>%
  left_join(div, by = "replicates")

# --- Environmental CV  ----
# Extract only environmental features for clustering (exclude label) ---
env_data_only <- pred %>% dplyr::select(-c("x", "y", "replicates", grouped_main_habitat)) %>%
  st_drop_geometry()


# --- Determine optimal number of clusters (k) using elbow method ----
wcss <- numeric(15)
for (i in 1:15) {
  kmeans_result <- kmeans(env_data_only, centers = i, nstart = 15)
  wcss[i] <- kmeans_result$tot.withinss
}
plot(1:15, wcss, type = 'b', xlab = 'Number of Clusters', ylab = 'Within-Cluster SS (WCSS)')


# --- Determine optimal number of clusters (k) vizualization  ----
# Vector of k values
k_vals <- 4:9  

# 1) PCA for visualization (on all variables)
pca <- prcomp(env_data_only, scale. = TRUE)
scores <- pca$x[, 1:2]  # first two principal components

# 2) Set up panel: 2 rows x 3 columns
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

# Color palette with enough colors for maximum k
cols <- rainbow(max(k_vals))

set.seed(123)  # for reproducible kmeans


# K-means PCA plots ---
for (k in k_vals) {
  km <- kmeans(env_data_only, centers = k, nstart = 25)
  
  plot(scores,
       col = alpha(cols[km$cluster], 0.5),
       pch = 19,
       xlab = "PC1",
       ylab = "PC2",
       main = paste("k =", k),
       asp = 1)
}

# K-means maps (panel) ---
for (k in k_vals) {
  # K-means clustering
  km <- kmeans(env_data_only, centers = k, nstart = 25)
  
  # Plot spatial map
  plot(pred$x, pred$y,
       col = alpha(cols[km$cluster], 0.5),
       pch = 19,
       xlab = "x",
       ylab = "y",
       main = paste("k =", k),
       asp = 1)
}

# (Optional) reset plotting parameters
par(op)


# K-means maps ---
# World polygons
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Turn your points into an sf object
# If x/y are longitude/latitude in degrees, CRS 4326 is correct.
# If they’re in another projection, adjust the crs argument accordingly.
pts_sf <- st_as_sf(pred, coords = c("x", "y"), crs = 4326)
pts_sf <- st_transform(pts_sf, crs = st_crs(world))  # Match world CRS

# Bounding box of your data
bbox <- st_bbox(pts_sf)

# Crop world to your data extent
world_crop <- st_crop(world, bbox)

set.seed(123)
k_vals <- 4:9

# K-means maps (distinct maps for k clusters) ---
# Convert your pred df to sf with correct CRS (Lambert-93)
pts_sf <- st_as_sf(pred, coords = c("x", "y"), crs = 2154)

# Transform to WGS84 for mapping
pts_wgs84 <- st_transform(pts_sf, 4326)

# Bounding box in WGS84
bbox <- st_bbox(pts_wgs84)

# World and crop to bbox
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_crop <- st_crop(world, bbox)


# Loop through k = 4 to 9

k_vals <- 4:9
set.seed(123)

# Make sure output directory exists
dir.create("./figures/Cross-val", recursive = TRUE, showWarnings = FALSE)

for (k in k_vals) {
  # Run k-means on env_data_only (same order as pts_wgs84)
  km <- kmeans(env_data_only, centers = k, nstart = 25)
  
  # Attach cluster to the sf object
  pts_plot <- pts_wgs84
  pts_plot$cluster <- factor(km$cluster)
  
  # Panel map: one facet per cluster
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "grey95", color = "grey80") +
    geom_sf(data = pts_plot, aes(color = cluster), size = 0.5) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"])
    ) +
    facet_wrap(~ cluster) +
    labs(
      title = paste("k =", k),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  
  # Save figure
  out_path <- sprintf("./figures/Cross-val/Map_k%d_predictors_T0.0.jpg", k)
  ggsave(out_path, p, width = 10, height = 8, dpi = 300)
  message("Saved: ", out_path)
}

















# --- Check : Geographical distance between points vs predictors PCA distance ----

# Pairwise distances in PCA space
d_pca <- dist(scores)  # Euclidean

# Pairwise geographical distances
d_xy <- dist(as.matrix(pred[, c("x", "y")]))

# Convert to numeric vectors
d_pca_vec <- as.numeric(d_pca)
d_xy_vec  <- as.numeric(d_xy)

# Sample a subset for plotting
set.seed(123)
n_pairs <- length(d_pca_vec) # 202 566 --> bad visibility in plot
sample_idx <- sample(seq_len(n_pairs), size = min(10000, n_pairs)) # sample 5000 pairs for plotting

# Build a clean data frame for the sampled points
df <- data.frame(
  xy = d_xy_vec[sample_idx],
  pca = d_pca_vec[sample_idx]
)

# Fit loess using simple names
lo <- loess(pca ~ xy, data = df)

# Create a grid for prediction
xgrid <- seq(min(df$xy), max(df$xy), length.out = 200)

# Predict on smooth grid
ygrid <- predict(lo, newdata = data.frame(xy = xgrid))

# Plot PCA distance vs. actual distance
plot(d_xy_vec[sample_idx], d_pca_vec[sample_idx],
     xlab = "Distance in actual space (x,y)",
     ylab = "Distance in PCA space (PC1–PC2)",
     main = "Comparison of distances: PCA vs. actual coordinates",
     pch  = 19, cex = 0.5,
     col = alpha("skyblue", 0.5))
lines(xgrid, ygrid, lwd = 3, col = "skyblue4", lty = 5)

# abline(0, 1, col = "red", lwd = 2)

# Report correlation
correlation <- cor(d_xy_vec, d_pca_vec, method = "spearman")
correlation # 0.3194651


# --- Fit K-means using chosen k (e.g., k=4 based on elbow plot) ----
set.seed(123)
k_optimal <- 6
kmeans_result <- kmeans(env_data_only, centers = k_optimal, nstart = 10)

# --- Assign each point to a cluster ----
cluster_assignments <- kmeans_result$cluster

# --- Initialize cluster list, enforcing class diversity ---
min_class_count <- 4
cluster_data <- list()
cluster_count <- max(cluster_assignments)
has_class_variable <- FALSE   # <-- CHANGE HERE depending on your data

for (i in 1:cluster_count) {
  cluster_i <- env_data_only[cluster_assignments == i, ]
  # --- Case 1: You DO have a categorical class variable ---
  if (has_class_variable) {
    class_distribution <- table(cluster_i$occurrenceStatus)
  
    # Keep if both classes are sufficiently represented
    if (sum(class_distribution >= min_class_count) >= 2) {
      cluster_data[[i]] <- cluster_i
    } else {
      # Try merging with a cluster that has the opposite class
      opposite_class <- ifelse(names(class_distribution)[which.max(class_distribution)] == "1", 0, 1)
      for (j in setdiff(1:cluster_count, i)) {
        cluster_j <- env_data_only[cluster_assignments == j, ]
        if (opposite_class %in% cluster_j$occurrenceStatus) {
          cluster_data[[i]] <- rbind(cluster_i, cluster_j)
          cluster_data[[j]] <- NULL
          break
        }
      }
    }
  } else {
    # --- Case 2: Numeric-only data (no categorical class variable) ---
    # Simply store clusters without class-based constraints
    cluster_data[[i]] <- cluster_i
  }
}

# --- Clean up any NULLs due to merges ---
cluster_data <- cluster_data[!sapply(cluster_data, is.null)]
cluster_count <- length(cluster_data)

cat("Final number of usable clusters: ", cluster_count, "\n")


# --- Export cluster assignments  ----
# Create a data frame with replicates and their cluster assignments
cluster_df <- data.frame(
  replicates = pred$replicates,
  cluster = cluster_assignments
)

# Save to CSV
write.csv(cluster_df, "./data/processed_data/Cross-val/cluster_assignments_k6_envCV_T0.0.csv", row.names = FALSE)

# --- Print number of points in each fold ----
# For each fold, print nb of points in the fold and in the remaining data
total_points <- nrow(env_data_only)
fold_points <- numeric(cluster_count)
for (i in 1:cluster_count) {
  fold_points <- nrow(cluster_data[[i]])
  remaining_points <- total_points - fold_points
  cat("Fold ", i, ": test : ", fold_points, ", train : ", remaining_points, "\n", sep = "")
}



# Fold 1: test : 79, train : 558
# Fold 2: test : 137, train : 500
# Fold 3: test : 39, train : 598
# Fold 4: test : 74, train : 563
# Fold 5: test : 123, train : 514
# Fold 6: test : 185, train : 452




