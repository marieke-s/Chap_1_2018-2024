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
  select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 1: Spherical")

range_data_2 <- results_df_2 %>% 
  select(variable, range, moran_i, moran_p) %>% 
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
  select(variable, range, moran_i, moran_p) %>% 
  mutate(method = "Method 1: Spherical")

range_data_2 <- results_df_2 %>% 
  select(variable, range, moran_i, moran_p) %>% 
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



