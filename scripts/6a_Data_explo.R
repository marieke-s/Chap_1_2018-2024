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
library(lubridate)
library(terra)
library(stringr)
library(pMEM)



# Load functions
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








# Temporal variability ----

# Data prep ----
# Compute year, month and season from date 
mtdt <- mtdt %>%
  mutate(
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    ))

# Replace month number by month name
mtdt <- mtdt %>%
  mutate(
    month = factor(month.abb[month], levels = month.abb)
  )

# Cols to plot 
cat_cols <- c("year", "month", "season")



# Map + hist season | month | year -----
map_categorical_plots(mtdt,
                      cols_to_plot = cat_cols,
                      version = "mtdt_7_temporal_var",
                      output_directory = "./figures/Mtdt/Map_Hist")








# Sampling effort variability ----

# Add area_km2 and dist_seabed_depthsampling to mtdt
mtdt <- mtdt %>%
  left_join(
    pred_raw %>%
      st_drop_geometry() %>%
      dplyr::select(replicates, area_km2, dist_seabed_depthsampling)
    , by = "replicates"
  )


# Map and hist of sampling effort 


#------------- Explore raw predictors ---------------------------------------------------------------------------------------------------------------
# Check NAs -----
sapply(pred_raw, function(x) sum(is.na(x))) 
# in predictors_raw_v1.2 : 5 NAs for habitat predictors

# Map + hist + summary : numerical data -----

cols_to_plot <- pred_raw[2:151] %>%
  st_drop_geometry() %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

map_index_plots(df = pred_raw, version = "predictors_raw_v1.2", output_directory = "./figures/Predictors/Map_Hist", cols_to_plot = cols_to_plot)


rm(cols_to_plot)



# Outliers ----











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

