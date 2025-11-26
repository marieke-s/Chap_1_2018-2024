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











