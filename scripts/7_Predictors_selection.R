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

#------------- Load and prep data ------------------
# Load transformed predictors v.0.0
pred <- st_read("./data/processed_data/predictors/predictors_tr_v1.0.gpkg")

#--------------- PREDICTORS SELECTION v.0.0 ---------------- 

colnames(pred)



sel <- pred %>%
  dplyr::select(c(
  ))


# 1. 
# sal, temp, chl : - 1 month mean
# habitat : main habitat
# mean bathy
# mean gravity 
# mean distance to shore 
# mean distance to canyon
# mean distance to port 
# mean distance to reserve
# reserve in / out 
# mean terrain indices

# 3. compute correlation matrix









#---------------- PRELIMINAR PREDICTORS SELECTION / CLEANING ----------------

# Remove range ----
# Remove predictors containing "range"
range_cols <- grepl(colnames(pred), "_range")

sel <- pred %>%
  dplyr::select(which(!range_cols))
















#--------------- PREDICTORS SELECTION PROCEDURE ----------------