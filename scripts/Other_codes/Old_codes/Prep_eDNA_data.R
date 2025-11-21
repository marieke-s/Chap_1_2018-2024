#------------- Description ---------------------
# Purpose: 
# This script aims to 

# The expected result is 

# Data source : 

# Author: Marieke Schultz

# Contributors : 

# Date script created: 23/10/2025

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
library(tidyr)
library(stringr)



# Load functions
source("./utils/Fct_Data-Prep.R")


#------------- Load and prep data ------------------
# Occurences
occ <- read.csv("./data/raw_data/eDNA/data_MED_teleo_pres_1824_V1.csv", sep = ";")

# Rename code_spygen 
# Rename pomatomus.saltatrix to Pomatomus.saltatrix
occ <- occ %>%
  rename(spygen_code = code_spygen) %>%
  rename(Pomatomus.saltatrix = pomatomus.saltatrix)

# Import PCR matrix 
pcr <- read.csv("./data/raw_data/eDNA/data_MED_teleo_rep_1824_V1.csv", sep = ";")

# Rename code_spygen
# Remove "variable" column
pcr <- pcr %>%
  rename(spygen_code = code_spygen) %>%
  dplyr::select(-variable)



# Metadata
mtdt_3 <- sf::st_read("./data/processed_data/eDNA/mtdt_3.gpkg")

# Replace replicates = "no" by "spygen_code"
mtdt_3$replicates[mtdt_3$replicates == "no"] <- mtdt_3$spygen_code[mtdt_3$replicates == "no"]

# When pool =/= "no" --> spygen_code = pool
mtdt_3$spygen_code[mtdt_3$pool != "no"] <- mtdt_3$pool[mtdt_3$pool != "no"]

#------------- Check and clean species assignations  ----------------
# Load species lists from different sources ----

sp_edna <- colnames(occ[-1])

outcelia <- read.csv("./data/raw_data/eDNA/out_of_range_species1.2.csv")
sp_celia <- outcelia$Species

sp_jeanne <- c(
  "Abramis brama",
  "Alburnus alburnus",
  "Barbatula quignardi",
  "Barbus barbus",
  "Gambusia holbrooki",
  "Gobio gobio",
  "Leucaspius delineatus",
  "Lepomis gibbosus",
  "Neogobius fluviatilis",
  "Neogobius melanostomus",
  "Notacanthus chemnitzii",
  "Ophisurus macrorhynchos",
  "Perca fluviatilis",
  "Pholis gunnellus",
  "Phoxinus phoxinus",
  "Planiliza macrolepis",
  "Rhodeus amarus",
  "Rutilus rutilus",
  "Scardinius erythrophthalmus",
  "Schedophilus velaini",
  "Silurus glanis",
  "Squalius cephalus",
  "Tinca tinca",
  "Acanthurus coeruleus",
  "Anguilla marmorata",
  "Artediellus neyelovi",
  "Aspidophoroides monopterygius",
  "Carangoides plagiotaenia",
  "Caranx ignobilis",
  "Cololabis saira",
  "Cyclopterus lumpus",
  "Electrophorus electricus",
  "Encrasicholina punctifer",
  "Epinephelus maculatus",
  "Gymnocanthus tricuspis",
  "Gymnosarda unicolor",
  "Hyperoplus lanceolatus",
  "Istiophorus platypterus",
  "Krobia guianensis",
  "Leporinus affinis",
  "Lutjanus bohar",
  "Lutjanus rivulatus",
  "Mallotus villosus",
  "Nansenia boreacrassicauda",
  "Oncorhynchus clarkii",
  "Oncorhynchus keta",
  "Pimelodella cristata",
  "Reinhardtius hippoglossoides",
  "Phoxinus septimaniae",
  "Babka gymnotrachelus"
)

# Modif Alicia ----

# # Manually check synonyms and find the species
# species[species == "Chelon auratus"] <- "Liza aurata"
# rownames(adne)[rownames(adne) == "Chelon_auratus"] <- "Liza_aurata"
# 
# species[species == "Chelon ramada"] <- "Liza ramada"
# rownames(adne)[rownames(adne) == "Chelon_ramada"] <- "Liza_ramada"
# 
# species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
# rownames(adne)[rownames(adne) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"
# 
# species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
# rownames(adne)[rownames(adne) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"
# 
# species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
# rownames(adne)[rownames(adne) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"



# Check differences between species lists ----
# Print species in sp_jeanne not in sp_celia
setdiff(sp_jeanne, sp_celia)

# Print species in sp_celia not in sp_jeanne
setdiff(sp_celia, sp_jeanne)
length(setdiff(sp_celia, sp_jeanne)) # 14 spp.

# Species to remove ----
# Count species to remove
outcelia %>% filter(Decision=="OUT") %>% dim() # 56 species to remove

# List species to remove
sp_to_remove <- outcelia %>%
  filter(Decision == "OUT") %>%
  pull(Species)

# Replace " " by "." to match occ column names
sp_to_remove <- gsub(" ", ".", sp_to_remove)

length(setdiff(sp_to_remove, colnames(occ[-1]))) # 55/56 species to remove are not in occ V1. --> only 1 species will need to be removed from occ V1.
intersect(sp_to_remove, colnames(occ[-1])) # "Cyprinus.carpio" (freshwater species)

# Remove species from occ and pcr 
occ <- occ %>%
  dplyr::select(-any_of(sp_to_remove))

pcr <- pcr %>%
  dplyr::select(-any_of(sp_to_remove))





# Clean 
rm(sp_edna, sp_celia, sp_jeanne, outcelia, sp_to_remove)



# [DO NOT RUN] Remove rare species ----
# # List rare species (i.e., species under a certain occurrence threshold)
# rare <- list_rare_sp(occ[,-1], sum = 5, exact = FALSE)  # sum = X sets the number of occurrences under which a species is considered rare
# print(length(rare)) # 60 spp < 5 occurences 
# # 81 spp < 10 occurences
# # 102 spp < 20 occurences
# 
# 
# # Remove rare species with less than 5 occurences
# # Explanation : 06/01/2025 : We start by removing species with less than 5 occurences : Loïc et Lola = ont retire les espèces avec moins de 5 occurences, Lola = test sur 4 et 10 :  ne change rien. Célia = retire les espèces avec moins de 10 occurences.
# 
# occ <- occ[, !colnames(occ) %in% rare] # 60 spp with < 5 occurences : 142 spp --> 82 spp
# rm(rare)
