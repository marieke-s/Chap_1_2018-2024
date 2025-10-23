#------------- Description ---------------------
# Purpose: 
# This script aims to 

# The expected result is 

# Data source: 

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


#------------- Load data ------------------
occ <- read.csv("./data/raw_data/eDNA/data_MED_teleo_pres_1824_V1.csv", sep = ";")




#------------- Check and clean species assignations  ----------------
sp_edna <- colnames(occ[-1])

sp_celia <- read.csv("./data/raw_data/eDNA/out_of_range_species1.1.csv")
sp_celia <- sp_celia$Species

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

# Modif Alicia

# Manually check synonyms and find the species
species[species == "Chelon auratus"] <- "Liza aurata"
rownames(adne)[rownames(adne) == "Chelon_auratus"] <- "Liza_aurata"

species[species == "Chelon ramada"] <- "Liza ramada"
rownames(adne)[rownames(adne) == "Chelon_ramada"] <- "Liza_ramada"

species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
rownames(adne)[rownames(adne) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"

species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
rownames(adne)[rownames(adne) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"

species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
rownames(adne)[rownames(adne) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"



# Print species in sp_jeanne not in sp_celia
setdiff(sp_jeanne, sp_celia)


# Print species in sp_celia not in sp_jeanne
setdiff(sp_celia, sp_jeanne)
length(setdiff(sp_celia, sp_jeanne))

# Clean 
rm(sp_edna, sp_celia, sp_jeanne)





#---------------------------------------------------------------------------- 1. UNDETECTED AND RARE SPECIES --------------------------------------

#---- CHECKS : Plot occurences and rarity ----
# Explanation : This part plots the number of occurence per species distribution. It highlights the species that have few occurences (under a chosen rarity threshold).

### 1. Rare species plot (Occurence / species)

# Define your threshold for highlighting species
threshold <- 5  # Set the threshold for highlighting (e.g., < 10 occurrences)

# Data preparation
plot_data <- occ %>%
  dplyr::select(-1) %>%  # Exclude the first column
  summarise(across(everything(), sum, na.rm = TRUE)) %>%  # Sum occurrences for each species
  tidyr::pivot_longer(cols = everything(), names_to = "species", values_to = "tot") %>%  # Reshape data
  mutate(
    rank = rank(-tot),  # Rank in descending order
    above_mean = ifelse(tot > mean(tot, na.rm = TRUE), 1, 0),  # Indicator for above mean
    below_mean = ifelse(tot < mean(tot, na.rm = TRUE), 1, 0)   # Indicator for below mean
  )

# Store variables for plotting 
mean_occ <- mean(plot_data$tot, na.rm = TRUE) # mean occurrence
tot <- plot_data$tot  # Total occurrences
percentage_above_mean <- (sum(plot_data$above_mean) / nrow(plot_data)) * 100
percentage_below_mean <- (sum(plot_data$below_mean) / nrow(plot_data)) * 100

# Calculate species below threshold
highlighted_species <- plot_data %>%
  filter(tot < threshold)

num_highlighted <- nrow(highlighted_species)  # Count species below threshold
percentage_highlighted <- (num_highlighted / nrow(plot_data)) * 100  # Percentage

# Plotting using ggplot2
# Data
ggplot(plot_data, aes(x = reorder(species, -tot), y = tot)) +
  geom_bar(stat = "identity", aes(fill = ifelse(tot < threshold, "Rare", "Common"))) +  # Use 'Highlighted' and 'Regular' for fill
  geom_hline(yintercept = threshold, linetype = "dashed", color = "darkslategrey", size = 0.5) +  # Add horizontal line for threshold
  geom_hline(yintercept = mean_occ, linetype = "dashed", color = "firebrick4", size = 0.5) +  # Add horizontal line for mean occurence
  
  # Annotations
  geom_text(x = max(seq_along(tot)) - 34, y = mean_occ,
            label = paste("Rarity threshold:", round(threshold, 0)),
            vjust = 2, hjust = 0, color = "darkslategrey", size = 4) +
  geom_text(x = max(seq_along(tot)) - 34, y = mean_occ,
            label = paste("Mean occurrence:", round(mean_occ, 0)),
            vjust = -1, hjust = 0, color = "firebrick4", size = 4) +
  # geom_text(aes(label = paste(round(percentage_above_mean, 1), "% above mean")),
  #           x = max(seq_along(tot)) - 33, y = mean_occ - 8,
  #           vjust = -5, hjust = 0, color = "firebrick", size = 3.5) +
  # geom_text(aes(label = paste(round(percentage_below_mean, 1), "% below mean")),
  #           x = max(seq_along(tot)) - 33, y = mean_occ - 10,
  #           vjust = -7, hjust = 0, color = "firebrick", size = 3.5) +
  scale_fill_manual(values = c("Rare" = "darkslategrey", "Common" = "darkslategray3"), name = "Species Type") +  
  theme_test() +
  
  # Custom fill with legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(color = "black", size = 10, hjust = 0.5),  # Center the caption
        legend.position = "bottom") +  # Position legend at the bottom
  labs(title = "Total Occurrences and Ranking for Each Species",
       x = "Species",
       y = "Total Occurrences",
       caption = paste("Rare species are species with no more than", threshold, "occurences. Total:", num_highlighted, "/",nrow(plot_data), "species (", round(percentage_highlighted, 1), "%)."))


# Cleanup
rm(plot_data, mean_occ, tot, percentage_above_mean, percentage_below_mean, highlighted_species, num_highlighted, percentage_highlighted, threshold)
#dev.off()


### 2. Occurence frequency distribution
hist(plot_data$tot, breaks = 10, col = "lightblue", xlab = "Total Occurrences", ylab = "Frequency", main = "Distribution of Total Occurrences per Species")







#---- RUN : Remove undetected species ----
# List undetected species (i.e., species with no occurrence)
undetected <- list_rare_sp(occ[,-1], sum = 0, exact = TRUE)  # Using occ[,-1] to exclude spygen_code column
print(length(undetected))  

# Remove undetected species
occ <- occ[, !colnames(occ) %in% undetected] # 96 undetected spp : 238 --> 142 spp # with Litto3D-removed points : 94 undetected spp : 238 --> 144 spp

# Cleanup
rm(undetected)



#---- [DO NOT RUN] Remove rare species ----
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


#---- Save results ----
# Explanation of file naming : P0 for Processing 0 (from data processing computed in analyses/data_prep/4_Data_Processing.R) _RX for Rarity set to threshold 0 (species with less than X occurences of removed).

# Merge back occ to tot
# remove occ columns from tot
result <- identify_columns(tot) 
occ1 <- colnames(tot)[result$response_columns]
# Remove occ1 columns from tot
tot <- tot[, -which(colnames(tot) %in% occ1)]

# Merge occ to tot
tot <- merge(tot, occ, by = "spygen_code", all.x = TRUE)

colnames(tot)

write_rds(occ, "./data/processed_data/data_prep/Med_TOT_2023_P0_R0.rds")

