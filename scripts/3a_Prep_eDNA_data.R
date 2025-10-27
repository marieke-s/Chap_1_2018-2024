#------------ TODO ---------------------
# Replace matrix with V2 version of the matrix (without species correction)


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







#------------- Pool occ by replicates ------------------
# CHECKS : Compare spygen_codes in occ and mtdt ---- ----
length(occ) - length(mtdt_3) # 199 rows more in occ than in mtdt_3

# spygen_codes in occ not in mtdt_3
setdiff(occ$spygen_code, mtdt_3$spygen_code)
length(setdiff(occ$spygen_code, mtdt_3$spygen_code)) # 447

# spygen_codes in mtdt_3 not in occ
setdiff(mtdt_3$spygen_code, occ$spygen_code)
length(setdiff(mtdt_3$spygen_code, occ$spygen_code)) # 40

# Check remaining missing spygen_codes 
missing_spygen <- mtdt_3 %>%
  dplyr::filter(!spygen_code %in% occ$spygen_code)

unique(missing_spygen$project) # All samples come from the project "Ange2mer" --> need to be removed because they were not analysed with Teleo markers (they were specifically analysed for Squatina squatina).


# Remove Ange2Mer spygen_codes from mtdt_3 ----
mtdt_3 <- mtdt_3 %>%
  dplyr::filter(project != "Ange2Mer")

rm(missing_spygen)

# CHECKS : Compare spygen_codes in pcr and mtdt ---- ----
length(pcr) - length(mtdt_3) # 199 rows more in pcr than in mtdt_3

# spygen_codes in pcr not in mtdt_3
setdiff(pcr$spygen_code, mtdt_3$spygen_code)
length(setdiff(pcr$spygen_code, mtdt_3$spygen_code)) # 447

# spygen_codes in mtdt_3 not in pcr
setdiff(mtdt_3$spygen_code, pcr$spygen_code)
length(setdiff(mtdt_3$spygen_code, pcr$spygen_code)) # 0

# Remove spygen_codes not in mtdt_3 from occ and pcr ----
occ <- occ %>%
  dplyr::filter(spygen_code %in% mtdt_3$spygen_code)

pcr <- pcr %>%
  dplyr::filter(spygen_code %in% mtdt_3$spygen_code)

# Pool per replicates ----
# Add replicates column to occ
occ <- occ %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")

pcr <- pcr %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")

# Identify species columns (excluding 'spygen_code', 'replict', 'geometry')
species_cols <- setdiff(colnames(occ), c("spygen_code", "replicates", "geom"))

# Nb of expected rows after pooling
length(unique(mtdt_3$replicates)) # 768 

# Presence/Absence Pooled -----
# Merge rows per replicates --> 0 or 1 
occ_pooled <- occ %>%
  group_by(replicates) %>%
  summarise(
    pooled_name = paste(sort(spygen_code), collapse = "_"),  # Combine spygen_code separated by "_"
    across(all_of(species_cols), ~ ifelse(any(. == 1), 1, 0)),  # Apply new merging rule
  ) %>%
  ungroup()

# Check nb of rows after pooling
nrow(occ_pooled)  # 768 rows after pooling
length(setdiff(unique(mtdt_3$replicates), occ_pooled$replicates)) # 0
length(setdiff(occ_pooled$replicates, unique(mtdt_3$replicates))) # 0 




#-------------- Incertitude matrix ----------------------
# Field replicate incertitude (occ) ----

# Count number of field replicates 
# Important methodological consideration : Pooled_samples are counted as 1 field replicate
occ$field_replicates <- sapply(occ$replicates, compute_field_replicates)



# Merge rows per replicates --> nb of "1" / nb of field replicates
occ_f_incertitude <- occ %>%
  group_by(replicates) %>%
  summarise(
    pooled_name = paste(sort(spygen_code), collapse = "_"),
    across(
      all_of(species_cols),
      ~ ifelse(
        first(field_replicates) > 0,
        sum(. == 1, na.rm = TRUE) / first(field_replicates),
        NA_real_
      )
    )
  ) %>%
  ungroup()


# PCR replicate incertitude  (pcr) ----

# Identify species columns 
species_cols <- setdiff(colnames(pcr), c("spygen_code", "replicates", "geom"))

# Count number of pcr replicates 
pcr$pcr_replicates <- sapply(pcr$replicates, compute_pcr_replicates)

# Merge rows per replicates --> nb >0 / nb of pcr replicates

occ_pcr_incertitude <- pcr %>%
  group_by(replicates) %>%
  summarise(
    pooled_name = paste(sort(spygen_code), collapse = "_"),
    across(
      all_of(species_cols),
      ~ ifelse(
        first(pcr_replicates) > 0,
        sum(. > 0, na.rm = TRUE) / first(pcr_replicates),
        NA_real_
      )
    )
  ) %>%
  ungroup()




# CHECKS pcr values -----

# compute max per columns of df
sort(sapply(pcr[ , species_cols], max, na.rm = TRUE))
max(sapply(pcr[ , species_cols], max, na.rm = TRUE))

# PCR values above 12 ? 
t <- pcr %>% filter(Gobius.xanthocephalus == '48' ) %>% print() # 48 PCR replicates 
t$spygen_code # "SPY210641" --> maximum should be 12

# make an hist of pcr across all columns 
hist(as.numeric(unlist(pcr[ , species_cols])), breaks = 50, main = "Histogram of PCR values", xlab = "PCR values", col = "lightblue")
boxplot(as.numeric(unlist(pcr[ , species_cols])), main = "Boxplot of PCR values", ylab = "PCR values")
min(as.numeric(unlist(pcr[ , species_cols])))
max(as.numeric(unlist(pcr[ , species_cols])))
median(as.numeric(unlist(pcr[ , species_cols])))
mean(as.numeric(unlist(pcr[ , species_cols])))

# count nb of as.numeric(unlist(pcr[ , species_cols])) > 12
sum(as.numeric(unlist(pcr[ , species_cols])) > 12, na.rm = TRUE) # 579 values above 12 --> WHY ? 

# Select rows containing species_cols values > 12
pcr_above_12 <- pcr %>%
  filter(if_any(all_of(species_cols), ~ . > 12))


# in mtdt_3 select spygen_code in pcr_above_12
mtdt_3_pcr_above_12 <- mtdt_3 %>%
  filter(spygen_code %in% pcr_above_12$spygen_code)

# Export mtdt_3_pcr_above_12 as csv
# Drop geometry
mtdt_3_pcr_above_12 <- st_drop_geometry(mtdt_3_pcr_above_12)
write.csv(mtdt_3_pcr_above_12, "./data/processed_data/eDNA/mtdt_3_pcr_above_12.csv", row.names = FALSE)





rm(species_cols)


# ????? Remove replicates column ----
occ_pooled <- occ_pooled %>% dplyr::select(-replicates)
occ_incertitude <- occ_incertitude %>% dplyr::select(-replicates)


















#------------- UNDETECTED AND RARE SPECIES --------------------------------------

# CHECKS : Plot occurences and rarity ----
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







# RUN : Remove undetected species ----
# List undetected species (i.e., species with no occurrence)
undetected <- list_rare_sp(occ[,-1], sum = 0, exact = TRUE)  # Using occ[,-1] to exclude spygen_code column
print(length(undetected))  

# Remove undetected species
occ <- occ[, !colnames(occ) %in% undetected] # 96 undetected spp : 238 --> 142 spp # with Litto3D-removed points : 94 undetected spp : 238 --> 144 spp

# Cleanup
rm(undetected)



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


