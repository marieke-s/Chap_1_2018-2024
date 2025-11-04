#------------- TODO ---------------------
# Redo code with new spygen data to be sent with reassigned samples deleting previous version to ensure 12 PCR per samples max. (Laure email 28/10/25)
# Handle species complex 
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
occ <- read.csv("./data/processed_data/eDNA/data_MED_teleo_pres_1824_V1.csv", sep = ";")

# Rename spygen_code 
occ <- occ %>%
  rename(spygen_code = sample)

# Import PCR matrix 
pcr <- read.csv("./data/processed_data/eDNA/data_MED_teleo_nb_rep_1824_V1.csv", sep = ";")

# Rename spygen_code
pcr <- pcr %>%
  rename(spygen_code = sample)

# Metadata
mtdt_3 <- sf::st_read("./data/processed_data/Mtdt/mtdt_3.gpkg")

# Replace replicates = "no" by "spygen_code"
mtdt_3$replicates[mtdt_3$replicates == "no"] <- mtdt_3$spygen_code[mtdt_3$replicates == "no"]

# Reoder pool values in increasing order 
mtdt_3$pool <- sapply(mtdt_3$pool, function(x) {
  # If not a pooled sample (e.g. "no"), return as-is
  if (x == "no") return(x)
  
  # Split into two SPY codes
  codes <- unlist(strsplit(x, "_"))
  
  # Extract numeric parts for ordering
  nums <- as.numeric(sub("SPY", "", codes))
  
  # Reorder by number and rebuild
  ordered_codes <- codes[order(nums)]
  paste(ordered_codes, collapse = "_")
})


# When pool =/= "no" --> spygen_code = pool
mtdt_3$spygen_code[mtdt_3$pool != "no"] <- mtdt_3$pool[mtdt_3$pool != "no"]






#------------- PREP DATA ----------------
# CHECKS : Compare spygen_codes in occ and mtdt ---- ----
length(occ) - length(mtdt_3) # 222 rows more in occ than in mtdt_3

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
length(pcr) - length(mtdt_3) # 222 rows more in pcr than in mtdt_3

# spygen_codes in pcr not in mtdt_3
setdiff(pcr$spygen_code, mtdt_3$spygen_code)
length(setdiff(pcr$spygen_code, mtdt_3$spygen_code)) # 447

# spygen_codes in mtdt_3 not in pcr
setdiff(mtdt_3$spygen_code, pcr$spygen_code)
length(setdiff(mtdt_3$spygen_code, pcr$spygen_code)) # 0

# Subset occ and pcr using mtdt_3 ----
occ <- occ %>%
  dplyr::filter(spygen_code %in% mtdt_3$spygen_code)

pcr <- pcr %>%
  dplyr::filter(spygen_code %in% mtdt_3$spygen_code)





# Add replicates column to occ and pcr ----
occ <- occ %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")

pcr <- pcr %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")



#-------------- SAMPLES WITH NO SPECIES DETECTED ----------------------
# Print samples with no species detected ----
species_cols <- setdiff(colnames(occ), c("spygen_code", "replicates", "geom"))
occ[species_cols] <- lapply(occ[species_cols], as.numeric)  # Ensure species columns are numeric

# Keep only samples with no species detected
d <- occ %>% dplyr:: filter(rowSums(across(all_of(species_cols)), na.rm = TRUE) == 0)

length(unique(d$spygen_code)) # 25 samples with 0 species detected 
print(unique(d$spygen_code))

no_sp <- unique(d$spygen_code)
no_sp <- c(
  "SPY182541", "SPY201282", "SPY201274", "SPY211141_SPY211142", 
  "SPY211169_SPY211181", "SPY211168_SPY211173", "SPY211172_SPY211177", 
  "SPY211178_SPY211184", "SPY211047", "SPY211285", "SPY231906", 
  "SPY211061", "SPY222967", "SPY210700", "SPY231871", "SPY233819", 
  "SPY233363", "SPY232133", "SPY232282", "SPY232289", "SPY233720", 
  "SPY233722", "SPY2401455", "SPY2401628", "SPY2401629"
)

# Subset mtdt_3 to check these samples
mtdt_empty_samples <- mtdt_3 %>%
  dplyr::filter(spygen_code %in% d$spygen_code) 

# SPY211285 : Comment : "No_filtration - Blanc terrain/test contamination" / component = "blank"
# 4 DeepHeart 2021-10-19 : field error pump was not properly set 

#---> Mail to Laure to check if this is normal for the other samples (03/11/2025)


# Remove samples with 0 species detected /!\ IF DOING SO WE NEED TO MODIFY REPLICATES AND SAMPLING EFFORT COLS ??----

rm(d, mtdt_empty_samples)
#------------- POOL OCC BY REPLICATES ------------------
# Presence/Absence Pooled -----
# Merge rows per replicates --> 0 or 1 
species_cols <- setdiff(colnames(occ), c("spygen_code", "replicates", "geom"))

occ_pooled <- occ %>%
  group_by(replicates) %>%
  summarise(
    pooled_name = paste(sort(spygen_code), collapse = "_"),  # Combine spygen_code separated by "_"
    across(all_of(species_cols), ~ ifelse(any(. == 1), 1, 0)),  # Apply new merging rule
  ) %>%
  ungroup()

# Check nb of rows after pooling
message("Number of rows after pooling: ", nrow(occ_pooled), " (Expected: ", length(unique(mtdt_3$replicates)), ")")

length(setdiff(unique(mtdt_3$replicates), occ_pooled$replicates)) # 0
length(setdiff(occ_pooled$replicates, unique(mtdt_3$replicates))) # 0 




#-------------- INCERTITUDES MATRICES ----------------------
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
species_cols <- setdiff(colnames(pcr), c("spygen_code", "replicates", "geom", "pcr_replicates"))

# Count number of pcr replicates 
pcr$pcr_replicates <- sapply(pcr$replicates, compute_pcr_replicates)

# Merge rows per replicates --> nb > 0 / nb of pcr replicates

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









# [ !!! PCR ISSUE !!! ] -----

#------------- CLEANUP ----------------------
# Clean pooled occurences 
occ_pooled <- occ_pooled %>%
  dplyr::select(-pooled_name)
occ_f_incertitude <- occ_f_incertitude %>%
  dplyr::select(-pooled_name)
occ_pcr_incertitude <- occ_pcr_incertitude %>%
  dplyr::select(-pooled_name)



#------------- UNDETECTED AND RARE SPECIES --------------------------------------
# Plot rare species frequencies ------
res <- plot_species_rarity(df = occ_pooled,
                           threshold = 10, 
                           col_to_del = c("replicates"),
                           plot_type = "both")
res$plots$bar
res$plots$hist








# Remove undetected species ----
# List undetected species (i.e., species with no occurrence)
undetected <- list_rare_sp(occ_pooled[,-1], sum = 0, exact = TRUE)  # Using occ[,-1] to exclude spygen_code column
print(length(undetected))  # 20 undetected species
print(undetected) #  6 chondrychtyens + 14 osteichthyens

# Remove undetected species
occ_pooled <- occ_pooled[, !colnames(occ_pooled) %in% undetected] # 242 species left
occ_f_incertitude <- occ_f_incertitude[, !colnames(occ_f_incertitude) %in% undetected]
occ_pcr_incertitude <- occ_pcr_incertitude[, !colnames(occ_pcr_incertitude) %in% undetected]

# Cleanup
rm(undetected)



# Clean
rm(res)
#-------------- REPLICATES WITH NO SPECIES DETECTED ----------------------
# Print replicates with no species detected ----
t <- occ_pooled %>% dplyr::filter(rowSums(across(where(is.numeric))) == 0)

print(t$replicates) # 6 replicates with 0 species detected :
# "SPY210700"                               "SPY211141_SPY211142"                     "SPY211168_SPY211173/SPY211177_SPY211172"
# "SPY211285"                               "SPY232282/SPY232289"                     "SPY2401628/SPY2401629"      

# Subset mtdt_3 to check these samples
mtdt_empty_samples <- mtdt_3 %>%
  dplyr::filter(replicates %in% t$replicates) 

# SPY211285 : Comment : "No_filtration - Blanc terrain/test contamination"
# 4 DeepHeart 2021-10-19 : field error pump was not properly set 

#---> Mail to Laure to check if this is normal for the other samples (03/11/2025)


# Remove replicates with 0 species detected ----
occ_pooled <- occ_pooled %>% dplyr::filter(rowSums(across(where(is.numeric))) != 0) 
occ_pcr_incertitude <- occ_pcr_incertitude %>% dplyr::filter(rowSums(across(where(is.numeric))) != 0)
occ_f_incertitude <- occ_f_incertitude %>% dplyr::filter(rowSums(across(where(is.numeric))) != 0)

rm(t, mtdt_empty_samples)

#------------- EXPORT POOLED DATA ----------------------
write.csv(occ_pooled, "./data/processed_data/eDNA/occ_pooled.csv", row.names = FALSE)
write.csv(occ_f_incertitude, "./data/processed_data/eDNA/occ_f_incertitude_pooled.csv", row.names = FALSE)
write.csv(occ_pcr_incertitude, "./data/processed_data/eDNA/occ_pcr_incertitude_pooled.csv", row.names = FALSE)


occ_pooled %>% dplyr::select("Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus") %>% filter("Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus" == 1) %>% dim()

sort(colnames(occ_pooled))
