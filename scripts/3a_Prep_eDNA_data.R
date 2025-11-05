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

# Reorder spygen_codes
occ$spygen_code <- sapply(occ$spygen_code, reorder_spygen_codes)




# PCR matrix 
pcr <- read.csv("./data/processed_data/eDNA/data_MED_teleo_nb_rep_1824_V1.csv", sep = ";")

# Rename spygen_code
pcr <- pcr %>%
  rename(spygen_code = sample)

# Reorder spygen_codes
pcr$spygen_code <- sapply(pcr$spygen_code, reorder_spygen_codes)




# Metadata
mtdt_3 <- sf::st_read("./data/processed_data/Mtdt/mtdt_3.gpkg")

# Replace replicates = "no" by "spygen_code"
mtdt_3$replicates[mtdt_3$replicates == "no"] <- mtdt_3$spygen_code[mtdt_3$replicates == "no"]

# Reorder pool
mtdt_3$pool <- sapply(mtdt_3$pool, reorder_spygen_codes)

# When pool =/= "no" --> spygen_code = pool
mtdt_3$spygen_code[mtdt_3$pool != "no"] <- mtdt_3$pool[mtdt_3$pool != "no"]

# Reorder pool within replicates
# Function to reorder SPY codes within "_" groups, keeping "/" structure
reorder_replicates <- function(x) {
  # Split the string by "/" to handle each group separately
  groups <- unlist(strsplit(x, "/"))
  
  # Process each group independently
  reordered_groups <- sapply(groups, function(g) {
    # Only reorder if multiple SPY codes separated by "_"
    if (grepl("_", g)) {
      codes <- unlist(strsplit(g, "_"))
      nums <- as.numeric(sub("SPY", "", codes))
      codes <- codes[order(nums)]
      paste(codes, collapse = "_")
    } else {
      g
    }
  })
  
  # Recombine the groups using "/"
  paste(reordered_groups, collapse = "/")
}

# Apply to your dataframe column
mtdt_3$replicates <- sapply(mtdt_3$replicates, reorder_replicates)






#--------------- HANDLE SPECIES COMPLEXES -----
# Check individual species of the complexes -----
## 1. Extract complex ie. colnames containing a "."
complex_species <- colnames(occ)[grepl("\\.", colnames(occ))]

## 2. Separate the species in the complexes and store in a list
separated_species <- unique(unlist(strsplit(complex_species, "\\.")))

## 3. Check which of these species are also present individually (ie not inside complexes)
separated_species[separated_species %in% colnames(occ)] # "Centrolabrus melanocercus" "Umbrina cirrosa" "Spicara smaris"   





# Resolve complexes -----
occ <- occ %>%
  
  
  ##--- 1. Species complex to delete :
  dplyr::select(-c("Coptodon_rendalli.Oreochromis_niloticus", # Both african fish that can be introduced in the Med. Common subfamily : Pseudocrenilabrinae. Oreochromis niloticus is in the traits db while Coptodon rendalli is not. --> DELETE 
  )) %>%

rename(
  
  
  ##--- 2. Species where geographical distributions help to choose :
  
  "Cheilopogon_heterurus" = "Cheilopogon_heterurus.Hirundichthys_speculiger", # Hirundichthys speculiger mostly found in tropical open water while Cheilopogon heterurus is a mediterranean fish.  (https://fishbase.se/summary/1029 and https://fishbase.se/summary/Hirundichthys-speculiger) --> KEEP Cheilopogon heterurus
  "Trachurus_mediterraneus" = "Trachurus_mediterraneus.Trachurus_trachurus", # This complex is present in 56% of the 2018-2024 samples detected on average on 4 PCR replicates (when present). Thus most probably Trachurus mediterraneus which is a least concerned fish that can widely spread in the Mediterranean sea (https://www.fishbase.se/summary/trachurus-mediterraneus) while Trachurus trachurus is a mostly atlantic species + is vulnerable (https://www.fishbase.se/summary/Trachurus-trachurus.html) --> KEEP Trachurus mediterraneus
  "Trisopterus_capelanus" = "Trisopterus_capelanus.Trisopterus_minutus", # Trisopterus minutus is not present in the Mediterranean sea (https://www.fishbase.se/summary/trisopterus-minutus) while Trisopterus capelanu is a mediterranean species (https://fishbase.se/summary/Trisopterus-capelanus.html) --> KEEP Trisopterus capelanus
  "Notoscopelus_elongatus" = "Notoscopelus_elongatus.Notoscopelus_kroyeri", # Notoscopelus kroyeri is endemic to the Atlantic sea (https://www.fishbase.se/summary/Notoscopelus-kroyeri) and Notoscopelus elongatus is found in the Mediterranean sea (https://fishbase.se/summary/841) --> KEEP Notoscopelus elongatus
  "Sphyraena_sphyraena" = "Sphyraena_chrysotaenia.Sphyraena_sphyraena", # Sphyraena chrysotaenia can be found in the Med as a Lessepsian migrant while Sphyraena sphyraena is commonly found in the Med (https://fishbase.se/Summary/SpeciesSummary.php?id=16905&lang=french and https://www.fishbase.se/summary/sphyraena-sphyraena) --> KEEP Sphyraena sphyraena
  
  
  ##--- 3. All species found in Med : keep closest common taxo
  
  "Parablennius_sp_1" = "Parablennius_tentacularis.Parablennius_zvonimiri", # Both can be found in our study area -->  Parablennius sp.
  "Parablennius_sp_2" = "Parablennius_incognitus.Parablennius_sanguinolentus", # Both can be found in our study area -->  Parablennius sp.
  "Gaidropsarus_sp" = "Gaidropsarus_biscayensis.Gaidropsarus_vulgaris", # Both can be found in our study area (https://www.fishbase.se/summary/1877 and https://doris.ffessm.fr/Especes/Gaidropsarus-vulgaris-Motelle-commune-2699) --->  Gaidropsarus sp.
  "Triglinae_sp" = "Chelidonichthys_lucerna.Lepidotrigla_dieuzeidei", # Both can be found in our study area (https://www.fishbase.se/summary/Chelidonichthys-lucerna.html and https://www.fishbase.se/summary/Lepidotrigla-dieuzeidei) -->  Triglinae sp.
  "Chelidonichthys_sp" = "Chelidonichthys_obscurus.Chelidonichthys_lastoviza", # Both can be found in our study area (fishbase) -->  Chelidonichthys sp.
  "Triglinae_sp_1" = "Eutrigla_gurnardus.Trigla_lyra", # Both can be found in our study area (fishbase) -->  Triglinae sp.
  "Sparidae_sp_2" = "Dentex_dentex.Pagrus_auriga.Pagrus_pagrus", # All 3 can be found in our study area (fishbase) -->  Sparidae sp.
  "Raja_sp" = "Raja_asterias.Raja_clavata.Raja_polystigma" # All 3 can be found in our study area (fishbase) -->  Raja sp.
) 




# Those where I don't know what to do : (check with reference data base on how they were detected)
# "Labrus merula.Labrus viridis", # also present in the complex Labrus merula.Labrus viridis.Centrolabrus melanocercus --> ??
# "Labrus merula.Labrus viridis.Centrolabrus melanocercus", # Centrolabrus melanocercus also present alone. Detected both alone and in the complex with MED_2025.  
# "Spicara flexuosum.Spicara smaris", # Spicara smaris also present alone. Detected both alone and in the complex with MED_2025.
# "Argyrosomus regius.Umbrina cirrosa", # Umbrina cirrosa also present alone. Detected alone with GenBank MIDORI_v264 and in the complex by MED_2025.


# For these complexes check if the species are detected in the same samples when alone and inside a complex
# Compare Labrus merula.Labrus viridis.Centrolabrus melanocercus and Centrolabrus melanocercus
occ %>%
  dplyr::filter(`Labrus merula.Labrus viridis.Centrolabrus melanocercus` == 1 & `Centrolabrus melanocercus` == 1) %>%
  dplyr::select(spygen_code, `Labrus merula.Labrus viridis.Centrolabrus melanocercus`, `Centrolabrus melanocercus`)

# Compare "Spicara flexuosum.Spicara smaris" and "Spicara smaris"
occ %>%
  dplyr::filter(`Spicara flexuosum.Spicara smaris` == 1 & `Spicara smaris` == 1) %>%
  dplyr::select(spygen_code, `Spicara flexuosum.Spicara smaris`, `Spicara smaris`)



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







#------------- HANDLE SAMPLES WITH NO SPECIES DETECTED ----------------------
# Print no detection samples  ----
species_cols <- setdiff(colnames(occ), c("spygen_code", "replicates", "geom"))
occ[species_cols] <- lapply(occ[species_cols], as.numeric)  # Ensure species columns are numeric

# Keep only samples with no species detected
nd <- occ %>% dplyr:: filter(rowSums(across(all_of(species_cols)), na.rm = TRUE) == 0)

length(unique(nd$spygen_code)) # 25 samples with 0 species detected 
print(unique(nd$spygen_code))


# Check no detection samples ----
# Check their mtdt
mtdt_empty_samples <- mtdt_3 %>%
  dplyr::filter(spygen_code %in% nd$spygen_code) 

# SPY211285 : Comment : "No_filtration - Blanc terrain/test contamination" / component = "blank"
# 4 DeepHeart 2021-10-19 : field error pump was not properly set 

# Compare with Laure's list of samples with no detection 
ndl <- read.csv("./data/raw_data/eDNA/NoDNA.csv", sep = ";") # includes 2018-2023 samples (no 2024)

length(setdiff(unique(nd$spygen_code), ndl$Sample)) # 10 samples are not in Laure's list

# Print 2018-2023 samples not in Laure's list
mtdt_empty_samples %>%
  dplyr::filter(spygen_code %in% setdiff(unique(nd$spygen_code), ndl$Sample)) %>%
  dplyr::filter(year(date) != 2024) %>%
  pull(spygen_code)

# 7 2023 samples : "SPY232133" "SPY232282" "SPY232289" "SPY233363" "SPY233720" "SPY233722" "SPY233819"
# Sent to Laure on 5/11/25 to check



# Check their buffer on QGIS
mtdt_2 <- sf::st_read("./data/processed_data/Mtdt/mtdt_2.gpkg")
mtdt_empty_samples_2 <- mtdt_2 %>%
  dplyr::filter(spygen_code %in% nd$spygen_code)
mtdt_empty_samples_3 <- mtdt_3 %>%
  dplyr::filter(spygen_code %in% nd$spygen_code)

sf::write_sf(mtdt_empty_samples_2, "./data/processed_data/Mtdt/mtdt_empty_samples_2.gpkg")
sf::write_sf(mtdt_empty_samples_3, "./data/processed_data/Mtdt/mtdt_empty_samples_3.gpkg")


## QGIS : Removing samples from replicates group will not to very slightly impact the buffer shape --> we can keep predictor's extractions already done on these shapes.





# Remove no detection samples----
occ <- occ %>%
  dplyr::filter(!spygen_code %in% nd$spygen_code)

pcr <- pcr %>%
  dplyr::filter(!spygen_code %in% nd$spygen_code)

mtdt_3 <- mtdt_3 %>%
  dplyr::filter(!spygen_code %in% nd$spygen_code)

# In mtdt_3, remove no detection samples from replicates ----
clean_replicates <- function(replicates, remove_codes) {
  # Split by "/" to separate replicate groups
  sapply(replicates, function(x) {
    groups <- unlist(strsplit(x, "/"))
    
    # Keep only those groups not fully matching any remove_codes
    kept_groups <- groups[!groups %in% remove_codes]
    
    # Within each group (like SPY211169_SPY211181), also remove
    # any subcodes that match any remove_code directly or as part of combined codes
    cleaned_groups <- sapply(kept_groups, function(g) {
      subcodes <- unlist(strsplit(g, "_"))
      # Remove any subcode present in remove_codes (unique codes)
      subcodes <- subcodes[!subcodes %in% remove_codes]
      paste(subcodes, collapse = "_")
    })
    
    # Combine cleaned groups back with "/"
    paste(cleaned_groups, collapse = "/")
  })
}

# Apply the function
mtdt_3 <- mtdt_3 %>% rename(replicates_old = replicates) 
mtdt_3$replicates <- clean_replicates(mtdt_3$replicates, nd$spygen_code)


# Cleanup
rm(mtdt_2, mtdt_empty_samples, mtdt_empty_samples_2, mtdt_empty_samples_3, nd, ndl, clean_replicates)




# Add replicates column to occ and pcr ----
occ <- occ %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")

pcr <- pcr %>%
  left_join(mtdt_3 %>% dplyr::select(c(spygen_code, replicates)), by = "spygen_code")


















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


colnames(occ)

# Check no detection replicates -----
occ_pooled %>%
  dplyr::filter(rowSums(across(all_of(species_cols)), na.rm = TRUE) == 0) %>%
  dim() # 0 replicates with no detection


#------------- COMPUTE SAMPLING EFFORTS IN MTDT ----------------------
# Nb of field replicates ----
# Compute nb of field replicates using the function compute_field_replicates
# Pooled_samples are counted as 1 field replicate
mtdt_3$field_replicates <- sapply(mtdt_3$replicates, compute_field_replicates)


# Nb of PCR replicates ----
# Compute nb of PCR replicates using the function compute_pcr_replicates
mtdt_3$PCR_replicates <- sapply(mtdt_3$replicates, compute_pcr_replicates)









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
occ <- occ %>% dplyr::select(-field_replicates)

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






#------------- EXPORT POOLED DATA ----------------------
write.csv(occ_pooled, "./data/processed_data/eDNA/occ_pooled.csv", row.names = FALSE)
write.csv(occ_f_incertitude, "./data/processed_data/eDNA/occ_f_incertitude_pooled.csv", row.names = FALSE)
write.csv(occ_pcr_incertitude, "./data/processed_data/eDNA/occ_pcr_incertitude_pooled.csv", row.names = FALSE)


occ_pooled %>% dplyr::select("Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus") %>% filter("Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus" == 1) %>% dim()

sort(colnames(occ_pooled))

#------------- EXPORT MODIFIED MTDT ----------------------
# In this script we made somes changes to the mtdt_3 that needs to be saved : removing of no detection samples, removing of Ange2Mer project, adding sampling effort columns. --> MTDT_6 (rows = sample)

# When pooled by replicates the mtdt should re-calculate the estimated_volume_total without the removed samples --> MTDT_7 (rows = replicates)



# MTDT_6 ----
mtdt_6 <- mtdt_3
sf::st_write(mtdt_6, "./data/processed_data/Mtdt/mtdt_6.gpkg", delete_dsn = TRUE)

# MTDT_7 ----
# Group by replicates ----
mtdt_7 <- mtdt_6 %>%
  group_by(replicates) %>%
  summarise(
    replicates_old = first(replicates_old),
    date = first(na.omit(date)),
    time_start = min(time_start, na.rm = TRUE),
    
    # Use first value if all non-NA values are identical; otherwise compute mean
    depth_sampling = if (n_distinct(na.omit(depth_sampling)) == 1) {
      first(na.omit(depth_sampling))
    } else {
      #print(cur_data_all())  # This prints the full data for the current group
      mean(depth_sampling, na.rm = TRUE)
    },
    
    # Use first value if all non-NA values are identical; otherwise take the first non-NA value
    depth_seafloor = if(n_distinct(na.omit(depth_seafloor)) == 1) {
      first(na.omit(depth_seafloor))
    } else {
      #print(cur_data_all())  # This prints the full data for the current group
      first(na.omit(depth_seafloor))
    },
    
    lockdown = first(na.omit(lockdown)),
    BiodivMed2023 = first(na.omit(BiodivMed2023)),
    method = first(na.omit(method)),
    country = first(na.omit(country)),
    region = first(na.omit(region)),
    site = first(na.omit(site)),
    subsite = first(na.omit(subsite)),
    component = first(na.omit(component)),
    habitat = first(na.omit(habitat)),
    protection = first(na.omit(protection)),
    mpa_name = first(na.omit(mpa_name)),
    replicates_geometry = first(na.omit(geom)),
    project = first(na.omit(project)),
    filter = first(na.omit(filter)),
    Tele01 = first(na.omit(Tele01)),
    Pleo = first(na.omit(Pleo)),  
    Mamm01 = first(na.omit(Mamm01)),
    Vert01 = first(na.omit(Vert01)),
    X16s_Metazoa = first(na.omit(X16s_Metazoa)),
    Bact02 = first(na.omit(Bact02)), 
    Euka02 = first(na.omit(Euka02)),
    estimated_volume_total = sum(na.omit(estimated_volume)),
    duration_total = sum(na.omit(duration)),
    field_replicates = first(na.omit(field_replicates)),
    PCR_replicates = first(na.omit(PCR_replicates)),
    
    # if comments differ combine them by copy pasting them with their associated replicates id :
    comments = paste(unique(na.omit(comments)), collapse = "; "),  # Combine unique comments with a semicolon
    .groups = "drop"
  )

# Check if mtdt_6 and mtdt_7 have the same columns and print columns that are only in one of the two datasets
cols_mtdt_6 <- colnames(mtdt_6)
cols_mtdt_7 <- colnames(mtdt_7)
only_in_mtdt_6 <- setdiff(cols_mtdt_6, cols_mtdt_7)
only_in_mtdt_7 <- setdiff(cols_mtdt_7, cols_mtdt_6)
print(paste("Columns only in mtdt_6:", paste(only_in_mtdt_6, collapse = ", ")))
print(paste("Columns only in mtdt_7:", paste(only_in_mtdt_7, collapse = ", ")))

# Export
sf::st_write(mtdt_7, "./data/processed_data/Mtdt/mtdt_7.gpkg", delete_dsn = TRUE)


