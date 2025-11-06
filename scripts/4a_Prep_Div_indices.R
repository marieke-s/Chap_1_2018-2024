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

# Libraries
library(ape)
library(dplyr)
library(fastDummies)
library(fishtree)
library(geiger)
library(ggplot2)
library(gridExtra)
library(igraph)
library(picante)
library(readr)
library(reshape2)
library(rlang)
library(sf)
library(svglite)
library(tidyr)
library(tibble)



# Load functions
source("./utils/Fct_Data-Prep.R")





#------------- COMPUTE DIVERSITY INDICES [ SCRIPT ALICIA-MK ] ----------------------

# See Dalongeville et al. 2022 (https://doi.org/10.1111/1365-2664.14276) 
# for descriptions of the indices.

## ------------------- Data Preparation -------------------
# Load occurrence data
occ <- read.csv("./data/processed_data/eDNA/occ_pooled.csv")

# Remove replicates with no species detected
occ <- occ[rowSums(occ[,-1]) > 0, ]

# Transpose occurrence data
adne <- as.data.frame(t(occ[, -1])) 
colnames(adne) <- occ[, 1]

# Replace dots with spaces in row names
rownames(adne) <- gsub("\\_", " ", rownames(adne))

# List of species and samples
species <- rownames(adne)
sample <- colnames(adne)

# Load species traits matrix
traits <- read.csv("./data/raw_data/traits/functional_data.csv")

# Create the results dataframe for diversity indices
indicators <- as.data.frame(matrix(NA, ncol = 12, nrow = ncol(adne),
                                   dimnames = list(colnames(adne),
                                                   c("spygen_code", "R", "FD", "LFI", "Crypto",
                                                     "DP_B_ratio", "RedList", "Chondri",
                                                     "Commercial", "High_commerc", "PD", "Vulner"))))

# Add sample IDs
indicators$spygen_code <- rownames(indicators)

## ------------------- 1. Species Richness (R) -------------------

indicators$R <- apply(adne, 2, sum)

## ------------------- 2. Functional Diversity (FD) -------------------

for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- adne %>%
    filter(adne[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Get the functional groups of these species
  fd_i <- as.factor(traits[which(traits$Species %in% s_i), "GF"])
  
  # Number of unique functional groups
  indicators[i, "FD"] <- nlevels(fd_i)
}

## ------------------- 3. Large Reef Fish Index (LFI) -------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "LFI"] <- sum(traits[which(traits$Species %in% s_i), "LRFI"]) 
}

## ------------------- 4. Ratio Demerso-pelagic / Benthic (DP_B_ratio) -------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  
  indicators[i, "DP_B_ratio"] <- sum(traits[which(traits$Species %in% s_i), "DP"]) / (sum(traits[which(traits$Species %in% s_i), "B"])+1)
}

## ------------------- 5. Chondrichthyan Species (Chondri) -------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "Chondri"] <- sum(traits[which(traits$Species %in% s_i), "SHarK"])
}

## ------------------- 6. Commercial Species (Commercial) -------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "Commercial"] <- sum(traits[which(traits$Species %in% s_i), "all_commercial_level"])
}

gc()

## ------------------- 7. Highly Commercial Species (High_commerc) -------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "High_commerc"] <- sum(traits[which(traits$Species %in% s_i), "highly_commercial_only"])
}

## ------------------- 8. Cryptobenthic Species (Crypto) -------------------

crypto_families <- c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                     "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                     "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

traits <- traits %>%
  mutate(crypto_Brandl = if_else(Family %in% crypto_families, 1, 0))

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "Crypto"] <- sum(traits[which(traits$Species %in% s_i), "crypto_Brandl"]) 
}

## ------------------- 9. Vulnerability Index (Vulner)-------------------

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  indicators[i, "Vulner"] <- mean(traits[which(traits$Species %in% s_i), "Vulnerability"], na.rm = TRUE)
}

gc()





# Checking rows with NA values 
# Identify rows with missing Vulnerability values
missing_rows <- which(is.na(indicators$Vulner))
print(missing_rows)

# Check which species were detected in these missing rows
for (i in missing_rows) {
  print(paste("Sample:", rownames(indicators)[i]))
  print("Species detected:")
  print(adne[, i][adne[, i] == 1])  # Shows detected species for that sample
}

# Check if those species exist in the traits dataset
is.na(traits$Vulnerability)

for (i in missing_rows) {
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  print(paste("Sample:", rownames(indicators)[i]))
  print("Species found in traits dataset:")
  print(traits[which(traits$Species %in% s_i), "Species"])
}

traits %>%
  filter(Species %in% c("Pegusa nasuta", "Pteroplatytrygon violacea", "Sardina pilchardus", 
                        "Sardinella aurita", "Sarpa salpa", "Sparus aurata", "Scomber scombrus")) %>%
  dplyr::select(Species, Vulnerability)







# Complete NA values with the average vulnerability of the family
# Compute average vulnerability per family
family_vuln <- traits %>%
  group_by(Family) %>%
  summarise(mean_vuln = mean(Vulnerability, na.rm = TRUE))

# Merge with original traits to fill missing values
traits <- traits %>%
  left_join(family_vuln, by = "Family") %>%
  mutate(Vulnerability = ifelse(is.na(Vulnerability), mean_vuln, Vulnerability)) %>%
  dplyr::select(-mean_vuln)

# Now re-run the original loop
for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  vuln_values <- traits[which(traits$Species %in% s_i), "Vulnerability"]
  
  if (length(vuln_values) == 0 || all(is.na(vuln_values))) {
    indicators[i, "Vulner"] <- median_vuln  
  } else {
    indicators[i, "Vulner"] <- mean(vuln_values, na.rm = TRUE)
  }
}

gc()







## ------------------- 10. Red List IUCN Status (RedList) -------------------
# Make a dummy variable for IUCN categories
traits <- fastDummies::dummy_cols(traits, select_columns = 'IUCN_Red_List_Category')

for (i in 1:nrow(indicators)) { 
  s_i <- adne %>%
    filter(adne[, i] == 1) %>%
    tibble::rownames_to_column(var = "Sp") %>%
    pull(Sp)
  
  # Number of species per category
  VU <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_VU"], na.rm = TRUE)
  EN <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_EN"], na.rm = TRUE)
  CR <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_CR"], na.rm = TRUE)
  
  # Calculate indicator
  indicators[i, "RedList"] <- VU + EN + CR
}

gc()

## ------------------- 11. Phylogenetic Diversity (PD) -------------------

# Retrieve the phylogeny of species across all three oceans
phy <- fishtree::fishtree_phylogeny(species = species)

# Plot Phylogeny of species
plot(phy, show.tip.label = TRUE, cex = 0.5)

# Format species names to match phylogeny
rownames(adne) <- gsub(" ", "_", rownames(adne), fixed = TRUE)

# Check that phylogeny and data have matching names
nc <- geiger::name.check(phy, adne) 
print(sort(nc$data_not_tree)) 
print(length(nc$data_not_tree)) # 26 species not in tree (chondrychtyans + synonyms)
print(sort(nc$data_not_tree))

# Manually check synonyms and find the species
species[species == "Chelon auratus"] <- "Liza aurata"
rownames(adne)[rownames(adne) == "Chelon_auratus"] <- "Liza_aurata"

species[species == "Chelon ramada"] <- "Liza ramada"
rownames(adne)[rownames(adne) == "Chelon_ramada"] <- "Liza_ramada"

species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
rownames(adne)[rownames(adne) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"

species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
rownames(adne)[rownames(adne) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"

# species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
# rownames(adne)[rownames(adne) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"

# Retrieve the missing phylogeny 
phy <- fishtree_phylogeny(species = species, type = "phylogram")
nc <- geiger::name.check(phy, adne) 
print(sort(nc$data_not_tree)) 
print(length(nc$data_not_tree)) # 22 species not in tree 

# Remove from the data the species that are not in the tree
adne2 <- adne[which(rownames(adne) %in% nc$data_not_tree == F),]

# Transpose the ADNe matrix 
adne2 <- t(adne2)

# prune the tree
prunedTree <- picante::prune.sample(adne2,phy)

# Calculate Faith's PD
pd.result <- pd(adne2, prunedTree, include.root=T)

# Add PD to indicator dataframe
indicators[,"PD"] <- pd.result$PD/pd.result$SR  # percentage of species richness with only species that are in the tree



rm(adne, adne2, phy, prunedTree, pd.result, nc, traits)



##-------------------- CHECKS : NA values -------------------
# Check if there NA values in indicators
na_values <- apply(indicators, 2, function(x) sum(is.na(x)))
print(na_values)






#------------- COMPUTE DIVERSITY INDICES [ SCRIPT MARIE ] ----------------------
#### Load  & format data ####

# Traits / taxo
traits <- read.csv("./data/raw_data/traits/species_traits_inferred.csv", sep = ";", header=T)

# Format colnames
traits <- traits %>% 
  dplyr::filter(!is.na(Genus)) %>% 
  rename_all(~ str_replace_all(., "\\.", " ")) %>% 
  rename_all(~ str_squish(.)) %>% 
  rename(log_depth_min = `log depth_min`, log_depth_max = `log depth_max`)

# Presence
presence <- read.csv("./data/processed_data/eDNA/occ_pooled_v1.0.csv", sep = ",", header=T)
names(presence) <- gsub("\\_", " ", names(presence))

# Metadata
meta <- sf::st_read("./data/processed_data/Mtdt/mtdt_6.gpkg")

meta <- st_drop_geometry(meta)




# CHECK : Missing species in traits -----
# Check species in presence that are not in traits
length(setdiff(colnames(presence[-1]), traits$species)) # 11 species not in traits matrix

setdiff(colnames(presence[-1]), traits$species) # Species from species complex (eg: "Parablennius sp 1" , "Parablennius sp 2") + 3 species

# 3 species :
# Zu cristatus : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1776 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=126529)

# Encrasicholina punctifer : no synonym found in traits db (https://www.fishbase.se/Nomenclature/558 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=219982)

# Stomias boa : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1806 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=234601)


# Filter traits to keep only our species  ----
traits <- traits %>%
  dplyr::filter(species %in% colnames(presence[-1]))



# CHECK and RESOLVE : NA in traits -----
# count Na per column
sapply(traits, function(x) sum(is.na(x)))


# NA in ICUN indicators and IUCN_Inferred_Loiseau23
# NA in Importance
# NA in Length
# NA worms_id
# NA body_width_ratio


# IUCN NA ----
# Checks
unique(traits$IUCN_category)
unique(traits$IUCN_inferred_Loiseau23)

# NA in both IUCN_category and IUCN_inferred_Loiseau23
sum(is.na(traits$IUCN_category) & is.na(traits$IUCN_inferred_Loiseau23)) # 3 NA

# Print these species 
traits$species[is.na(traits$IUCN_category) & is.na(traits$IUCN_inferred_Loiseau23)] # "Argentina sphyraena" "Isurus oxyrinchus"   "Oedalechilus labeo" 

# Combine both
traits <- traits %>%
  mutate(IUCN_combined = ifelse(is.na(IUCN_inferred_Loiseau23), IUCN_category, IUCN_inferred_Loiseau23)) 

# Manually complete these 3 species
traits <- traits %>%
  mutate(IUCN_combined = ifelse(species == "Argentina sphyraena", "NE", IUCN_combined)) %>% # Not evaluated
  mutate(IUCN_combined = ifelse(species == "Isurus oxyrinchus", "EN", IUCN_combined)) %>% # Endangered
  mutate(IUCN_combined = ifelse(species == "Oedalechilus labeo", "NE", IUCN_combined)) # Not evaluated



sum(is.na(traits$IUCN_combined)) # 0 NA 














# Length NA ----
sum(is.na(traits$Length)) # 1

# print NA species
traits$species[is.na(traits$Length)] # "Aetomylaeus bovinus"

# Fishbase : 
# Length at first maturity / Size / Weight / Age
# Maturity: Lm 90.0, range 35 - 148 cm
# Max length : 222.0 cm WD (female); common length : 150 cm WD male/unsexed; (Ref. 57025); max. published weight: 116.0 kg (Ref. 85836) 

# Complete with max length
traits <- traits %>%
  mutate(Length = ifelse(species == "Aetomylaeus bovinus", 222.0, Length))




# Importance NA [TO RESOLVE]----
# check NAs
sum(is.na(traits$Importance)) # 54 NA

# print NA species
NA_imp_sp <- traits$species[is.na(traits$Importance)]

# Check if present in Alycia db 
traits_ad <- read.csv("./data/raw_data/traits/functional_data.csv")


traits_ad_NA_imp <- traits_ad %>%
  filter(Species %in% NA_imp_sp) 


unique(traits$Importance)


# Resolve 
# 1. join traits_ad$all_commercial_level and $highly_commercial_only to traits
# 2. mutate combined_importance = ifelse(is.na(Importance),ifelse(highly_commercial_only == 1, "highly commercial", ifelse(all_commercial_level == 1, "commercial", "of no interest")), Importance)

traits <- traits %>%
  left_join(traits_ad %>% dplyr::select(Species, all_commercial_level, highly_commercial_only), 
            by = c("species" = "Species")) %>%
  mutate(Importance_combined = ifelse(is.na(Importance),
                                      ifelse(highly_commercial_only == 1, "highly commercial",
                                             ifelse(all_commercial_level == 1, "commercial", "of no interest")),
                                      Importance))

sum(is.na(traits$Importance_combined)) # 3 NA

# print NA species
traits$species[is.na(traits$Importance_combined)] # "Facciolella oxyrhynchus" "Knipowitschia caucasica" "Tripterygion melanurus"

# Manually complete these 3 species //TO DO 
traits <- traits %>%
  mutate(Importance_combined = ifelse(species == "Facciolella oxyrhynchus", "of no interest", Importance_combined)) %>%
  mutate(Importance_combined = ifelse(species == "Knipowitschia caucasica", "of no interest", Importance_combined)) %>%
  mutate(Importance_combined = ifelse(species == "Tripterygion melanurus", "of no interest", Importance_combined))




# Schooling NA [TO RESOLVE]----
# body_depth_ratio NA [TO RESOLVE]----
# body_width_ratio NA [TO RESOLVE]----
# log_depth_min NA [TO RESOLVE]----
# log_depth_max NA [TO RESOLVE]----
#### Biodiversity indicators ####

# Create indicator matrix
indicators <- matrix(NA,nrow(presence), 16,
                     dimnames = list(rownames(presence),
                                     c("code_spygen",
                                       "R", 
                                       "Crypto",     # Number of Cryptobenthic species
                                       "Elasmo",     # Number of Elasmobranch species
                                       "DeBRa",      # Demersal Benthic Ratio or Crypto
                                       "RedList",    # IUCN
                                       "LRFI",       # Large Fish Species
                                       "TopPred",    # Trophic
                                       "Vulner",     # Vulnerable
                                       "Clim_vuln",  # Vulnerable to climate change
                                       "Commercial", # Commercial and Highly commercial species (Importance)
                                       "Phill",      # Phylogenetic metrics
                                       "Fhill",      # Functional metric
                                       "Thill",      # Taxonomic metric
                                       "Grouper",    # Epinephelus marginatus (0/1)
                                       "Angel_shark" # Squatina squatina (0/1)
                                     )))             


indicators[,1] <- presence %>% pull(replicates)


#### 1  - Species richness - R ####

indicators[,2] <- rowSums(presence[,-1])


#### 2  - Cryptobenthic species - Crypto ####
# Definition Brandl et al 2018
# Brandl SJ, Goatley CHR, Bellwood DR, Tornabene L. 2018 The hidden half: 
# ecology and evolution of cryptobenthic fishes on coral reefs. Biol. Rev. 93 
# doi:10.1111/brv.12423

crypto_families = c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                    "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                    "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

traits <- traits %>%
  mutate(crypto = if_else(Family %in% crypto_families, 1,0))
## crée une nouvelle colonne dans TRAITS et y ajoute la valeur 0 ou 1 selon présence/absence de la famille dans traits


# Check NAs
sum(is.na(traits$crypto)) # 0 NA 







#### 3  - Elasmobranch species - Elasmo ####

traits <- traits %>%
  mutate(elasmo = if_else(Class == "Elasmobranchii", 1, 0))


# Check NAs
sum(is.na(traits$elasmo)) # 0 NA 







#### 4  - Demerso-pelagic species ratio - DeBRa ####

traits <- traits %>%
  mutate(Benthic = ifelse(DemersPelag %in% c("reef-associated", "demersal", "bathydemersal"), 1, 0)) %>% 
  mutate(DP = ifelse(DemersPelag %in% c("pelagic-oceanic", "pelagic-neritic", "bathypelagic", "pelagic"), 1, 0)) # Create a new Dermerso-Pelagic column which will take the value 1 if in the DemersPelag column it is represented by pelagic-oceanic, pelagic-neritic, bathypelagic. Otherwise it will take the value 0.


# Check NAs
sum(is.na(traits$DemersPelag)) # 0 NA
sum(is.na(traits$Benthic)) # 0 NA 
sum(is.na(traits$DP)) # 0 NA







#### 5  - Number of species listed VU, EN or CR on the IUCN Red List - RedList ####

traits <- traits %>%
  mutate(RedList = if_else(IUCN_combined %in% c("CR", "EN", "VU", "NT", "Threatened"), 1, 0))

# Check NA
sum(is.na(traits$IUCN_combined)) # 0 NA 





# Check NAs
sum(is.na(traits$RedList)) # 0 NA 

# print NA species for IUCN_category
traits$species[is.na(traits$IUCN_category)]

# Check if these species as inside presence
length(setdiff(colnames(presence[-1]), traits$species[is.na(traits$IUCN_category)])) # /!\ 224 species with no IUCN status are inside our presence matrice






#### 6  - Large Reef Fish Indicator (max length > 20cm) - LRFI ####

traits <- traits %>%
  mutate(LRFI = if_else(Length >= 20, 1,0)) 

# check NAs
sum(is.na(traits$Length)) # 0
sum(is.na(traits$LRFI)) # 0











#### 7  - Top Predators (piscivore & max length > 50cm) - TopPred ####

traits <- traits %>%
  mutate(TopPred = if_else(Length >= 50 & Troph >= 4, 1,0))

# check NAs
sum(is.na(traits$Troph)) # 0
sum(is.na(traits$TopPred)) # 0









#### 8  - Commercial species - Commercial ####

# A modifier avec la richesse des espèces commerciales

traits <- traits %>%
  mutate(Commercial = if_else(Importance == "subsistence fisheries" | 
                                Importance == "minor commercial" | 
                                Importance == "of no interest" | 
                                Importance == "of potential interest", 0,
                              if_else(Importance == "commercial" |
                                        Importance == "highly commercial" , 1, NA_real_
                              )))

# check NAs
sum(is.na(traits$Importance)) # 54 NA
sum(is.na(traits$Commercial)) # 0 NA









#### Compute indicators ####

for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- presence %>%
    dplyr::select_if(presence[i,] == 1)
  s_i <- names(s_i)
  s_i <- gsub("_", " ", s_i)
  
  length(s_i)
  
  # Calculate indicator
  indicators[i,"Crypto"] <- sum(traits[which(traits$species %in% s_i), "crypto"], na.rm = T)
  indicators[i,"Elasmo"] <- sum(traits[which(traits$species %in% s_i), "elasmo"], na.rm = T)
  indicators[i,"DeBRa"] <- sum(traits[which(traits$species %in% s_i), "DP"], na.rm = T) / (sum(traits[which(traits$species %in% s_i), "Benthic"], na.rm = T) + 1) 
  indicators[i,"RedList"] <- sum(traits[which(traits$species %in% s_i), "RedList"], na.rm = T)
  indicators[i,"LRFI"] <- sum(traits[which(traits$species %in% s_i), "LRFI"], na.rm = T) 
  indicators[i,"TopPred"] <- sum(traits[which(traits$species %in% s_i), "TopPred"], na.rm = T) 
  indicators[i,"Vulner"] <- mean(traits[which(traits$species %in% s_i), "Vulnerability"], na.rm = T) 
  indicators[i,"Clim_vuln"] <- mean(traits[which(traits$species %in% s_i), "ClimVuln_SSP585"], na.rm = T)
  indicators[i,"Commercial"] <- sum(traits[which(traits$species %in% s_i), "Commercial"], na.rm = T) # Commercial
  
}


# Important methodological considerations : NA values are removed for indicators compution. In some case up to 





#### Build data for next indicators  ####

# list species present in the sample
s_i <- colnames(presence[-1])

# select traits
traits_fd <- traits %>%
  dplyr::select(species,
                a, 
                b, 
                K, 
                Schooling,
                Troph,
                body_depth_ratio,
                body_width_ratio,
                log_depth_min,
                log_depth_max,
                DemersPelag,
  ) %>%
  dplyr::filter(species %in% s_i) %>% # remove species that are not in the eDNA data
  distinct(species, .keep_all = TRUE) %>% 
  tibble::column_to_rownames(var = "species") %>%
  mutate_if(sapply(., is.character), as.factor) %>%
  mutate_if(sapply(., is.integer), as.factor)

str(traits_fd)

# check species in presence not in traits
length(setdiff(s_i, rownames(traits_fd))) # 21 species not in traits matrix
setdiff(s_i, rownames(traits_fd)) # Species complex + 3 species

# 3 species :
# Zu cristatus : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1776 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=126529)

# Encrasicholina punctifer : no synonym found in traits db (https://www.fishbase.se/Nomenclature/558 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=219982)

# Stomias boa : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1806 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=234601)


# Build the traits category matrix
traits_cat <- cbind.data.frame(trait_name = colnames(traits_fd),
                               trait_type = c("Q", "Q", "Q", "N", "Q","Q","Q", "Q", "Q", "N"))

# Species traits summary:
traits_summ <- mFD::sp.tr.summary(traits_cat, traits_fd, stop_if_NA = F)
traits_summ$"tr_types"
traits_summ$"mod_list"

####################################################@
# functional distance
sp_dist <- mFD::funct.dist(
  sp_tr         = traits_fd,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = F)



save(sp_dist, file = "./data/processed_data/Traits/functional_distances.Rdata")

### Remove the species that are not in traits from the presence matrix
# and reorder

presence_fd <- presence %>% 
  # select(all_of(s_i)) %>%
  dplyr::select(all_of(rownames(traits_fd))) %>% 
  as.matrix(.) 

presence_fd[,order(colnames(presence_fd))]




#### 9  - functional - Fhill  ####
# q=0 ; thau = mean 

# Compute alpha FD Hill with q = 0:
FD0mean <- mFD::alpha.fd.hill(
  asb_sp_w = presence_fd, 
  sp_dist  = sp_dist, 
  tau      = "mean", 
  q        = 0)

indicators[, "Fhill"] <- round(FD0mean$"asb_FD_Hill", 2)


#### 10 - Phylogenetic Diversity - Phill ####

# Retrieve the phylogeny of only native reef species across all three oceans.
phy <- fishtree::fishtree_phylogeny(species = s_i)                                         # Requested 241 but only found 167 species, 74 not in tree (Phylo does not contain Elasmobranch)

## check that phylogeny and data have matching names (same species tree / dataset)
# Format the data
presence_phy <- t(presence) 
rownames(presence_phy) <- gsub(" ", "_", rownames(presence_phy))

# Remove from the data the species that are not in the tree
nc <- geiger::name.check(phy, presence_phy)                                     # 74 species not in tree
missing <- gsub("_", " ", nc$data_not_tree, fixed = TRUE)

presence_pd <- presence %>% 
  dplyr::select(- all_of(missing)) 

# presence is converted to probability vector, summing to 1: 
ps <- presence_pd %>%
  mutate(across(everything(.), ~ entropart::as.ProbaVector(.x))) %>%
  as.matrix(.)
colnames(ps) <- gsub(" ", "_", colnames(ps), fixed = TRUE)

## Compute PD for each survey 
Phill <- pbapply::pbapply(ps, 1, ChaoPD,
                          q = 1, PhyloTree = phy, Normalize = TRUE, Prune = FALSE, CheckArguments = F)

indicators[, "Phill"] <- round(Phill, 2)

#### 11 - Taxonomic - Thill ####
# q=0 ; thau = min
# Compute alpha FD Hill with q = 0:

FD0min <- mFD::alpha.fd.hill(
  asb_sp_w = presence_fd, 
  sp_dist  = sp_dist, 
  tau      = "min", 
  q        = 0)

indicators[, "Thill"] <- round(FD0min$"asb_FD_Hill", 2)



#### 12 - Grouper ####

indicators[, "Grouper"] <- ifelse(presence$`Epinephelus marginatus` != 0, 1, 0)
summary(factor(indicators[,"Grouper"]))

## au cas ou pas de grouper dans la zone mettre des 0 
# indicators[, "Grouper"] <- 0


#### 13 - Angel shark ####

indicators[, "Angel_shark"] <- ifelse(presence$`Squatina squatina` != 0, 1, 0)
summary(factor(indicators[,"Angel_shark"]))

indicators[, "Angel_shark"] <- 0


#### Export results ####


indicators <- indicators %>% 
  data.frame(.) %>% 
  tibble::rownames_to_column(var = "ID") %>% 
  dplyr::select(-ID)

# Merge indicators and metadata 
indicators <- indicators %>%
  rename(replicates = code_spygen)
ind_meta <- merge(indicators, meta, by = "replicates")
ind_meta <- ind_meta %>% distinct(code_spygen, .keep_all = TRUE)

write.table(indicators, "./data/processed_data/Traits/indicators-with-NA.csv", row.names = F, dec = ".", sep = ";")


#---------- CHECKS OF DIV INDICES --------------
# Check NAs for each column ----
sapply(indicators, function(x) sum(is.na(x))) # No NA exept for Fhill ad Thill 

# Plots ----

# Hist + summary

## Make a numeric copy of indicators (except 'replicates')
indicators_num <- indicators %>%
  mutate(across(-replicates, ~ as.numeric(.)))

## List of columns to plot
cols_to_plot <- names(indicators_num)[names(indicators_num) != "replicates"]

## Filter out columns that are all NA
cols_to_plot <- cols_to_plot[sapply(indicators_num[cols_to_plot], function(x) any(!is.na(x)))]

## Compute histogram plots safely
plots_list <- lapply(cols_to_plot, function(col_name) {
  gg_hist_summary(indicators_num[[col_name]], col_name = col_name)
})

## Name the list for easier reference
names(plots_list) <- cols_to_plot

# Display plots
plots_list










