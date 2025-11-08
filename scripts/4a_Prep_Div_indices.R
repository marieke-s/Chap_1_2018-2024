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






#------------- Load  & format data ---------
# Traits / taxo
traits <- read.csv("./data/raw_data/traits/species_traits_inferred.csv", sep = ";", header=T)

# Format colnames
traits <- traits %>% 
  dplyr::filter(!is.na(Genus)) %>% 
  rename_all(~ str_replace_all(., "\\.", " ")) %>% 
  rename_all(~ str_squish(.)) %>% 
  rename(log_depth_min = `log depth_min`, log_depth_max = `log depth_max`) %>%
  rename("log_geographic_range_Albouy19" = "log geographic_range_Albouy19")

# Presence
presence <- read.csv("./data/processed_data/eDNA/occ_pooled_v1.1.csv", sep = ",", header=T)
names(presence) <- gsub("\\_", " ", names(presence))

# Metadata
meta <- sf::st_read("./data/processed_data/Mtdt/mtdt_6.gpkg")

meta <- st_drop_geometry(meta)




#------------- Handle NA in traits ---------
# Check Missing species in traits -----
# Check species in presence that are not in traits
length(setdiff(colnames(presence[-1]), traits$species)) # 14 species not in traits matrix

setdiff(colnames(presence[-1]), traits$species) # Complex species + 3 species

# 3 species :
# Zu cristatus : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1776 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=126529)

# Encrasicholina punctifer : no synonym found in traits db (https://www.fishbase.se/Nomenclature/558 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=219982) # Species detected for the 1rst time in the Med by Ciftci et al. 2017 (https://doi.org/10.1111/jai.13276)

# Stomias boa : no synonym found in traits db (https://www.fishbase.se/Nomenclature/1806 ; https://www.marinespecies.org/aphia.php?p=taxdetails&id=234601)


# Fill missing species  -----
# Manually complete traits for missing species
traits <- traits %>%
  add_row(species = "Zu cristatus",
          IUCN_category = "LC",
          Length = 118,
          Importance = "of no interest",
          PriceCateg = NA,
          Aquarium = NA,
          Troph = 4, # piscivore
          a = NA,
          b = NA,
          K = NA,
          TempPrefMean = NA,
          Schooling = "no",
          DemersPelag = "bathypelagic",
          body_depth_ratio = NA,
          body_width_ratio = NA,
          IUCN_inferred_Loiseau23 = NA,
          log_depth_min = NA,
          log_depth_max = NA,
          ClimVuln_SSP126 = NA,
          ClimVuln_SSP585 = NA,
          log_geographic_range_Albouy19 = NA,
          spec_code = NA,
          worms_id = 126529,
          fishbase_name = "Zu cristatus",
          Phylum = "Chordata",
          Class = "Teleostei",
          Order = "Lampriformes",
          Family = "Trachipteridae",
          Genus = "Zu") %>%
  add_row(species = "Encrasicholina punctifer",
          IUCN_category = "LC",
          Length = 13,
          Importance = "minor commercial",
          PriceCateg = NA,
          Aquarium = NA,
          Troph = NA,
          a = NA,
          b = NA,
          K = NA,
          TempPrefMean = NA,
          Schooling = "yes",
          DemersPelag = "reef-associated",
          Vulnerability = NA,
          body_depth_ratio = NA,
          body_width_ratio =NA,
          IUCN_inferred_Loiseau23 = NA,
          log_depth_min = NA,
          log_depth_max = NA,
          ClimVuln_SSP126 = NA,
          ClimVuln_SSP585 = NA,
          log_geographic_range_Albouy19 = NA,
          spec_code = NA,
          worms_id = 219982,
          fishbase_name = "Encrasicholina punctifer",
          Phylum = "Chordata",
          Class = "Teleostei",
          Order = "Clupeiformes",
          Family = " Engraulidae",
          Genus = "Encrasicholina") %>%
  add_row(species = "Stomias boa",
          IUCN_category = "LC",
          Length = 32.2,
          Importance = "minor commercial",
          PriceCateg = NA,
          Aquarium = NA,
          Troph = 4, # piscivore
          a = NA,
          b = NA,
          K = NA,
          TempPrefMean = NA,
          Schooling = "yes",
          DemersPelag = " bathypelagic",
          Vulnerability = NA,
          body_depth_ratio = NA,
          body_width_ratio =NA,
          IUCN_inferred_Loiseau23 = NA,
          log_depth_min = NA,
          log_depth_max = NA,
          ClimVuln_SSP126 = NA,
          ClimVuln_SSP585 = NA,
          log_geographic_range_Albouy19 = NA,
          spec_code = NA,
          worms_id = 234601,
          fishbase_name = "Stomias boa" ,
          Phylum = "Chordata",
          Class = "Teleostei",
          Order = " Stomiiformes",
          Family = "  Stomiidae",
          Genus = "Stomias") 
          
# Fill traits for complex species -----


# Check if each separated species is in traits
cp <- setdiff(colnames(presence[-1]), traits$species) # Complex species + 3 species
separated_species <- unique(unlist(strsplit(cp, "\\.")))
setdiff(separated_species, traits$species) # All separated species are in traits

# For each complex species, check if the traits differ 
for (complex_sp in cp) {
  species_list <- unlist(strsplit(complex_sp, "\\."))
  traits_subset <- traits %>% filter(species %in% species_list)
  
  # Print the complex species and its traits
  print(paste("-----------------", complex_sp, "---------------------------------------------"))
  cat("\n")
  print(traits_subset)
  cat("\n")
  cat("\n")
  cat("\n")

}


# Mean numerical traits and merge chr traits
# ---- CONFIG 
# columns that should be forced to NA in the composite row
force_na_cols <- c("fishbase_name", "worms_id", "spec_code", "X")

# merge chr var
resolve_char_col <- function(values, colname) {
  uniq <- unique(values)
  uniq <- uniq[!is.na(uniq)]
  
  if (length(uniq) == 0) {
    return(NA_character_)
  } else if (length(uniq) == 1) {
    return(uniq[1])
  } else if (length(uniq) == 2) {
    message(sprintf(
      "Column '%s' has conflicting values: %s | %s",
      colname, paste0("'", uniq[1], "'"), paste0("'", uniq[2], "'")
    ))
  } else {
    message(sprintf(
      "Column '%s' has >2 conflicting values: %s",
      colname, paste(paste0("'", uniq, "'"), collapse = " | ")
    ))
  }
  
  # manual choice
  cat(sprintf("Choose a value for '%s':\n", colname))
  choices <- c(uniq, "<set NA>")
  for (i in seq_along(choices)) cat(sprintf("  [%d] %s\n", i, choices[i]))
  repeat {
    sel <- suppressWarnings(as.integer(readline("Enter number: ")))
    if (!is.na(sel) && sel >= 1 && sel <= length(choices)) break
    cat("Invalid choice. Try again.\n")
  }
  if (sel == length(choices)) return(NA_character_) else return(uniq[sel])
}

# build a single composite row from a subset of rows
make_composite_row <- function(traits_subset, composite_species_name = NULL, force_na_cols = force_na_cols) {
  if (nrow(traits_subset) == 0) abort("traits_subset is empty.")
  
  # Prepare a named list to become a one-row tibble with same column order
  out <- vector("list", length = ncol(traits_subset))
  names(out) <- names(traits_subset)
  
  for (col in names(traits_subset)) {
    vals <- traits_subset[[col]]
    
    # Force selected columns to NA (type-aware)
    if (col %in% force_na_cols) {
      if (is.character(vals)) {
        out[[col]] <- NA_character_
      } else if (is.integer(vals)) {
        out[[col]] <- NA_integer_
      } else if (is.numeric(vals)) {
        out[[col]] <- NA_real_
      } else {
        out[[col]] <- NA
      }
      next
    }
    
    if (is.numeric(vals)) {
      # mean of numeric columns; NA if all NAs
      if (all(is.na(vals))) {
        out[[col]] <- NA_real_
      } else {
        out[[col]] <- mean(vals, na.rm = TRUE)
      }
    } else if (is.character(vals)) {
      # species column handled below if composite_species_name provided
      if (!is.null(composite_species_name) && identical(col, "species")) {
        out[[col]] <- composite_species_name
      } else {
        out[[col]] <- resolve_char_col(vals, col)
      }
    } else {
      # Fallback: carry first non-NA; else NA
      nn <- vals[!is.na(vals)]
      out[[col]] <- if (length(nn)) nn[1] else NA
    }
  }
  
  # Ensure species column is set, if we have a label
  if (!is.null(composite_species_name) && "species" %in% names(out)) {
    out[["species"]] <- composite_species_name
  }
  
  # Return as a one-row tibble with original column order and types where possible
  as_tibble(out)
}

# ---- append_composites_for_cp
# This will create and append one composite row per complex to `traits`.
append_composites_for_cp <- function(traits, cp) {
  for (complex_sp in cp) {
    species_list <- unlist(strsplit(complex_sp, "\\."))
    traits_subset <- traits %>% filter(species %in% species_list)
    
    cat("\n------------------------------------------------------------\n")
    cat("Complex:", complex_sp, "\n")
    print(traits_subset)
    
    composite_row <- make_composite_row(
      traits_subset = traits_subset,
      composite_species_name = complex_sp,
      force_na_cols = force_na_cols
    )
    
    # Append to traits
    traits <- bind_rows(traits, composite_row)
    cat("-> Added composite row for:", complex_sp, "\n")
  }
  traits
}

# ---- RUNs
traits2 <- append_composites_for_cp(traits, cp)

traits <- traits2


rm(traits2)
rm(traits_subset, cp, complex_species, species_list, separated_species)
rm(make_composite_row, resolve_char_col, append_composites_for_cp, force_na_cols, complex_sp)


# Filter traits to keep only our species  ----
traits <- traits %>%
  dplyr::filter(species %in% colnames(presence[-1]))



# CHECK and RESOLVE : NA in traits -----
sapply(traits, function(x) sum(is.na(x)))





# IUCN NA ----
# Checks
unique(traits$IUCN_category)
unique(traits$IUCN_inferred_Loiseau23)

# NA in both IUCN_category and IUCN_inferred_Loiseau23 
sum(is.na(traits$IUCN_category) & is.na(traits$IUCN_inferred_Loiseau23)) # 6 NA

# Print these species 
traits$species[is.na(traits$IUCN_category) & is.na(traits$IUCN_inferred_Loiseau23)] # "Argentina sphyraena" "Isurus oxyrinchus"   "Oedalechilus labeo" 
# Complexes : "Labrus merula.Labrus viridis"              "Argyrosomus regius.Umbrina cirrosa"        "Dentex dentex.Pagrus auriga.Pagrus pagrus"

# Combine both IUCN : when no category use inferred by Loiseau et al. 2024 (https://doi.org/10.1371/journal.pbio.3002773)
traits <- traits %>%
  mutate(IUCN_combined = ifelse(is.na(IUCN_category), IUCN_inferred_Loiseau23, IUCN_category)) 

# Manually complete these 3 species
traits <- traits %>%
  mutate(IUCN_combined = ifelse(species == "Argentina sphyraena", "NE", IUCN_combined)) %>% # Not evaluated
  mutate(IUCN_combined = ifelse(species == "Isurus oxyrinchus", "EN", IUCN_combined)) %>% # Endangered
  mutate(IUCN_combined = ifelse(species == "Oedalechilus labeo", "NE", IUCN_combined)) # Not evaluated



sum(is.na(traits$IUCN_combined)) # 3 NA 














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




# Importance NA ----
# Check NAs
sum(is.na(traits$Importance)) # 54 NA

# print NA species
NA_imp_sp <- traits$species[is.na(traits$Importance)]

# Check if present in Alicia db 
traits_ad <- read.csv("./data/raw_data/traits/functional_data.csv")

traits_ad_NA_imp <- traits_ad %>%
  filter(Species %in% NA_imp_sp) 




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

unique(traits$Importance_combined)

# Manually complete these 3 species (checked on fish base + FAO) 
traits <- traits %>%
  mutate(Importance_combined = ifelse(species == "Facciolella oxyrhynchus", "of no interest", Importance_combined)) %>%
  mutate(Importance_combined = ifelse(species == "Knipowitschia caucasica", "of no interest", Importance_combined)) %>%
  mutate(Importance_combined = ifelse(species == "Tripterygion melanurus", "of no interest", Importance_combined))

sum(is.na(traits$Importance_combined)) # 0 NA


# Save NA resolved traits -----
# Export traits2
write.table(traits, "./data/processed_data/Traits/species_traits_NA-resolved_v1.0.csv", row.names = F, dec = ".", sep = ";")

# Reload traits -------
traits <- read.csv("./data/processed_data/Traits/species_traits_NA-resolved_v1.0.csv", sep = ";", header=T)


#------------- Biodiversity indicators ---------

#### Setup ####

# Number of cols
n_cols <- 11

# Create indicator matrix
indicators <- matrix(NA, nrow(presence), n_cols,
                     dimnames = list(rownames(presence),
                                     c("replicates",
                                       "R", 
                                       "Crypto",     # Number of Cryptobenthic species
                                       "Elasmo",     # Number of Elasmobranch species
                                       "DeBRa",      # Demersal Benthic Ratio or Crypto
                                       "RedList",    # IUCN
                                       "LRFI",       # Large Fish Species
                                       "TopPred",    # Trophic
                                       "Commercial", # Commercial and Highly commercial species (Importance)
                                       "Grouper",    # Epinephelus marginatus (0/1)
                                       "AngelShark" # Squatina squatina (0/1)
                                     )))             


indicators[,1] <- presence %>% pull(replicates)

rm(n_cols)

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
sum(is.na(traits$DemersPelag)) # 3 NA (complexe species)
sum(is.na(traits$Benthic)) # 0 NA 
sum(is.na(traits$DP)) # 0 NA







#### 5  - Number of species listed VU, EN or CR on the IUCN Red List - RedList ####

traits <- traits %>%
  mutate(RedList = if_else(IUCN_combined %in% c("CR", "EN", "VU", "NT", "Threatened"), 1, 0))

# Check NA
sum(is.na(traits$IUCN_combined)) # 3 NA (complexe species)
sum(is.na(traits$RedList)) # 0 NA 






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
sum(is.na(traits$Troph)) # 1
sum(is.na(traits$TopPred)) # 0

# Print na species
traits$species[is.na(traits$Troph)] # "Encrasicholina punctifer"







#### 8  - Commercial species - Commercial ####

traits <- traits %>%
  mutate(Commercial = if_else(Importance_combined == "subsistence fisheries" | 
                                Importance_combined == "minor commercial" | 
                                Importance_combined == "of no interest" | 
                                Importance_combined == "of potential interest", 0,
                              if_else(Importance_combined == "commercial" |
                                        Importance_combined == "highly commercial" , 1, NA_real_
                              )))

# check NAs
sum(is.na(traits$Importance_combined)) # 0 NA
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
  indicators[i,"Commercial"] <- sum(traits[which(traits$species %in% s_i), "Commercial"], na.rm = T) # Commercial
  
}


#### 12 - Grouper ####

indicators[, "Grouper"] <- ifelse(presence$`Epinephelus marginatus` != 0, 1, 0)
summary(factor(indicators[,"Grouper"]))



#### 13 - AngelShark ####

indicators[, "AngelShark"] <- ifelse(presence$`Squatina squatina` != 0, 1, 0)
summary(factor(indicators[,"AngelShark"]))

#### Export results ####


indicators <- indicators %>% 
  data.frame(.) %>% 
  tibble::rownames_to_column(var = "ID") %>% 
  dplyr::select(-ID) 


# v1.0 ----
# Based on occ_pooled_v1.1 and species_traits_NA-resolved_v1.0.csv (MO traits 2023 extracted by Ulysse)
# 9 Div indices : "R"          "Crypto"     "Elasmo"     "DeBRa"      "RedList"    "LRFI"       "TopPred"    "Commercial" "Grouper"   "AngelShark"
write.table(indicators, "./data/processed_data/Traits/div_indices_v1.0.csv", row.names = F, dec = ".", sep = ";")




#---------- CHECKS OF DIV INDICES --------------
# Load data ----
indicators <- read.csv("./data/processed_data/Traits/div_indices_v1.0.csv", sep = ";", header=T)

# Check NAs for each column ----
sapply(indicators, function(x) sum(is.na(x))) 

# Plots ----

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

# Map indices ----

# Add coordinates to indicators
coords <- sf::st_read("./data/processed_data/predictors/predictors_raw_v1.1.gpkg") %>%
  dplyr::select(c(x,y,replicates)) %>%
  st_drop_geometry() %>%
  left_join(indicators, by = "replicates") %>%
  filter(replicates %in% indicators$replicates)

indicators <- coords
rm(coords)

# Load world map data using rnaturalearth
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# List of indices columns 
index_columns <- indicators %>%
  dplyr::select(-c(replicates, x, y)) %>%
  colnames()


# Specify the output directory to save the maps
output_directory <- "./figures/Div_indices/Maps/"

# Loop over each index column to create the map + frequency distribution 
for (col_name in index_columns) {
  
  # Calculate stats for the current indicator column
  mean_val <- mean(indicators[[col_name]], na.rm = TRUE)
  sd_val   <- sd(indicators[[col_name]], na.rm = TRUE)
  n_samp   <- sum(!is.na(indicators[[col_name]]))
  
  # Map for each indicator column
  map_plot <- ggplot() +
    geom_sf(data = world, fill = "white", color = "lightblue") +  # Background world map
    geom_point(
      data = indicators,
      aes(x = x, y = y, color = !!sym(col_name)),
      size = 1, alpha = 0.8
    ) +
    scale_color_gradient(low = "yellow", high = "red") +          # Color gradient
    ggtitle(paste("Map of", col_name)) +                          # Title
    xlab("Longitude") + ylab("Latitude") +                        # Axis labels
    coord_sf(
      xlim = c(min(indicators$x, na.rm = TRUE) - 0.1,
               max(indicators$x, na.rm = TRUE) + 0.0),
      ylim = c(min(indicators$y, na.rm = TRUE) - 0.0,
               max(indicators$y, na.rm = TRUE) + 0.0)
    ) +                                                           # Map limits
    theme_test() +                                                # Theme
    theme(panel.background = element_rect(fill = "#e0edff"))      # Panel background
  
  # Histogram (using Frequency instead of Density)
  hist_plot <- ggplot(indicators, aes(x = !!sym(col_name))) +
    geom_histogram(
      aes(y = after_stat(count), fill = after_stat(count)),       # <-- frequency
      bins = 50, color = "black", linewidth = 0.01
    ) +
    scale_fill_gradient(low = "#80ffdb", high = "#03045e") +
    ggtitle(paste("Distribution of", col_name)) +
    xlab(col_name) + ylab("Frequency") +                          # <-- label updated
    theme_minimal() +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste(
        "Mean:", round(mean_val, 2),
        "\nSD:", round(sd_val, 2),
        "\nN:", n_samp
      ),
      hjust = 1.1, vjust = 1.5, size = 5, color = "black"
    )
  
  # Combine map and histogram (map on top, histogram below)
  combined_plot <- grid.arrange(
    map_plot, hist_plot,
    layout_matrix = rbind(c(1), c(2)),
    heights = c(2, 1)  # Map 2/3 height, histogram 1/3
  )
  
  # Print combined plot
  print(combined_plot)
  
  # Save combined plot (taller aspect)
  ggsave(
    filename = paste0(output_directory, "map_", col_name, ".png"),
    plot = combined_plot, width = 8, height = 10, dpi = 300
  )
}






