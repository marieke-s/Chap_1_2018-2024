# ──────────────────────────────────────────────────────
# Script: Clean and format eDNA detection dataset (SPYGEN)
# Author: Célia Bertrand
# Date: 10/27/2025
# Purpose: Clean and prepare raw data from SPYGEN for analysis
# ──────────────────────────────────────────────────────

# ─── Load Required Libraries ───────────────────────────
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# Remove all env exept "raw"
rm(list = setdiff(ls(), "raw"))


# ─────────────────────────────────────────────
# ─── Step 1: Load Raw SPYGEN Data ──────────────────────
# ─────────────────────────────────────────────
# Load raw xlsx file WITHOUT headers
raw <- readxl::read_xlsx(
  path = "./data/raw_data/eDNA/Results_Teleo_2018-2024_with_IPOCOM.xlsx",
  col_names = FALSE
)

# Find the row containing column headers (i.e., "scientific_name")
header_row <- which(
  apply(raw, 1, function(x) any(grepl("^scientific_name$", x, ignore.case = TRUE)))
)
if(length(header_row) == 0) stop("Header row with 'scientific_name' not found.")

# Extract species rows (below header)
species_rows <- (header_row + 1):nrow(raw)
species_names_raw <- as.character(raw[[4]][species_rows])


# ─────────────────────────────────────────────
# ─── Step 2: Identify and Remove Problematic Taxa ───────
# ─────────────────────────────────────────────

# Define cleaning mask:
# Remove:
#  - Entries with only ONE word (e.g. "Mugilidae" or "Raja")
#  - NA or empty entries
# Keep:
#  - Normal binomials (e.g. "Scomber scombrus")
#  - Complex detections with underscores (e.g. "C. heterurus_H. speculiger")

bad_mask <- (
  is.na(species_names_raw) |                         # Missing
    species_names_raw == "" |                          # Empty
    # grepl("dae", species_names_raw, ignore.case = TRUE) |  # Families
    sapply(strsplit(species_names_raw, "[ _]"), length) < 2  # Only one token (no space or underscore)
)

# Optionally print examples of removed taxa
if(any(bad_mask)) {
  cat("Examples of problematic taxon rows:\n")
  print(unique(species_names_raw[bad_mask]))
}

# Remove problematic taxon rows from raw data
bad_species_rows <- species_rows[bad_mask]
raw_clean <- raw[-bad_species_rows, ]


# ─────────────────────────────────────────────
# ─── Step 3: Format Dataset for Analysis ────────────────
# ─────────────────────────────────────────────

# Suppress 4 first lines 
raw_clean <- raw_clean[-c(1:4), ]

# Transpose dataset so rows become columns
data_all <- data.frame(t(raw_clean), stringsAsFactors = FALSE)

# Fill column 1, lines 2–5, with "sample"
if (nrow(data_all) >= 5) {
  data_all[2:5, 1] <- "sample"
}
# Suppress lines 1 to 3 and 5 (keeping others)
rows_to_remove <- c(1:3, 5)
rows_to_remove <- rows_to_remove[rows_to_remove <= nrow(data_all)]
data_all <- data_all[-rows_to_remove, ]

# Assign species names as column headers
## ── Promote first row to header & rename first two columns
hdr <- as.character(unlist(data_all[1, ]))
# Force the first two names to what they really are
hdr[1:2] <- c("sample", "metric")

# Replace any empty headers with placeholders
hdr[hdr == "" | is.na(hdr)] <- paste0("unnamed_", seq_len(sum(hdr == "" | is.na(hdr))))

# Make safe/unique column names 
colnames(data_all) <- hdr

# Drop the header row from the data
data_all <- data_all[-1, , drop = FALSE]

## ── Coerce types
# Keep sample/metric as character
data_all$sample <- as.character(data_all$sample)
data_all$metric <- as.character(data_all$metric)

# Everything else should be numeric detections
species_cols <- setdiff(names(data_all), c("sample", "metric"))
data_all[species_cols] <- lapply(
  data_all[species_cols],
  function(x) (as.numeric(stringr::str_replace_all(x, "\\s+", "")))
)

# Replace NA with 0
data_all[is.na(data_all)] <- 0


# ─────────────────────────────────────────────
# ─── Step 4: Clean dataset ───────────────────────────────
# ─────────────────────────────────────────────

# (1) 
# ─────────────────────────────────────────────
# ─────── Merge duplicated species ────────
# ─────────────────────────────────────────────
# Get species columns (skip metadata)
species_cols <- setdiff(names(data_all), c("sample", "metric"))
species_base <- gsub("\\.\\d+$", "", species_cols)

# Merge duplicates (issue is that duplicated species column are not numerical)
merged_count <- 0

for (sp in unique(species_base)) {
  dupes <- grep(paste0("^", sp, "(\\.|$)"), names(data_all))
  
  if (length(dupes) > 1) {
    # ---- Extract duplicate columns ----
    sub <- data_all[, dupes, drop = FALSE]
    
    # ---- Force conversion to numeric ----
    sub[] <- lapply(sub, function(x) {
      x <- as.character(x)
      x <- stringr::str_replace_all(x, "[[:space:]]|,", "")
      suppressWarnings(as.numeric(x))
    })
    
    # ---- Merge duplicates safely ----
    merged_vals <- rowSums(sub, na.rm = TRUE)
    
    # ---- Remove duplicates by position (base R, no dplyr) ----
    data_all <- data_all[, -dupes, drop = FALSE]
    
    # ---- Add back the merged numeric column ----
    data_all[[sp]] <- merged_vals
    
    merged_count <- merged_count + 1
    message("Merged duplicated columns for species: ", sp)
  }
}

# All species columns should now be numeric
non_numeric_cols <- names(data_all)[!sapply(data_all, is.numeric)]
non_numeric_cols <- setdiff(non_numeric_cols, c("sample", "metric"))
print(non_numeric_cols)

# Check if any duplicates remain
anyDuplicated(gsub("\\.\\d+$", "", names(data_all)))


# (2)
# ─────────────────────────────────────────────
# ─── Merge Dasyatis marmorata → Dasyatis tortonesei ─────
# ─────────────────────────────────────────────

# Check that the columns exist
marmo_col <- names(data_all)[grepl("^Dasyatis marmorata$", names(data_all), ignore.case = TRUE)]
torto_col <- names(data_all)[grepl("^Dasyatis tortonesei$", names(data_all), ignore.case = TRUE)]

if (length(marmo_col) > 0 & length(torto_col) > 0) {
  message("Merging detections of 'Dasyatis marmorata' into 'Dasyatis tortonesei'.")
  
  # Merge numeric values safely (handle NAs)
  data_all[[torto_col]] <- rowSums(data_all[, c(marmo_col, torto_col)], na.rm = TRUE)
  
  # Remove the duplicate / outdated column
  data_all <- data_all %>% select(-all_of(marmo_col))
  
} else if (length(marmo_col) > 0 & length(torto_col) == 0) {
  message("Renaming 'Dasyatis marmorata' → 'Dasyatis tortonesei' (no existing D. tortonesei column).")
  names(data_all)[names(data_all) == marmo_col] <- "Dasyatis tortonesei"
} else {
  message("No 'Dasyatis marmorata' column found — nothing to merge.")
}



# (3) 
# ────────────────────────────────────────────────────────
# ────── From metadata, clean order of pooled sites  ─────
# ────────────────────────────────────────────────────────
# Load
metadatas <- read.csv("./data/raw_data/Mtdt/Med_metadonnees_ADNe - v1.2_2018-2024.csv", sep = ",", dec = ".") %>% 
  rename(sample = spygen_code) %>% 
  dplyr::filter(Tele01 == 1) %>% 
  mutate(pool = gsub("\\.", "_", pool)) %>% 
  mutate(pool = ifelse(grepl("no", pool, ignore.case = TRUE),
                       sample,
                       pool))%>%
  
  # ── Ensure pooled IDs are sorted alphabetically ("croissant")
  mutate(pool = sapply(pool, function(x) {
    # Split the pool name if it contains an underscore separating two SPY codes
    parts <- unlist(strsplit(x, "_"))
    
    # If it's a pooled filter with two parts, sort them alphabetically
    if (length(parts) == 2 && all(grepl("^SPY", parts))) {
      parts <- sort(parts)
      x <- paste(parts, collapse = "_")
    }
    
    # Return the corrected pool name
    return(x)
  }))


# Verify that all pooled filters are correctly ordered
check_pools <- metadatas %>%
  filter(grepl("SPY", pool) & grepl("_", pool)) %>%
  mutate(is_sorted = sapply(pool, function(x) {
    parts <- unlist(strsplit(x, "_"))
    identical(parts, sort(parts))
  }))

table(check_pools$is_sorted)

# ────────────────────────────────────────
# ─── Align data_all with metadata ──────
# ────────────────────────────────────────
data_all <- data_all %>%
  mutate(sample = sapply(sample, function(x) {
    parts <- unlist(strsplit(x, "_"))
    if (length(parts) >= 2 && all(grepl("^SPY", parts))) {
      paste(sort(parts), collapse = "_")
    } else x
  }))


# Check that all filters are also present in metadata
data_all_meta <- data_all %>%
  filter(sample %in% metadatas$pool)
##### all samples present


# (4)
# ─────────────────────────────────────────────
# ────── Remove Angelshark contamination ─────
# ─────────────────────────────────────────────
ls_conta_shark <- metadatas %>% 
  dplyr::filter(stringr::str_detect(comments, "contamination Ange de mer"))  

data_all <- data_all %>% 
  mutate(`Squatina squatina` = ifelse(sample %in% ls_conta_shark$sample, 0, `Squatina squatina`))



# (5)
# ─────────────────────────────────────────────
# ────── Remove empty species columns ─────
# ─────────────────────────────────────────────
# Remove empty species columns (i.e. total sum = 0)
#    We exclude the first columns (metadata/sample info)
species_cols <- setdiff(names(data_all_meta), c("sample", "metric"))
non_empty_species <- names(data_all_meta[, species_cols])[
  colSums(data_all_meta[, species_cols], na.rm = TRUE) != 0
]

data_all_meta <- data_all_meta %>%
  select(any_of(c("sample", "metric", non_empty_species)))

message(length(species_cols) - length(non_empty_species),
        " empty species columns removed.")
#### zero species removed # MARIEKE : 3 spp removed

# [OPTIONAL] Print the names of the species that were removed
setdiff(species_cols, non_empty_species) 



# (6) 
# ─────────────────────────────────────────────
# ────────── Rename species complex   ─────────
# ─────────────────────────────────────────────

# Named vector of old → new names
species_rename_map <- c(
  "C. heterurus_H. speculiger"                   = "Cheilopogon_heterurus.Hirundichthys_speculiger",
  "P_tentacularis_P_zvonimiri"                   = "Parablennius_tentacularis.Parablennius_zvonimiri",
  "Parablennius_incognitus_P_sanguinolentus"     = "Parablennius_incognitus.Parablennius_sanguinolentus",
  "Trachurus_mediterraneus_T_trachurus"          = "Trachurus_mediterraneus.Trachurus_trachurus",
  "Copt. rendalli_Oreo. niloticus"               = "Coptodon_rendalli.Oreochromis_niloticus",
  "T_capelanus_T_minutus"                        = "Trisopterus_capelanus.Trisopterus_minutus",
  "G_biscayensis_G_vulgaris"                     = "Gaidropsarus_biscayensis.Gaidropsarus_vulgaris",    #(Gaidropsarus_biscayensis is junior synonym of Gaidropsarus macrophthalmus (Günther, 1867) ref:  Barros-García, D., A.S. Comesaña, R. Bañón, F. Baldó, C.J. Arronte, E. Froufe and A. de Carlos, 2022. Multilocus species delimitation analyses show junior synonyms and deep‑sea unknown species of genus Gaidropsarus (Teleostei: Gadiformes) in the North Atlantic/Mediterranean Sea area. Mar. Biol. 169(131):1-10. )
  "Labrus_merula_L_viridis"                      = "Labrus_merula.Labrus_viridis",
  "Labrus_merula_L_viridis_C_melanocernus"       = "Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus",
  "A_regius_U_cirrosa"                           = "Argyrosomus_regius.Umbrina_cirrosa",
  "S_chrysotaenia_S_sphyraena"                   = "Sphyraena_chrysotaenia.Sphyraena_sphyraena",
  "C. lucerna_L. dieuzeidei"                     = "Chelidonichthys_lucerna.Lepidotrigla_dieuzeidei",
  "C. obscurus_T. lastoviza"                     = "Chelidonichthys_obscurus.Chelidonichthys_lastoviza",
  "E. gurnardus_T. lyra"                         = "Eutrigla_gurnardus.Trigla_lyra",
  "Spicara_flexuosa_S_smaris"                    = "Spicara_flexuosum.Spicara_smaris",
  "D. dentex_P. auriga_P. pagrus"                = "Dentex_dentex.Pagrus_auriga.Pagrus_pagrus",
  "Raja_asterias_R_clavata_R_polystigma"         = "Raja_asterias.Raja_clavata.Raja_polystigma",
  "Notoscopelus elongatus kroyeri"               = "Notoscopelus_elongatus.Notoscopelus_kroyeri"   # sub-species
)

# Check which of the old names actually exist in data_all
existing_rename <- intersect(names(species_rename_map), names(data_all))

if (length(existing_rename) > 0) {
  message("Found ", length(existing_rename), " species names to rename:")
  print(existing_rename)
  
  # Perform renaming only for those that exist
  for (old_name in existing_rename) {
    new_name <- species_rename_map[[old_name]]
    names(data_all)[names(data_all) == old_name] <- new_name
  }
  
  message("Species renaming completed successfully.")
} else {
  message("ℹ️ None of the listed species names found in data_all — nothing to rename.")
}

# Check renamed columns
names(data_all)[names(data_all) %in% species_rename_map]


# (7) 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ─────────────── Rename species with accepted names from Fishbase (misspelling or synonyms)   ────────────────
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Create a named vector of corrections
species_corrections <- c(
  "Facciolella oxyrhyncha"   = "Facciolella oxyrhynchus",
  "Tripterygion melanurum"    = "Tripterygion melanurus",
  "Symphodus melanocercus"    = "Centrolabrus melanocercus",
  "Pomatomus saltator"        = "Pomatomus saltatrix",
  "Oblada melanura"           = "Oblada melanurus"
)

# Check which names exist in your data_all
existing <- intersect(names(species_corrections), names(data_all))

# Rename columns safely (only if they exist)
for (old_name in existing) {
  new_name <- species_corrections[[old_name]]
  names(data_all)[names(data_all) == old_name] <- new_name
}

# Check that replacements worked
print(names(data_all)[names(data_all) %in% species_corrections])




### LATER : Pleuronectes platessa  → possible Platichthys flesus (check % with reference database because species use to be present in Med historicaly)




# (8) 
# ────────────────────────────────────────────────────
# ────────── Check for out of range species   ─────────
# ────────────────────────────────────────────────────
# Preselect using FishBase datasets 

# Species in the Mediterranean Sea n= 788 (Complete)                                      / download date : 10/21/2025
# from https://www.fishbase.se/trophiceco/FishEcoList.php?ve_code=13
Med_fishbase <- readr::read_delim("./data/raw_data/eDNA/Med_fishbase.csv", 
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)

Med_fishbase <- Med_fishbase %>%
  dplyr::mutate(
    Species = ifelse(Species == "Bathytoshia centroura", "Bathytoshia lata", Species),   # see : https://www.fishbase.se/summary/Bathytoshia-centroura
    Status = ifelse(Species == "Bathytoshia lata", "native", Status),  
    Species = ifelse(Species == "Gobius ophiocephalus", "Zosterisessor ophiocephalus", Species) # Change name of Gobius ophiocephalus with real validated name under WORMS : Zosterisessor ophiocephalus
  )
  
# Add missing species
# add Encrasicholina punctifer (https://doi.org/10.1111/jai.13276)
# add Centrophorus squamosus  (not defined in Med but globally distributed deep sea shark that could be in Med as well because detected)
extra_species <- data.frame(
  Species = c("Encrasicholina punctifer", "Centrophorus squamosus"),
  Status = c("introduced", "questionable")
)
Med_fishbase <- bind_rows(Med_fishbase, extra_species)

# Keep only valid status categories
valid_status <- c("endemic", "native", "introduced")
Med_fishbase_valid <- Med_fishbase %>%
  filter(tolower(Status) %in% valid_status)

# Extract species names to compare
species_in_med <- Med_fishbase_valid$Species

# List of Marine Fishes reported from France n = 763                                       / download date : 10/21/2025
# from https://www.fishbase.se/country/CountryChecklist.php?c_code=250&vhabitat=saltwater&csub_code=&cpresence=present 
France_fishbase <- readr::read_csv("./data/raw_data/eDNA/France_fishbase.csv")
species_in_france <- France_fishbase$Species

# Combine unique species names to have all possible fish species in the area 
all_species <- sort(unique(c(species_in_med, species_in_france)))

# Check for miss-match with data_all 
# Extract detected species (exclude metadata columns)
detected_species <- setdiff(names(data_all), c("sample", "metric"))

# Define multi-species complexes to exclude
species_complex <- c(
  "Cheilopogon_heterurus.Hirundichthys_speculiger",
  "Parablennius_tentacularis.Parablennius_zvonimiri",
  "Parablennius_incognitus.Parablennius_sanguinolentus",
  "Trachurus_mediterraneus.Trachurus_trachurus",
  "Coptodon_rendalli.Oreochromis_niloticus",
  "Trisopterus_capelanus.Trisopterus_minutus",
  "Gaidropsarus_biscayensis.Gaidropsarus_vulgaris",
  "Labrus_merula.Labrus_viridis",
  "Labrus_merula.Labrus_viridis.Centrolabrus_melanocercus",
  "Argyrosomus_regius.Umbrina_cirrosa",
  "Sphyraena_chrysotaenia.Sphyraena_sphyraena",
  "Chelidonichthys_lucerna.Lepidotrigla_dieuzeidei",
  "Chelidonichthys_obscurus.Chelidonichthys_lastoviza",
  "Eutrigla_gurnardus.Trigla_lyra",
  "Spicara_flexuosum.Spicara_smaris",
  "Dentex_dentex.Pagrus_auriga.Pagrus_pagrus",
  "Raja_asterias.Raja_clavata.Raja_polystigma",
  "Notoscopelus_elongatus.Notoscopelus_kroyeri"
)

# Species detected by eDNA but not in reference list
species_not_in_ref <- setdiff(detected_species, all_species)
# Remove complex names (those we want to keep)
species_not_in_ref_clean <- setdiff(species_not_in_ref, species_complex)
##### 13 species 



# For species not in Fishbase datasets check with our home dataset :
out_of_range_species1_2 <- readr::read_csv("./data/raw_data/eDNA/out_of_range_species1.2.csv")


# Keep only rows where Decision == "OUT"
out_species <- out_of_range_species1_2 %>%
  filter(toupper(Decision) == "OUT") %>%
  pull(Species) %>%
  unique()

# Check which "species_not_in_ref_clean" are listed as OUT
species_confirmed_out <- intersect(species_not_in_ref_clean, out_species)

# Species not in ref but NOT listed as OUT (potentially valid or uncertain)
species_not_flagged_out <- setdiff(species_not_in_ref_clean, out_species)

# Summary
cat("\nCross-check summary:\n")
cat("• Total species not in reference list: ", length(species_not_in_ref_clean), "\n")
cat("• Of these, marked as 'OUT' in local dataset: ", length(species_confirmed_out), "\n")
cat("• Remaining (not flagged as OUT): ", length(species_not_flagged_out), "\n")

if (length(species_confirmed_out) > 0) {
  cat("\n⚠️ Species confirmed as OUT:\n")
  print(species_confirmed_out)
}

if (length(species_not_flagged_out) > 0) {
  cat("\nℹ️ Species not in reference but not flagged as OUT:\n")
  print(species_not_flagged_out)
}


# (9)
# ─────────────────────────────────────────────────────────────────────────────────────
# ──────── Final cleanup and export replicate, sequence & presence data frames ────────────────
# ─────────────────────────────────────────────────────────────────────────────────────

# Identify which columns correspond to OUT species
cols_to_remove <- intersect(names(data_all), species_confirmed_out)

# Remove those columns
data_all_final <- data_all %>%
  select(-all_of(cols_to_remove))

if (length(cols_to_remove) > 0) {
  cat("Removed species:\n")
  print(cols_to_remove)
}

# Replace spaces with underscores in their names
species_cols <- setdiff(names(data_all_final), c("sample", "metric"))

names(data_all_final)[names(data_all_final) %in% species_cols] <- 
  gsub(" ", "_", names(data_all_final)[names(data_all_final) %in% species_cols])



# (10)
# ───────────────────────────────────────────────────────────────────────────
# ──────── Export replicate, sequence & presence data frames ────────────────
# ───────────────────────────────────────────────────────────────────────────
# Split by metric type

# Number of PCR replicates
data_rep <- data_all_final %>%
  filter(metric == "nb_rep") %>%
  select(-metric)   # remove metric column

# Number of reads
data_seq <- data_all_final %>%
  filter(metric == "nb_seq") %>%
  select(-metric)   # remove metric column

# Presence/absence (based on nb_rep)
data_pres <- data_rep %>%
  mutate(across(-sample, ~ ifelse(. > 0, 1, 0)))   # convert detections to 1/0


# Save all 3 datasets to CSV
write.table(data_rep,
            "./data/processed_data/eDNA/data_MED_teleo_nb_rep_1824_V1.csv",
            row.names = FALSE, dec = ".", sep = ";")

write.table(data_seq,
            "./data/processed_data/eDNA/data_MED_teleo_nb_seq_1824_V1.csv",
            row.names = FALSE, dec = ".", sep = ";")

write.table(data_pres,
            "./data/processed_data/data_MED_teleo_pres_1824_V1.csv",
            row.names = FALSE, dec = ".", sep = ";")



 # ─── DONE ───────────────────────────────────────────────

