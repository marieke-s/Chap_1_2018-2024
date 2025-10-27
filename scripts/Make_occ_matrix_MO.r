#### Bibliothèques et espace de travail ####

# Clean env
rm(list = ls())

library(ape)
library(dplyr)
library(entropart)
library(fastDummies)
library(fishtree)
library(geiger)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(margins)
library(mFD)
library(pbapply)
library(picante)
library(purrr)
library(rfishbase)
library(rsq)
library(scales)
library(stringr)
library(tidyr)
library(tidyverse)
library(viridis)
library(readxl)

library(readxl)
library(dplyr)
library(tibble)


#1) Charger la donnée ####

#Charger la checklist Med - Atlantique
check <- read.csv("./data/raw_data/eDNA/checklist_Med_Atl_cleaned.csv", sep = ";")
`%!in%` = Negate(`%in%`)  # Créer un opérateur pour la négation de `%in%`


#Lire le fichier brut SPYGEN sans noms de colonnes (pour travailler sur les lignes)
raw <- readxl::read_xlsx("./data/raw_data/eDNA/Results_Teleo_2018-2024_with_IPOCOM.xlsx",
                 col_names = FALSE)

#trouver la ligne qui contient "scientific_name" 
header_row <- which(apply(raw, 1, function(r) any(grepl("^scientific_name$",
                                                        r, ignore.case = TRUE))))
if(length(header_row) == 0) stop("Impossible de trouver la ligne contenant 'scientific_name'.")

#Indices des lignes qui correspondent aux taxons (juste après la ligne header)
species_rows <- (header_row + 1):nrow(raw)

#Prendre la colonne des noms scientifiques (4ème colonne dans ton fichier)
sp_col <- as.character(raw[[4]][species_rows])

#Détecter les "mauvaises" lignes PARMI ces seules lignes d'espèces
is_sp_word     <- grepl("\\bsp\\.?\\b", sp_col, ignore.case = TRUE)  # sp ou sp.
is_underscore  <- grepl("_", sp_col, fixed = TRUE)                   # contient "_"
is_family      <- grepl("dae", sp_col, ignore.case = TRUE)           # contient "dae"

bad_mask <- is_sp_word | is_underscore | is_family
bad_species_rows <- species_rows[bad_mask]

#(Optionnel) afficher les valeurs qui seront supprimées
if(length(bad_species_rows) > 0){
  message("Taxons détectés comme à supprimer (exemples) :")
  print(unique(sp_col[bad_mask]))
} else {
  message("Aucun taxon problématique détecté parmi les lignes d'espèces.")
}

#créer un raw_clean sans ces lignes d'espèces problématiques
raw_clean <- raw[-bad_species_rows, ]


#2) Formater le dataframe ####
data_all <- data.frame(t(raw_clean))                                                     # transposes the data frame  
data_all <- tibble::rownames_to_column(data_all, var = "region")                                # add the row names as a new column 
names(data_all) <- data_all[4,]                                                         # replace the Column names by the species names
data_all <- data_all[,-1:-5] # colonnes
names(data_all)[1] <- "code_spygen"
data_all <- data_all[-1:-5,]  #lignes
names(data_all)[2] <- "variable"

# Remplacer les NA par 0
data_all[is.na(data_all)] <- 0

# delete the firts 5 lines (adjust for each dataframe to process)
data_all[,3:length(data_all)] <- apply(data_all[,3:length(data_all)], 2, 
                                       function(x)gsub('\\s+', '',x))                    # remove blank spaces in excess
data_all[,3:length(data_all)] <- sapply(data_all[,3:length(data_all)], as.numeric)              # make sure the data are numeric


#----- Check PCR replicates numbers [Marieke 27/10/2025] -----
pcr_rep <- data_all %>%
  dplyr::filter(variable == "nb_rep") 

# Print max values per species column
colnames(pcr_rep)
sapply(pcr_rep[,3:ncol(pcr_rep)], max)  

# Select rows with values > 12
pcr_above_12 <- pcr_rep %>%
  filter(if_any(-c(code_spygen, variable), ~ . > 12))

# Compare with mtdt_3_pcr_above_12 from 3a_Prep_eDNA_data.R
mtdt_3_pcr_above_12 <- read.csv("./data/processed_data/eDNA/mtdt_3_pcr_above_12.csv", sep = ",")

# samples in pcr_above_12 (raw data) not in mtdt_3_pcr_above_12 (matrix data)
length(setdiff(pcr_above_12$code_spygen, mtdt_3_pcr_above_12$pool))  

# samples in mtdt_3_pcr_above_12 (matrix data) not in pcr_above_12 (raw data)
setdiff(mtdt_3_pcr_above_12$pool, pcr_above_12$code_spygen)

length(pcr_above_12$code_spygen)

#3) CLEAN LA DATA (1.)#### 
# Correct the identification of 'Dasyatis marmorata' into 'Dasyatis tortonesei'
which(grepl("Dasyatis marmorata", names(data_all))) 
which(grepl("Dasyatis tortonesei", names(data_all))) 

data_all[length(data_all)+1] <- data_all[,278] + data_all[,280]

data_all <- data_all %>%  
  dplyr:: select(! c("Dasyatis marmorata", "Dasyatis tortonesei")) %>%                  # remove the columns with the duplicated species
  rename("Dasyatis tortonesei" = "V300" )                                       # rename the merged column

which(grepl("Notoscopelus elongatus", names(data_all))) 
which(grepl("Dasyatis SP", names(data_all))) 





#4) Matching w/ Metadatas ####
# Load metadata
metadatas <- read.csv("Desktop/ADNe_MATRICES_sept25/Med_metadonnees_ADNe - v1.2_2018-2024.csv", sep = ",", dec = ".") %>% 
  rename(code_spygen = spygen_code) %>% 
  dplyr::filter(Tele01 == 1) %>% 
  mutate(pool = gsub("\\.", "_", pool)) %>% 
  mutate(pool = ifelse(pool == "no", code_spygen, pool)) %>%                 # Dans la colonne "pool" remplacer les valeurs "no" par le "code_spygen" (utile pour matching avec data_all par la suite)

  mutate(pool_inv = case_when(nchar(pool) > 10 ~ 
                                paste(substr(pool,11,19),
                                      substr(pool, 1, 9), sep = "-"), 
                              TRUE ~ pool)) 

# Remove Angelshark contamination
ls_conta_shark <- metadatas %>% 
  dplyr::filter(str_detect(comments, "contamination Ange de mer"))  

data_all <- data_all %>% 
  mutate(`Squatina squatina` = ifelse(code_spygen %in% ls_conta_shark$code_spygen, 0, `Squatina squatina`))

# Remove filters absent from the metadata and remove empty species columns
data_all_meta <- data_all %>% 
  mutate(code_spygen = gsub("\\-", "_", code_spygen)) %>% 
  dplyr::filter(code_spygen %in% metadatas$pool) %>%
  mutate(across(-c(1:2), ~ if (sum(.) != 0) . else NULL)) %>%                   # check for empty columns # 104 species removed
  dplyr::select(where(~ !is.null(.)))


# List of filters present in the results but absent from metadata
data_out_meta <- subset(data_all, !(data_all$code_spygen %in% metadatas$pool))  
filter_out_meta <- data_out_meta %>%                                        
  pull(code_spygen) %>% 
  unique()


# List of filters in the metadata but absent from the results
meta_out_data <- subset(metadatas, !(metadatas$pool %in% data_all$code_spygen))  
filter_out_data <- meta_out_data %>% 
  pull(code_spygen) %>% 
  unique()



# ----------------------- CLEAN LA DATA 2. -------------------------
#5)  DOUBLONS : Identifier toutes les colonnes se terminant par ".1"  ####
cols_dot1 <- grep("\\.1$", names(data_all_meta), value = TRUE)

#Déduire les noms "de base" (sans le .1)
base_names <- sub("\\.1$", "", cols_dot1)

# Boucle sur chaque espèce pour fusionner les colonnes
for (species in base_names) {
  
  # Vérifier que la colonne sans .1 existe bien
  if (species %in% names(data_all_meta)) {
    
    # Créer une nouvelle colonne fusionnée = somme des deux
    data_all_meta[[species]] <- data_all_meta[[species]] + data_all_meta[[paste0(species, ".1")]]
    
    # Supprimer la colonne doublon
    data_all_meta <- data_all_meta %>%
      select(!all_of(paste0(species, ".1")))
    
  } else {
    message("Colonne ", species, " n’a pas de doublon clair, sautée.")
  }
}


#6)  Matching avec la checklist MED/ATL ####
# Check for species absent from the checklist
s_i <- data_all_meta %>%
  dplyr::select(-code_spygen, -variable) %>% 
  dplyr::select(names(.)[colSums(.) > 0])
s_i <- names(s_i)
s_i <- data.frame(s_i)

sp_in_checklist <- s_i %>% dplyr::filter(s_i %in% check$species) 
sp_out_checklist <- s_i %>% dplyr::filter(s_i %!in% check$species) 


#REMOVE INCORRECT SP NAMES
sp_out_notOK <- sp_out_checklist %>%
  filter(!grepl("\\s", s_i))  # Filtre les noms qui n'ont pas d'espace, donc composés d'un seul mot

# Afficher ou sauvegarder l'objet
sp_out_notOK #plot

sp_out_OK <- sp_out_checklist %>%
  filter(grepl("\\s", s_i)) 

sp_out_OK #plot


#6.2) ou : les identifier d'abord > faire le merged des colonnes doublons une par une ####
## METHODE 2 ## Lister les noms d'espèces avec .1

#sps_dot1 <- grep("\\.1$", sp_out_OK$s_i)
#sp_out_OK[sps_dot1, ]

### Merged les colonnes x + x.1 
#Tripterygion tripteronotum
#which(grepl("Tripterygion tripteronotum", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,44] + data_all_meta[,45]
#data_all_meta <- data_all_meta %>%  
 # dplyr::select(! c("Tripterygion tripteronotum", "Tripterygion tripteronotum.1")) %>%                  # remove the columns with the duplicated species
  #rename("Tripterygion tripteronotum" = "V299" )                     

#Trachinotus ovatus
#which(grepl("Trachinotus ovatus", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,48] + data_all_meta[,49]
#data_all_meta <- data_all_meta %>%  
 # dplyr::select(! c("Trachinotus ovatus", "Trachinotus ovatus.1")) %>%                  
 # rename("Trachinotus ovatus" = "V298" )                     

#Gobius cobitis
#which(grepl("Gobius cobitis", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,82] + data_all_meta[,83]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Gobius cobitis", "Gobius cobitis.1")) %>%                  
  #rename("Gobius cobitis" = "V297" )

#Pomatoschistus minutus
#which(grepl("Pomatoschistus minutus", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,95] + data_all_meta[,96]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Pomatoschistus minutus", "Pomatoschistus minutus.1")) %>%                 
  #rename("Pomatoschistus minutus" = "V296" )

#Ceratoscopelus maderensis
#which(grepl("Ceratoscopelus maderensis", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,130] + data_all_meta[,131]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Ceratoscopelus maderensis", "Ceratoscopelus maderensis.1")) %>%                  
  #rename("Ceratoscopelus maderensis" = "V295" )

#Dicentrarchus labrax
#which(grepl("Dicentrarchus labrax", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,138] + data_all_meta[,139]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Dicentrarchus labrax", "Dicentrarchus labrax.1")) %>%                  
  #rename("Dicentrarchus labrax" = "V294" )

#Syngnathus acus
#which(grepl("Syngnathus acus", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,244] + data_all_meta[,245]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Syngnathus acus", "Syngnathus acus.1")) %>%                  
 # rename("Syngnathus acus" = "V293" )

#Torpedo marmorata
#which(grepl("Torpedo marmorata", names(data_all_meta))) 
#data_all_meta[length(data_all_meta)+1] <- data_all_meta[,282] + data_all_meta[,283]
#data_all_meta <- data_all_meta %>%  
  #dplyr::select(! c("Torpedo marmorata", "Torpedo marmorata.1")) %>%                  
  #rename("Torpedo marmorata" = "V292" )



#Identifier toutes les colonnes se terminant par ".1"
cols_dot1 <- grep("\\.1$", names(data_all_meta2), value = TRUE)

#Déduire les noms "de base" (sans le .1)
base_names <- sub("\\.1$", "", cols_dot1)

# Boucle sur chaque espèce pour fusionner les colonnes
for (species in base_names) {
  
  # Vérifier que la colonne sans .1 existe bien
  if (species %in% names(data_all_meta2)) {
    
    # Créer une nouvelle colonne fusionnée = somme des deux
    data_all_meta2[[species]] <- data_all_meta2[[species]] + data_all_meta2[[paste0(species, ".1")]]
    
    # Supprimer la colonne doublon
    data_all_meta2 <- data_all_meta2 %>%
      select(!all_of(paste0(species, ".1")))
    
  } else {
    message("Colonne ", species, " n’a pas de doublon clair, sautée.")
  }
}





#7) Suppression des colonnes de data_all_meta qui ont des noms présents dans noms_invalides ####
data_all_meta_clean <- data_all_meta %>%
  dplyr::select(-all_of(sp_out_notOK$s_i))

# trier les colonnes par ordre alphabétique
data_clean <- data_all_meta_clean %>%
  dplyr::select(1:2, sort(names(.)[3:ncol(.)]))

# supprimer les espèces exogènes détectées par Alice (SPYGEN) + sp_out_OK après checking
data_clean <- data_clean %>%
  select(-c(
    "Ophisurus macrorhynchos",
    "Nansenia boreacrassicauda",
    "Cololabis saira",
    "Encrasicholina punctifer",
    "Lophius litulon",
    "Planiliza macrolepis",
    "Diaphus anderseni",
    "Johnius belangerii",
    "Larimichthys crocea",
    "Cataetyx rubrirostris",
    "Stomias boa",
    "Centrophorus squamosus",
    "Gobiusculus flavescens",
    "Bolinichthys supralateralis",
    "Notacanthus chemnitzii",
    "Euthynnus affinis",
    "Vinciguerria nimbaria",
    "Hyperoplus lanceolatus"
  ))

# supprimer les espèces exogènes restantes présentées par Celia 
data_clean <- data_clean %>%
  select(-c(
    "Istiophorus platypterus" ))

# renommer Pomatomus saltator selon nom correct > WORMS
names(data_clean)[names(data_clean) == "pomatomus saltatrix"] <- "Pomatomus saltatrix"


#8) Create 3 data frames : number of PCR replicates, number of reads, presence/absence of species ####

data_rep <- data_clean %>% 
  dplyr::filter(variable == "nb_rep") #%>% 
#distinct(code_spygen, keep_all = TRUE)


data_seq <- data_clean %>% 
  dplyr::filter(variable == "nb_seq") #%>% 
#distinct(code_spygen, keep_all = TRUE)

data_pres <- data_rep %>% 
  dplyr:: select(-variable ) %>% 
  mutate_at(vars(-code_spygen), ~ifelse(. > 0, 1, .)) #%>%                        # replace values by 0/1
# distinct(code_spygen, keep_all = TRUE)                           

write.table(data_pres, "Desktop/ADNe_MATRICES_sept25/data_MED_teleo_pres_1824_V.csv", row.names = F, dec = ".", sep = ";")
write.table(data_rep, "Desktop/ADNe_MATRICES_sept25/data_MED_teleo_rep_1824_V1.csv", row.names = F, dec = ".", sep = ";")
write.table(data_seq, "Desktop/ADNe_MATRICES_sept25/data_MED_teleo_seq_1824_V1.csv", row.names = F, dec = ".", sep = ";")

