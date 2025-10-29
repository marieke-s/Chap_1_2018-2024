# title: BIODIVMED - eDNA in the Atlantic Ocean - Ecological indicators
# author: Jean-Baptiste Juhel modified by Amandine Avouac & Marie ORBLIN
# date: 2024-05-14

#### Libraries  & workspace ####

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
# library(rfishbase)
library(rsq)
library(scales)
library(stringr)
library(tidyr)
library(tidyverse)
library(viridis)

setwd("/Users/marieorblin/Desktop/JB/INDICATEURS/")


#### Load  & format data ####

# Traits / taxo
traits <- read.csv("/Users/marieorblin/Desktop/JB/INDICATEURS Teleo/species_traits_inferred.csv", sep = ";", header=T)
traits <- traits %>% 
  dplyr::filter(!is.na(Genus)) %>% 
  rename_all(~ str_replace_all(., "\\.", " ")) %>% 
  rename_all(~ str_squish(.)) %>% 
  rename(log_depth_min = `log depth_min`, log_depth_max = `log depth_max`)

# Presence
presence <- read.csv("/Users/marieorblin/Desktop/REMY/selection_presences.csv", header = TRUE, sep = ",")  
names(presence) <- gsub("\\.", " ", names(presence))
presence <- presence %>% dplyr::filter(rowSums(across(where(is.numeric))) != 0) # remove filters with 0 species 

# Metadata
meta <- read.csv("/Users/marieorblin/Desktop/ADNe NOVEMBRE 2024/metadatas/dec24/Med_metadonnees_ADNe - v2.1_2018-2023.csv", sep = ",", dec = ".") %>% 
  rename(code_spygen = spygen_code) %>% 
  dplyr::filter(Tele01 == 1) %>% 
  mutate(pool = gsub("\\.", "_", pool)) %>% 
  mutate(pool = ifelse(pool == "", "No", pool)) %>% 
  mutate(code_spygen = ifelse(pool == "No", code_spygen, pool))  
# Modify the Spygen code to add the pooled filters



#### Biodiversity indicators ####

code_Spygen <- presence %>% pull(code_spygen)
# presence <- presence %>% select(-c("code_spygen"))
presence <- presence %>% column_to_rownames("code_spygen")


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


indicators[,1] <- code_Spygen


#### 1  - Species richness - R ####

indicators[,2] <- rowSums(presence)
summary(indicators[,"R"])

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


#### 3  - Elasmobranch species - Elasmo ####

traits <- traits %>%
  mutate(elasmo = if_else(Class == "Elasmobranchii", 1, 0))


#### 4  - Demerso-pelagic species ratio - DeBRa ####

traits <- traits %>%
  mutate(Benthic = ifelse(DemersPelag %in% c("reef-associated", "demersal", "bathydemersal"), 1, 0)) %>% 
  mutate(DP = ifelse(DemersPelag %in% c("pelagic-oceanic", "pelagic-neritic", "bathypelagic", "pelagic"), 1, 0)) # Create a new Dermerso-Pelagic column which will take the value 1 if in the DemersPelag column it is represented by pelagic-oceanic, pelagic-neritic, bathypelagic. Otherwise it will take the value 0.


summary(is.na(traits$DemersPelag)) # Nombre de valeurs pour la colonne.


#### 5  - Number of species listed VU, EN or CR on the IUCN Red List - RedList ####

traits <- traits %>%
  mutate(RedList = if_else(IUCN_category %in% c("CR", "EN", "VU", "NT"), 1, 0))

summary(is.na(traits$IUCN_category)) # 280 TRUE NA
summary(is.na(traits$RedList)) # 


#### 6  - Large Reef Fish Indicator (max length > 20cm) - LRFI ####

traits <- traits %>%
  mutate(LRFI = if_else(Length >= 20, 1,0)) 

summary(is.na(traits$Length)) # 48
summary(is.na(traits$LRFI)) # 48

#### 7  - Top Predators (piscivore & max length > 50cm) - TopPred ####

traits <- traits %>%
  mutate(TopPred = if_else(Length >= 50 & Troph >= 4, 1,0))

summary(is.na(traits$Troph)) # 0
summary(is.na(traits$TopPred)) # 7

#### 8  - Commercial species - Commercial ####

# A modifier avec la richesse des espèces commerciales

print(traits$Importance)

traits <- traits %>%
  mutate(Commercial = if_else(Importance == "subsistence fisheries" | 
                                Importance == "minor commercial" | 
                                Importance == "of no interest" | 
                                Importance == "of potential interest", 0,
                              if_else(Importance == "commercial" |
                                        Importance == "highly commercial" , 1, NA_real_
                              )))

summary(is.na(traits$Commercial)) # 819 TRUE NA

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

summary(indicators)



#### Build data for next indicators  ####

# list species present in the sample
s_i <- presence %>%
  select(names(.)[colSums(.) > 0])
s_i <- names(s_i)
s_i <- gsub("_", " ", s_i)

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
  column_to_rownames(var = "species") %>%
  mutate_if(sapply(., is.character), as.factor) %>%
  mutate_if(sapply(., is.integer), as.factor)

str(traits_fd)

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



save(sp_dist, file = "functional_distances.Rdata")

### Remove the species that are not in traits from the presence matrix
# and reorder

presence_fd <- presence %>% 
  # select(all_of(s_i)) %>%
  select(all_of(rownames(traits_fd))) %>% 
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
phy <- fishtree_phylogeny(species = s_i)                                         # Requested 236 but only found 183 species, 53 not in tree (Phylo does not contain Elasmobranch)

## check that phylogeny and data have matching names (same species tree / dataset)
# Format the data
presence_phy <- t(presence) 
rownames(presence_phy) <- gsub(" ", "_", rownames(presence_phy))

# Remove from the data the species that are not in the tree
nc <- geiger::name.check(phy, presence_phy)                                     # 53 species not in tree
missing <- gsub("_", " ", nc$data_not_tree, fixed = TRUE)

presence_pd <- presence %>% 
  select(- all_of(missing)) 

# presence is converted to probability vector, summing to 1: 
ps <- presence_pd %>%
  mutate(across(everything(.), ~ as.ProbaVector(.x))) %>%
  as.matrix(.)
colnames(ps) <- gsub(" ", "_", colnames(ps), fixed = TRUE)

## Compute PD for each survey 
Phill <- pbapply(ps, 1, ChaoPD,
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
  rownames_to_column(var = "ID") %>% 
  select(-ID)

# Merge indicators and metadata 
ind_meta <- merge(indicators, meta, by = "code_spygen")
ind_meta <- ind_meta %>% distinct(code_spygen, .keep_all = TRUE)

write.table(indicators, "/Users/marieorblin/Desktop/REMY/selection_indicateurs.csv", row.names = F, dec = ".", sep = ";")


#### PARAGRAPHE AMANDINE > ajout de la colone protection pour le fichier metadatas (manquante) (ajout manuel provisoire avant mise a jour fichier ATLANTIQUE ) ####

# Par type de protection_reglementation : no_fishing ; with_fishing ; in_windfarm
# On prend le même code couleur que celui utilisé par JB pour la dbRDA
# Copié-collé à partir de ce que JB a codé et envoyé par e-mail (14.10.24) 

# Travail à partir du fichier indicateur généré à l'étape précédente contenant : scores indicateurs + métadonnées = "ind_meta"
# - sauf que métadonnées colonne protection_reglementation n'apparait pas dedans
# donc ajout de la colonne protection_reglementation provenant de "metadatas_ATL_AMS"

ind_meta_ord <- ind_meta %>% 
  arrange(code_spygen)

metadatas_ATL_AMS_ord <- metadatas_ATL_AMS %>% 
  arrange(code_spygen) %>%
  # On ne garde que les deux colonnes d'intérêt pour les box_plots
  select(c(code_spygen, protection_reglementation)) %>%
  # On ne garde dans metadatas_ATL_AMS_ord que les code_spygen qui apparaissent dans ind_meta_ord
  filter(code_spygen %in% ind_meta_ord$code_spygen) 

# Il manque deux échantillons dans metadatas_ATL_AMS_ord : SPY234565 et SPY234568
# Echantillons a supprimer des indicateurs car concernent zone de Yeu pour le projet de Pierre parc Yeu-Noirmoutier NectonYN
ind_meta_ord_modified <- ind_meta_ord %>%
  mutate(code_spygen = ifelse(code_spygen == "SPY234565", "toremove", code_spygen),
         code_spygen = ifelse(code_spygen == "SPY234568", "toremove", code_spygen)) %>%
  filter(!code_spygen == "toremove")

# Maintenant qu'on a tous les échantillons et les mêmes entre ind_meta_ord_modified et metadatas_ATL_AMS_ord
# On va ajouter la colonne protection_reglementation au fichier ind_meta_ord_modified
ind_meta_prot <- merge(ind_meta_ord_modified, metadatas_ATL_AMS_ord, by.x = 'code_spygen', by.y = 'code_spygen') %>%
  arrange(protection_reglementation) %>%
  relocate(protection_reglementation, .after = code_spygen)

# Petit check manuel qui a été fait 

# Je génère le fichiers contenant calculs indicateurs + métadonnées + protection_reglementation
write.csv(ind_meta_prot, "... /indicators_meta_protection_yeu-guerande.csv")

# A partir de ce fichier où j'ai toutes les informations nécessaires, je réalise le plot (copié-collé e-mail JB du 14.10.24))
# Par type de protection_reglementation : no_fishing ; with_fishing ; in_windfarm
# On prend le même code couleur que celui utilisé par JB pour la dbRDA


results <- read.csv("indicateurs_yeu_guerande_results_1510/indicators_meta_protection_yeu-guerande.csv")



###### BOXPLOT####
## Richesse spécifique

box_richness <- ggplot(results) +
  aes(x = factor(protection_reglementation ,levels = c("no_fishing", "with_fishing", "in_windfarm"), labels = c("NO FISHING", "WITH FISHING", "IN WINDFARM")), y = R) +
  geom_boxplot(fill = c("royalblue", "#FF6F61", "mediumaquamarine"), color = "black") +
  scale_color_manual(values = c("Y" = "orangered", "N" = "black")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("TOTAL SPECIES RICHNESS")


## ElasmoBranches

box_elasmo <- ggplot(results) +
  aes(x = factor(protection_reglementation ,levels = c("no_fishing", "with_fishing", "in_windfarm"), labels = c("NO FISHING", "WITH FISHING", "IN WINDFARM")), y = Elasmo) +
  geom_boxplot(fill = c("royalblue", "#FF6F61", "mediumaquamarine"), color = "black") +
  scale_color_manual(values = c("Y" = "orangered", "N" = "black")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("TOTAL ELASMOBRANCHS")


## IUCN

box_redlist <- ggplot(results) +
  aes(x = factor(protection_reglementation ,levels = c("no_fishing", "with_fishing", "in_windfarm"),labels = c("NO FISHING", "WITH FISHING", "IN WINDFARM")), y = RedList ) +
  geom_boxplot(fill = c("royalblue", "#FF6F61", "mediumaquamarine"), color = "black") +
  scale_color_manual(values = c("Y" = "orangered", "N" = "black")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("TOTAL IUCN REDLIST")

## SPS Commerciales

box_commercials <- ggplot(results) +
  aes(x = factor(protection_reglementation ,levels = c("no_fishing", "with_fishing", "in_windfarm"),  labels = c("NO FISHING", "WITH FISHING", "IN WINDFARM")), y = Commercial ) +
  geom_boxplot(fill = c("royalblue", "#FF6F61", "mediumaquamarine"), color = "black") +
  scale_color_manual(values = c("Y" = "orangered", "N" = "black")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("TOTAL COMMERCIAL SPECIES")

