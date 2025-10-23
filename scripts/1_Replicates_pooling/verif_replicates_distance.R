# Script pour vérifier les réplicates (distance entre les deux coord géog)
setwd("/Users/amandineavouac/Documents/Metadatas_Med_2024")
library(tidyverse)
library(sf)

# Chargement des métadonnées
metadata_med_2018_2024 <- read.csv("/Users/amandineavouac/Documents/Metadatas_Med_2024/datas/Med_metadonnees_ADNe_2018_2024_check_component.csv", sep = ";", dec =",") %>%
  select(-X) %>%
  mutate(pool = ifelse(pool == "SPY211186_SPY211185", "SPY211185_SPY211186", pool)) %>%
  mutate(replicates = ifelse(replicates == "SPY211265_SPY223409_SPY211270_SPY211271", "SPY211265_SPY223409/SPY211270_SPY211271", replicates))

# Passer les pools dans la colonne spygen_code lorsque pool != "no"
metadata_rep1 <- metadata_med_2018_2024 %>%
  mutate(spygen_code = case_when(pool != "no" ~ pool, TRUE ~ spygen_code)) %>%
  
# Supprimer les doublons (car un même pool apparaît sur deux lignes = on ne veut qu'une ligne par pool)
  distinct(spygen_code, .keep_all = TRUE) %>%
# .keep_all = TRUE pour garder les autres colonnes (variables) sinon on n'aurait plus que la colonne spygen_code dans le dt

# Supprimer replicates == "no" car on ne vérifie les distance entres les coord géog que pour les replicates (logique :))
  filter(replicates != "no")

# Création d'une variable replicates (sans doublons)
replicates <- metadata_rep1 %>%
  distinct(replicates) %>%
  unlist()

# Création d'une boucle qui va tourner encore et encore tant que toutes les combinaisons de coord géog n'auront pas été vérifiées et mises dans nouvelle matrice
dist_rep_min = dist_rep_max = vector()

for (i in 1:length(replicates)) {
# pour faire tourner la boucle i de 1 jusqu'à la fin de replicates (on aurait pu écrire length = 769)
# j'assigne à rep les valeurs SPY apparaissant dans replicates on les dissociant en fonction du "/" qui les sépare
  rep <- str_split(replicates[i], "/" ) %>%
    unlist()

  
  
# LES STARTS  
mtd_temp_start <- metadata_rep1 %>%
# j'assigne à mtd_temp_start (en ne gardant que les spygen_code qui apparaissent dans metadata_rep1)
# st_as_sf pour convertir un objet R classique (comme df) en objet sf (simple features) permettant la manip de données géospatiales
  filter(spygen_code %in% rep) %>%
  st_as_sf(coords = c("latitude_start_DD", "longitude_start_DD"), crs = "+proj=longlat +datum=WGS84 +no_defs", remove = FALSE) 
# crs = permet de définir le système de référence de coord géog en degrés (longitude/latitude) et en WGS84 
# remove = FALSE indique que les colonnes coord (latitude_start_DD et longitude_start_DD) ne seront pas supprimées après  conversion en sf



# LES ENDS  
# ajouter une condition en if : si j'ai du NA dans les ends alors ne fait pas la boucle avec les end
if (sum(is.na(mtd_temp_start$latitude_end_DD) | is.na(mtd_temp_start$longitude_end_DD))>0) {
# Si NA dans latitude_end_DD ou dans longitude_end_DD alors ignorer cette ligne et passer à i suivant
  #next
  #passage au i suivant
  
  #d1 <- st_distance(mtd_temp_start[1,], mtd_temp_start[2,])
  
  d <- st_distance(mtd_temp_start) %>%
    replace(., col(.) == row(.), NA)
  
  dist_min <- min(d, na.rm=TRUE)
  dist_max <- max(d, na.rm=TRUE)
  dist_rep_min[i] <- dist_min 
  dist_rep_max[i] <- dist_max
  
} else {
  
  mtd_temp_end <- metadata_rep1 %>%
    filter(spygen_code %in% rep) %>%
    st_as_sf( coords = c("latitude_end_DD", "longitude_end_DD"), crs = "+proj=longlat +datum=WGS84 +no_defs", remove = FALSE) 
 
  #d1 <- st_distance(mtd_temp_start[1,], mtd_temp_start[2,])
  #d2 <- st_distance(mtd_temp_start[1,], mtd_temp_end[2,])
  #d3 <- st_distance(mtd_temp_start[2,], mtd_temp_end[1,])
  #d4 <- st_distance(mtd_temp_end[1,], mtd_temp_end[2,])
  
  #dist <- min(d1,d2,d3,d4)
  mtd_temp <- rbind(mtd_temp_start, mtd_temp_end)
  
  d <- st_distance(mtd_temp) %>%
    replace(., col(.) == row(.), NA)
  
  dist_min <- min(d, na.rm=TRUE)
  dist_max <- max(d, na.rm=TRUE)
  dist_rep_min[i] <- dist_min 
  dist_rep_max[i] <- dist_max
}

# J'assigne à d1 d2 d3 d4 les distance entre les start/end/1/2/3/4 selon le nbre de spygen_code intervenant dans rep (replicates)
  
}


df <- cbind(as.vector(replicates), dist_rep_min, dist_rep_max)
write.csv(df, "datas/dist_min_max.csv")

# Ajout nouvelle colonne pour identifier les spygen_code où on a besoin de vérifier les distances
df2 <- data.frame(df) %>%
  add_column(., "to_check") %>%
  rename(to_check = '"to_check"') %>%
  mutate(to_check = ifelse(to_check == to_check, NA, to_check))

# Je convertis en numérique les colonnes sur lesquelles nous travaillons
df2$dist_rep_max <- as.numeric(df2$dist_rep_max)
df2$dist_rep_min <- as.numeric(df2$dist_rep_min)

# Je mets en évidence "to_check" les colonnes où dist_rep_max >= 1000m (=1km)
# Pour les vérifier par la suite
df3 <- df2 %>%
  mutate(to_check = ifelse(dist_rep_max >= 1000, "to_check", to_check)) %>%
  rename(replicates = V1)

# J'ajoute les métadonnées (site, date, project) pour avoir un peu plus d'infos sur les éch concernés à vérifier
metadata_med_2018_2024_2 <- metadata_med_2018_2024 %>%
  select(replicates, site, subsite, date, project, latitude_start_DD, longitude_start_DD, latitude_end_DD, longitude_end_DD) %>%
  distinct(replicates, .keep_all = TRUE)
# préciser distinct 
# Je fusionne les métadonnées d'intérêt avec le df3
df4 <- merge(df3, metadata_med_2018_2024_2, by = "replicates")
# Le problème c'est que je passe de 769 obs (ou replicates) de df3 à 814 en df4 et je ne comprends pas pq 

# jusqu'à 2 km d'écart entre start et end c'est bon 
# au-delà : to_check, partir des + distants pour faire verif manuelle Google Earth (voir si tjrs le même projet qui revient par ex.) 
# voir si QGIS possible (avec table d'attribut pour avoir direct info du projet)

write.csv(df4, "datas/df4.csv")

df4 <- df3 %>%
  left_join(metadata_med_2018_2024_2, by = "replicates")


left_join(metadata_med_2018_2024_2, by = c("replicate" = "replicates"))

# Essais : (à supprimer par la suite car infructueux)
# dans df3 quels replicates de df4 sont absents ?
test <- df4 %>%
  filter(replicates %in% df3$replicates)

# dans df4 quels réplicates de df3 sont absents ?
test <- df3 %>%
  filter(replicates %in% df4$replicates)
  
# = dans FDT_17_Corse quels échantillons de Online_Corse sont absents ?
  Missing_Corsica_in_FDT_17  <- Online_Corse %>%
  filter(!spygen_code %in% FDT_17_Corse$spygen_code) 


df4 <- df3 %>%
  inner_join(metadata_med_2018_2024_2, by = "replicates")
  
df4 <- df3 %>%
  left_join(metadata_med_2018_2024_2, by = "replicates")




df4 <- metadata_med_2018_2024_2 %>%
  left_join(df3, by = "replicates")
  





