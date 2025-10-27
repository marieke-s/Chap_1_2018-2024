# Bout de script pour Marieke - métadonnées Med 

#### Ajout colonne replicates ####
# colonne replicates pour indiquer quel code SPY a été prélevé en parallèle d'un autre
# (sachant que pour 2018-2023 il y avait parfois 4 réplicats/site car deux bateaux en même temps et deux filtrations en même temps en aller + en retour sur chq bateau)
# reprendre 2023 nomenclature replicates : séparer par des / (avant 2023 : deux pompes = 4 réplicats pour un même site )
raw_10 <- raw_8 %>%
  add_column(., "replicates") %>%
  rename(replicates = '"replicates"') %>%
  mutate(replicates = ifelse(replicates == replicates, NA, replicates)) %>%
  mutate(replicates = ifelse(spygen_code == "SPY214171", "SPY214171/SPY214175", replicates),
         replicates = ifelse(spygen_code == "SPY214175", "SPY214171/SPY214175", replicates),
         replicates = ifelse(spygen_code == "SPY231866", "SPY231866/SPY231867", replicates),
         replicates = ifelse(spygen_code == "SPY231867", "SPY231866/SPY231867", replicates),
         replicates = ifelse(spygen_code == "SPY231868", "SPY231868/SPY231876", replicates),
         replicates = ifelse(spygen_code == "SPY231876", "SPY231868/SPY231876", replicates),
         replicates = ifelse(spygen_code == "SPY2401065", "SPY2401065/SPY2401066", replicates),
         replicates = ifelse(spygen_code == "SPY2401066", "SPY2401065/SPY2401066", replicates),
         replicates = ifelse(spygen_code == "SPY2401034", "SPY2401034/SPY2402614", replicates),
         replicates = ifelse(spygen_code == "SPY2402614", "SPY2401034/SPY2402614", replicates),
         replicates = ifelse(spygen_code == "SPY2402598", "SPY2402598/SPY2402639", replicates),
         replicates = ifelse(spygen_code == "SPY2402639", "SPY2402598/SPY2402639", replicates),
         replicates = ifelse(spygen_code == "SPY2403584", "SPY2403584/SPY2403589", replicates),
         replicates = ifelse(spygen_code == "SPY2403589", "SPY2403584/SPY2403589", replicates),
         replicates = ifelse(spygen_code == "SPY2403574", "SPY2403574/SPY2403591", replicates),
         replicates = ifelse(spygen_code == "SPY2403591", "SPY2403574/SPY2403591", replicates),
         replicates = ifelse(spygen_code == "SPY2401063", "SPY2401063/SPY2401104", replicates),
         replicates = ifelse(spygen_code == "SPY2401104", "SPY2401063/SPY2401104", replicates),
         replicates = ifelse(spygen_code == "SPY2401064", "SPY2401064/SPY2401078", replicates),
         replicates = ifelse(spygen_code == "SPY2401078", "SPY2401064/SPY2401078", replicates)) %>%
  
# Là je me rends compte qu'à la main c'est source d'erreur et promet d'être EXTREMEMENT long 
# Alors je fais appel à ChatGPT, après plusieurs tentatives car il combinait des réplicats même lorsque NA dans time_start... j'obtiens le code ci-dessous : 
# Regroupement par site, date, time_start, et duration
  group_by(site, date, time_start, duration) %>%
# Créer la colonne replication_group contenant les spygen_code triés
  mutate(replication_group = list(sort(spygen_code))) %>%
# Créer la colonne replicates : on vérifie si time_start est NA, si oui on met NA, sinon on combine les spygen_code
  mutate(replicates = if_else(is.na(time_start), 
                              NA_character_,  # Si time_start est NA, on met NA dans replicates
                              if_else(length(replication_group[[1]]) > 1,
                                      paste(replication_group[[1]][1], replication_group[[1]][2], sep = "/"),
                                      NA_character_))) %>%
# Supprimer la colonne temporaire replication_group
  select(-replication_group) %>%
  relocate(replicates, .after = spygen_code) %>%
# Maintenant je dois aller chercher les replicates pour tous les spygen_code où is.na(replicates)
# Créer une colonne 'replication_group' pour grouper par site, date, subsite
  group_by(site, date, subsite) %>%
# Mettre à jour 'replicates' si is.na(replicates) et respecter les conditions
# S'applique pour tous les sites sauf marseille et pelagos car ne peuvent être combinés par rapport aux date / site / subsite
  mutate(replicates = if_else(is.na(replicates) & site != "marseille" & site != "pelagos" & n_distinct(spygen_code) > 1,  
# Concaténer les spygen_codes dans l'ordre croissant
  paste(sort(spygen_code)[1], sort(spygen_code)[2], sep = "/"), replicates)) %>%  # Si la condition n'est pas remplie, garder l'ancien replicates
  
# Supprimer les groupes pour que le dataframe soit "dégrouper"
  ungroup()
# Pour contrôler : raw_10 %>% filter(is.na(replicates)) %>% group_by(spygen_code, replicates, site, date, time_start, duration) %>% arrange(site) %>% count() %>% print(n=500)

#### Vérification si 'ungroup' est vraiment nécessaire ####
# Vérification si 'ungroup' est vraiment nécessaire
if (n_groups(raw_10) > 0) {
  raw_10 <- raw_10 %>%
    ungroup()  # Appliquer ungroup uniquement si le dataframe est encore groupé
}



#### Modifications pour site == "palermo" et site == "marseille"car réplicats incorrects (duration = 15min, 4 fois time_start = 08:17) ####
# : raw_10 %>% filter(site == "palermo") %>% group_by(spygen_code, pool, replicates, site, duration, time_start) %>% count() %>% arrange(time_start) %>% print(n=500)
raw_11 <- raw_10 %>%
#### idem pour les sites == "marseille", il faut utiliser infos de "pool" pour compléter "replicates" ####
mutate(replicates = case_when(site == "palermo" & pool != "no" ~ pool, TRUE ~ replicates)) %>%
# modifie la colonne replicates en y mettant la valeur de "pool" lorsque site == "palermo" et pool différent de "no" 
# (pour pas modifier les replicates quand pool = no sinon ça met le non dans replicates et ça enlève les replicats qui étaient corrects)
# là j'ai récupéré les infos des pools pour mettre dans replicates
  mutate(replicates = case_when(site == "marseille" & pool != "no" ~ pool, TRUE ~ replicates)) %>%
  mutate(replicates = str_replace_all(replicates, "\\_", "/")) %>%
# uniquement pour les site == "palermo" (car pour site == "marseille ça ne peut pas marcher car date/duration pas complété)
# j'ajoute des replicates les pools (que j'ai déjà) avec AUSSI les autres code_spygen réplicats du site/date/duration !
# = 4 code_spygen pour certains éch (2 x 2 pools pour 1 site car 4 réplicats de 15min par exemple = 2 code_spy de 15min + 2 autres code_spy de 15min)
# uniquement pour les site == "palermo" (car pour site == "marseille ça ne peut pas marcher car date/duration pas complété)
  filter(site == "palermo") %>%  
  group_by(time_start, duration) %>%  # Regrouper par time_start et duration
  mutate(replicates = if_else(duration == "15" & n() > 1,  
                              paste(sort(spygen_code), collapse = "/"),
                              replicates)) %>% 
  ungroup()  %>%
  bind_rows(filter(raw_10, site != "palermo"))

# Je dois finir de compléter les replicates pour site == marseille (deepheart) à la main car dépend de la fiche terrain envoyée à Spygen

# Réplicats qui vont ensemble : 
# DH 22/01/2024
# OK : Moyades 20m SPY232659 et SPY232668 c tout
# OK : Imperial 20m SPY232670 et SPY232661 et SPY232666 et SPY232652
# OK : Moyades 40m SPY232669 et SPY232651
# OK : Imperial 80m SPY232658 et SPY232667 et SPY232643 et SPY232649
raw_12 <- raw_11 %>%
  mutate(replicates = case_when(spygen_code == "SPY232659" | spygen_code == "SPY232668" ~ "SPY232659/SPY232668", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY232670" | spygen_code == "SPY232661" | spygen_code == "SPY232666" | spygen_code == "SPY232652" ~ "SPY232652/SPY232661/SPY232666/SPY232670", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY232669" | spygen_code == "SPY232651" ~ "SPY232651/SPY232669", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY232658" | spygen_code == "SPY232667" | spygen_code == "SPY232643" | spygen_code == "SPY232649" ~ "SPY232643/SPY232649/SPY232658/SPY232667", TRUE ~ replicates)) %>%
# Modification de certains pools DH qui étaient incorrects et se sont mélanges avec les replicates
  mutate(pool = case_when(spygen_code == "SPY232659" | spygen_code == "SPY232668" ~ "SPY232659_SPY232668", TRUE ~ pool),
         pool = case_when(spygen_code == "SPY232669" | spygen_code == "SPY232651" ~ "SPY232651_SPY232669", TRUE ~ pool),
         pool = case_when(spygen_code == "SPY232670" | spygen_code == "SPY232661" ~ "SPY232661_SPY232670", TRUE ~ pool),
         pool = case_when(spygen_code == "SPY232666" | spygen_code == "SPY232652" ~ "SPY232652_SPY232666", TRUE ~ pool),
         pool = case_when(spygen_code == "SPY232658" | spygen_code == "SPY232667" ~ "SPY232658_SPY232667", TRUE ~ pool),
         pool = case_when(spygen_code == "SPY232643" | spygen_code == "SPY232649" ~ "SPY232643_SPY232649", TRUE ~ pool)) %>%
# Vérifications des replicates pour site == "pelagos" = filtrations de 60 min qui ne se suivent pas (time_start) = pas de réplicats pour un même endroit
# raw_12 %>% filter(site == "pelagos") %>% group_by(spygen_code, pool, replicates, date, duration, time_start) %>% count() %>% arrange(date, duration) 
  mutate(replicates = ifelse(site == "pelagos", "no", replicates)) %>%
# Complétion d'infos pour les éch DH du 22/01/24 car infos manquantes dans les métadonnées
  mutate(duration = case_when(date == "2024-01-22" & project == "DeepHeart" ~ "15", TRUE ~ duration),
         estimated_volume = case_when(date == "2024-01-22" & project == "DeepHeart" ~ "15", TRUE ~ estimated_volume)) %>%
  
# Données DeepHeart dates pr lesquels il y a des NA dans replicates : 
# 2024-06-20 : pools et replicates ok dans métadonnées
# 2024-09-30 : pools et replicates ok dans métadonnées
# 2024-10-01 : pools et replicates ok dans métadonnées
# mettre les replicats ensemble 
# Replicates
  mutate(replicates = case_when(spygen_code == "SPY2401239" | spygen_code == "SPY2401255" | spygen_code == "SPY2401226" | spygen_code == "SPY2401251" ~ "SPY2401226/SPY2401239/SPY2401251/SPY2401255", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2401236" | spygen_code == "SPY2401238" | spygen_code == "SPY2401248" | spygen_code == "SPY2401230" ~ "SPY2401230/SPY2401236/SPY2401238/SPY2401248", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2401237" | spygen_code == "SPY2401241" | spygen_code == "SPY2401232" | spygen_code == "SPY2401254" ~ "SPY2401232/SPY2401237/SPY2401241/SPY2401254", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2401253" | spygen_code == "SPY2401245" | spygen_code == "SPY2401244" | spygen_code == "SPY2401252" ~ "SPY2401244/SPY2401245/SPY2401252/SPY2401253", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2402767" | spygen_code == "SPY2401148" | spygen_code == "SPY2401163" | spygen_code == "SPY2401168" ~ "SPY2401148/SPY2401163/SPY2401168/SPY2402767", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2401185" | spygen_code == "SPY2402759" | spygen_code == "SPY2401249"  | spygen_code == "SPY2401191" ~ "SPY2401185/SPY2401191/SPY2401249/SPY2402759", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2404178" | spygen_code == "SPY2404174" | spygen_code == "SPY2404202" | spygen_code == "SPY2404181" ~ "SPY2404174/SPY2404178/SPY2404181/SPY2404202", TRUE ~ replicates),
         replicates = case_when(spygen_code == "SPY2401235" | spygen_code == "SPY2404175" | spygen_code == "SPY2402757" | spygen_code == "SPY2404190" ~ "SPY2401235/SPY2402757/SPY2404175/SPY2404190", TRUE ~ replicates))

#### Reprendre la colonne replicates ####
# Travail sur les replicats pour y incorporer le _ provenant des pools dans la colonne replicates
# et les classer en mettant en premier le pool dont code SPY est le plus petit 
# Classer les replicates par ordre croissant si deux pools mettre le premier plus petit code SPY en premier (avant le /)
raw_15 <- raw_14 %>%
  mutate(replicates = str_replace(replicates, "SPY232643/SPY232649", "SPY232643_SPY232649"),
         replicates = str_replace(replicates, "SPY232658/SPY232667", "SPY232658_SPY232667"),
         replicates = str_replace(replicates, "SPY232651/SPY232669", "SPY232651_SPY232669"),
         replicates = str_replace(replicates, "SPY232652/SPY232661", "SPY232652_SPY232666"),
         replicates = str_replace(replicates, "SPY232666/SPY232670", "SPY232661_SPY232670"),
         replicates = str_replace(replicates, "SPY232659/SPY232668", "SPY232659_SPY232668"),
         replicates = str_replace(replicates, "SPY2401148/SPY2401163", "SPY2401148_SPY2402767"),
         replicates = str_replace(replicates, "SPY2401168/SPY2402767", "SPY2401163_SPY2401168"),
         replicates = str_replace(replicates, "SPY2401179/SPY2401231", "SPY2401179_SPY2402764"),
         replicates = str_replace(replicates, "SPY2402764/SPY2402766", "SPY2401231_SPY2402766"),
         replicates = str_replace(replicates, "SPY2401185/SPY2401191", "SPY2401185_SPY2402759"),         
         replicates = str_replace(replicates, "SPY2401249/SPY2402759", "SPY2401191_SPY2401249"),
         replicates = str_replace(replicates, "SPY2401187/SPY2401246", "SPY2401187_SPY2401246"),
         replicates = str_replace(replicates, "SPY2401250/SPY2401257", "SPY2401250_SPY2401257"),
         replicates = str_replace(replicates, "SPY2401226/SPY2401239", "SPY2401226_SPY2401251"),
         replicates = str_replace(replicates, "SPY2401251/SPY2401255", "SPY2401239_SPY2401255"),
         replicates = str_replace(replicates, "SPY2401227/SPY2401256", "SPY2401227_SPY2401256"),
         replicates = str_replace(replicates, "SPY2402756/SPY2402771", "SPY2402756_SPY2402771"),
         replicates = str_replace(replicates, "SPY2401230/SPY2401236", "SPY2401230_SPY2401248"),
         replicates = str_replace(replicates, "SPY2401238/SPY2401248", "SPY2401236_SPY2401238"),
         replicates = str_replace(replicates, "SPY2401232/SPY2401237", "SPY2401232_SPY2401254"),
         replicates = str_replace(replicates, "SPY2401241/SPY2401254", "SPY2401237_SPY2401241"),
         replicates = str_replace(replicates, "SPY2401233/SPY2401234", "SPY2401233_SPY2401240"),
         replicates = str_replace(replicates, "SPY2401240/SPY2401242", "SPY2401234_SPY2401242"),
         replicates = str_replace(replicates, "SPY2401235/SPY2402757", "SPY2401235_SPY2404175"),
         replicates = str_replace(replicates, "SPY2404175/SPY2404190", "SPY2402757_SPY2404190"),
         replicates = str_replace(replicates, "SPY2401244/SPY2401245", "SPY2401244_SPY2401252"),
         replicates = str_replace(replicates, "SPY2401252/SPY2401253", "SPY2401245_SPY2401253"),
         replicates = str_replace(replicates, "SPY2404174/SPY2404178", "SPY2404174_SPY2404181"),
         replicates = str_replace(replicates, "SPY2404181/SPY2404202", "SPY2404178_SPY2404202")) %>%
  add_column(., "protection", .after = "habitat") %>%
  rename(protection = '"protection"') %>%
  mutate(protection = ifelse(protection == protection, NA, protection)) %>%
  relocate(protection, .after = habitat)%>%
  
  add_column(., "mpa_name", .after = "protection") %>%
  rename(mpa_name = '"mpa_name"') %>%
  mutate(mpa_name = ifelse(mpa_name == mpa_name, NA, mpa_name)) %>%
  relocate(mpa_name, .after = protection) %>%
  
  add_column(., "mpa_dist", .after = "protection") %>%
  rename(mpa_dist = '"mpa_dist"') %>%
  mutate(mpa_dist = ifelse(mpa_dist == mpa_dist, NA, mpa_dist)) %>%
  relocate(mpa_dist, .after = protection) %>%
  
  relocate(spygen_code,	pool,replicates,	date,	country,	region,	site,	subsite,subsite_andromede,	component,	habitat,	protection,
           mpa_name,	mpa_dist,	latitude_start_DD,	longitude_start_DD,	latitude_end_DD,	longitude_end_DD,	depth_seafloor,
           depth_sampling,	method,	time_start,	duration,	filter,	estimated_volume, flow_sensor_volume,	project,	lockdown,	BiodivMed2023,	Tele01,
           Pleo,	Mamm01,	Vert01,	X16s_Metazoa,	Bact02,	Euka02,	comments) %>%
  arrange(site, subsite, date) %>%
# Je convertis en numérique les colonnes nécessaires sinon bind_rows ci-dessous pas possible
  mutate(latitude_start_DD = as.numeric(latitude_start_DD)) %>%
  mutate(longitude_start_DD = as.numeric(longitude_start_DD)) %>%
  mutate(latitude_end_DD = as.numeric(latitude_end_DD)) %>%
  mutate(longitude_end_DD = as.numeric(longitude_end_DD)) %>%
  mutate(depth_seafloor = as.numeric(depth_seafloor)) %>%
  mutate(depth_sampling = as.numeric(depth_sampling))  %>%
  mutate(duration = as.numeric(duration)) %>%
  mutate(filter = as.numeric(filter)) %>%
  mutate(estimated_volume = as.numeric(estimated_volume)) %>%
  mutate(flow_sensor_volume = as.numeric(flow_sensor_volume)) %>%
  mutate(BiodivMed2023 = as.numeric(BiodivMed2023))

# Si dans la colonne replicates il n'y a qu'un seul "_" qui apparaît alors remplace par "no"
# Car certains pools (SPYA_SPYB) ont été remplis par défaut dans replicates, alors qu'ils n'ont pas de réplicat
#replicates = ifelse(str_count(replicates, "_") == 1, "no", replicates),
# NE SURTOUT PAS METTRE LA LIGNE CI-DESSUS : car pour les duplicats site qui ne sont pas poolés (et séparés par un / on laisse dans replicat)

# Je génère le fichier en .csv
write.csv(raw_15, "/Users/amandineavouac/Documents/Metadatas_Med_2018-2023/2024_datas/raw_15.csv")
# contient 332 obs et 37 var

# Je combine les métadonnées retravaillées de 2018-2023 et celles de 2024
# Lecture du fichier de métadonnées MED 2018-2023 retravaillé en 2024
metadatas_med_2018_2023_clean <- read.csv("/Users/amandineavouac/Documents/Metadatas_Med_2018-2023/2024_datas/raw_6.csv") %>%
  select(-X)
# contient 1520 obs et 37 var

raw_16 <- bind_rows(metadatas_med_2018_2023_clean, raw_15) %>%
  # j'ai bien 1852 obs (332+1520) et 37 var donc correct j'ai tous les éch de 2018-2024 inclus
  arrange(site, subsite, date) 

# Je génère le fichier 
write.csv(raw_16, "/Users/amandineavouac/Documents/Metadatas_Med_2018-2023/2024_datas/Med_metadonnees_ADNe_2018_2024_V1.csv")
write.csv(raw_16, "/Users/amandineavouac/Documents/Metadatas_Med_2024/datas/raw_16.csv")

#### Quelques coquilles ont été relevées donc j'en prends cpte : replicats Calvi 2020, time_start éch Calvi 2023, Majorque/Minorque banyuls ####
# je repars de la dernière version des métadonnées V3 ci-dessus "/Users/amandineavouac/Documents/Metadatas_Med_2024/datas/Med_metadonnees_ADNe_2018_2024_V3.csv" (dont les component checkés)
raw_18 <- read.csv2("datas/Med_metadonnees_ADNe_2018_2024_V2.csv", sep = ",") %>%
# email Laure : pour faire suite à notre discussion d'hier voici les réplicats pour l'échantillonnage Calvi 2020
# SPY201221/SPY201227/SPY201238/SPY201264 (réserve)
# SPY201217/SPY201233/SPY201234/SPY201284 (5 km)
# SPY201214/SPY201219/SPY201235/SPY201290 (10km)
  select(-X.1) %>%
  mutate(replicates = ifelse(spygen_code == "SPY201221", "SPY201221/SPY201227/SPY201238/SPY201264", replicates),
         replicates = ifelse(spygen_code == "SPY201227", "SPY201221/SPY201227/SPY201238/SPY201264", replicates),
         replicates = ifelse(spygen_code == "SPY201238", "SPY201221/SPY201227/SPY201238/SPY201264", replicates),
         replicates = ifelse(spygen_code == "SPY201264", "SPY201221/SPY201227/SPY201238/SPY201264", replicates),
         replicates = ifelse(spygen_code == "SPY201217", "SPY201217/SPY201233/SPY201234/SPY201284", replicates),
         replicates = ifelse(spygen_code == "SPY201233", "SPY201217/SPY201233/SPY201234/SPY201284", replicates),
         replicates = ifelse(spygen_code == "SPY201234", "SPY201217/SPY201233/SPY201234/SPY201284", replicates),
         replicates = ifelse(spygen_code == "SPY201284", "SPY201217/SPY201233/SPY201234/SPY201284", replicates),
         replicates = ifelse(spygen_code == "SPY201214", "SPY201214/SPY201219/SPY201235/SPY201290", replicates),
         replicates = ifelse(spygen_code == "SPY201219", "SPY201214/SPY201219/SPY201235/SPY201290", replicates),
         replicates = ifelse(spygen_code == "SPY201235", "SPY201214/SPY201219/SPY201235/SPY201290", replicates),
         replicates = ifelse(spygen_code == "SPY201290", "SPY201214/SPY201219/SPY201235/SPY201290", replicates)) %>%
  # time_start des éch "SPY233816" "SPY233817" "SPY233818" "SPY233819" "SPY233824" "SPY233825" pas au bon format
  mutate(time_start = ifelse(time_start == "10h", "10:00", time_start),
         time_start = ifelse(time_start == "11h", "11:00", time_start),
         time_start = ifelse(time_start == "14h", "14:00", time_start),
         time_start = ifelse(time_start == "15h", "15:00", time_start)) %>%
  
  # site ne concorde pas avec region Minorque/banyuls ; Majorque/banyuls
  mutate(site = ifelse(region == "Minorque", NA, site),
         site = ifelse(region == "Majorque", NA, site)) %>%
  
#### Question MedPortADNe Stéphanie/Bastien quel code SPY est réplicat de quel autre ####
mutate(replicates = case_when(spygen_code == "SPY202505" | spygen_code == "SPY202519" ~ "SPY202505/SPY202519", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY202520" ~ "no", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY213512" | spygen_code == "SPY213516" ~ "SPY213512/SPY213516", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY213515" ~ "no", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY202505" | spygen_code == "SPY202519" ~ "SPY202505/SPY202519", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY202520" ~ "no", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY213521" | spygen_code == "SPY213522" ~ "SPY213521/SPY213522", TRUE ~ replicates),
       replicates = case_when(spygen_code == "SPY213518" ~ "no", TRUE ~ replicates)) %>%
  select(-X) %>%
  arrange(region, site) %>%
 
#### GENERER NOUVELLE VERSION FICHIER ####
write.csv(raw_18, "datas/Med_metadonnees_ADNe_2018_2024_V4.csv")

#### INFOS ####
# je peux regrouper les replicats pour les site == "marseille" d'après les fiches terrain (faut que je les aies)

# pour calanques je peux avoir : SPYA_SPYB/SPYC_SPYD ; SPYA_SPYB/SPYC pas de D car pompe n'a pas marchée 

# voir comment faire pour gérer les réplicats 4 filtres : réserve benefits / eREF (en 2020 : confinement, / CERCA (Calvi trois réplicats car 1 pompe avait lâché) 

# Question MedPortADNe Stéphanie/Bastien
# Bastien : "Je n'étais pas présent lors de cet échantillonnage, mais il me semble que dans ces cas là, deux réplicats ont été fait en bateau (transect), et pour le biohut il s'agit d'un point à part qui n'a pas été fait en transect (point fixe).
# J'aurais donc tendance à dire que le biohut n'est pas un réplicat des deux autres comme il ne s'agit pas d'un transect, et que donc la technique est différente".

  
  
