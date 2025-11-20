# CHECKS pcr values [ ISSUE !! CF MAIL LAURE VELEZ 28/10/2025] -----

# ---- Load data ----
pcr <- read.csv("./data/processed_data/eDNA/data_MED_teleo_nb_rep_1824_V1.csv", sep = ";")

# Identify species columns
colnames(pcr)
species_cols <- setdiff(colnames(pcr), c("sample"))

# mtdt_3data
mtdt_3 <- sf::st_read("./data/processed_data/Mtdt/mtdt_3.gpkg")

mtdt_3 <- st_drop_geometry(mtdt_3)

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

# Remove "Ange2mer" project
mtdt_3 <- mtdt_3 %>% dplyr::filter(project != "Ange2mer")


# ---- Check max PCR values ----

# compute max per columns of df
sort(sapply(pcr[ , species_cols], max, na.rm = TRUE))
max(sapply(pcr[ , species_cols], max, na.rm = TRUE))

# PCR values above 12 ? 
t <- pcr %>% filter(Gobius_xanthocephalus == '48' ) %>% print() # 48 PCR replicates 
t$spygen_code # "SPY210641" --> maximum should be 12

rm(t)

# ---- CHECK Laure  ----





# Check whether all PCR values in a sample are multiples of (max / 12)
pcr_check <- pcr %>%
  rowwise() %>%
  mutate(
    min_val_above_0 = min(c_across(all_of(species_cols))[c_across(all_of(species_cols)) != 0], na.rm = TRUE),
    max_val = max(c_across(all_of(species_cols)), na.rm = TRUE),
    diviseur = max_val / 12,
    all_multiple = all((c_across(all_of(species_cols)) %% diviseur == 0) | is.na(c_across(all_of(species_cols))))
  ) %>%
  ungroup() %>%
  dplyr::select(sample, min_val_above_0, max_val, diviseur, all_multiple)


rm(pcr_check)




# ---- Plots ----



# make an hist of pcr across all columns for values > 0
# Extract and filter values (>0, non-NA)
pcr_values <- as.numeric(unlist(pcr[, species_cols]))
pcr_values <- pcr_values[pcr_values > 0 & !is.na(pcr_values)]



hist(pcr_values,
     breaks = 50,
     main = "Histogram of PCR values (>0)",
     xlab = "PCR values",
     col = "lightblue")

boxplot(pcr_values,
        main = "Boxplot of PCR values (>0)",
        ylab = "PCR values",
        col = "lightblue")

summary(pcr_values)


rm(pcr_values)

# ---- Select samples with PCR values > 12 ----

# count nb of as.numeric(unlist(pcr[ , species_cols])) > 12
sum(as.numeric(unlist(pcr[ , species_cols])) > 12, na.rm = TRUE) # 626 values above 12 --> WHY ? 

# Select rows containing species_cols values > 12
pcr_above_12 <- pcr %>%
  filter(if_any(all_of(species_cols), ~ . > 12))

# Rename sample to spygen_code
colnames(pcr_above_12)[colnames(pcr_above_12) == "sample"] <- "spygen_code"

# in mtdt_3 select spygen_code in pcr_above_12
mtdt_3_pcr_above_12 <- mtdt_3 %>%
  filter(spygen_code %in% pcr_above_12$spygen_code)

# Export mtdt_3_pcr_above_12 as csv
# Drop geometry
mtdt_3_pcr_above_12 <- st_drop_geometry(mtdt_3_pcr_above_12)
write.csv(mtdt_3_pcr_above_12, "./data/processed_data/eDNA/mtdt_3_pcr_above_12.csv", row.names = FALSE)

t <- pcr_above_12 %>% filter(spygen_code == "SPY2402598") %>% data.table::transpose() 

rm(t)



# ---- Check R ----

# Compute richness for filter with pcr > 12. Count nb of species with pcr > 0 per row
species_cols <- setdiff(colnames(pcr_above_12), c("spygen_code", "replicates", "geom", "pcr_replicates"))
pcr_above_12$R <- rowSums(pcr_above_12[ , species_cols] >0, na.rm = TRUE)

# Compute richness for all filters with pcr > 12
pcr$R <- rowSums(pcr[ , species_cols] >0, na.rm = TRUE)

summary(pcr_above_12$R)
summary(pcr$R)

# Statistically compare the richness of filters with pcr > 12 and all filters
t.test(pcr_above_12$R, pcr$R) # Significant difference (p-value < 0.05) --> filters with pcr > 12 have a significantly higher richness than all filters.

# Remove R columns
pcr_above_12 <- pcr_above_12 %>% dplyr::select(-R)
pcr <- pcr %>% dplyr::select(-R)

rm(t, pcr_above_12, mtdt_3_pcr_above_12, pcr_values)




