mtdt_3 <- st_read("./data/processed_data/Mtdt/mtdt_3.gpkg")
mtdt_i <- mtdt_3 %>% filter(project == "IPOCOM")

pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v1.1.gpkg")
pred_raw_i <- pred_raw %>%
  filter(replicates %in% unique(mtdt_i$replicates))

length(unique(pred_raw_i$replicates))

# export 
st_write(mtdt_i, "./scripts/Old_codes/mtdt_ipocom.gpkg", delete_dsn = TRUE)
st_write(pred_raw_i, "./scripts/Old_codes/predictors_ipocom.gpkg", delete_dsn = TRUE)
