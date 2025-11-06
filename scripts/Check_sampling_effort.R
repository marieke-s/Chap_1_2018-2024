#------------- Check data for homogenous sampling effort ##############

buff <- as.data.frame(buff)
buff %>%
  group_by(PCR_replicates) %>%
  summarise(count = n()) %>%
  print()


buff %>%
  group_by(estimated_volume_total) %>%
  summarise(count = n()) %>%
  print()

buff %>%
  filter(estimated_volume_total > 55 & estimated_volume_total < 67 & PCR_replicates == 24) %>%
  dim()


summary(buff$area_km2)
hist(buff$area_km2, breaks = 10)
sd(buff$area_km2, na.rm = TRUE)


buff$area_km2 <- as.numeric(buff$area_km2)
buff %>%
  filter(estimated_volume_total > 55 & estimated_volume_total < 67 & PCR_replicates == 24 & area_km2 < 2) %>%
  dim()






