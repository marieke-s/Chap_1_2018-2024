#------------- Description ----
# The aim of this script is to extract predictor values for the predictions cells. 





#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(compositions)
library(dplyr)
library(exactextractr)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(terra)
library(stringr)
library(pMEM)
library(purrr)



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- Load Predictions ------------------
# T1.33 ----
bottom <- readRDS("./output/predictions/rf/T1.33/T1.33_R_predictions_bottom_2023-av-sept.rds")
surf <- readRDS("./output/predictions/rf/T1.33/T1.33_R_predictions_surface_2023-av-sept.rds")


# Check
dt <- bottom$`2023-05-01`
dt <- surf$`2023-05-01`
summary(dt$prediction_R)

# Loop through all months of surf + bottom and print summary of prediction_R
for (date in names(bottom)) {
  cat("--- Bottom - Date:", date, "\n")
  print(summary(bottom[[date]]$prediction_R))
  
  cat("--- Surface - Date:", date, "\n")
  print(summary(surf[[date]]$prediction_R))
}

# Violin plots of predictions
# 1) Gather bottom + surface predictions into one data.frame
bottom_df <- map_dfr(names(bottom), ~ {
  bottom[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),   # e.g. "mai 2023"
      layer = "Fond",
      prediction_R
    )
})

surf_df <- map_dfr(names(surf), ~ {
  surf[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),
      layer = "Surface",
      prediction_R
    )
})

tot_pred <- bind_rows(bottom_df, surf_df) %>%
  mutate(
    month = factor(month, levels = unique(month))  # keep chronologic order
  )

# 2) Violin + boxplot + jitter, Surface vs Fond for each month
# Manual colors
col_bottom  <- "#08306B"   # dark blue
col_surface <- "#6BAED6"   # light blue
col_points  <- "grey80"

pos_dodge <- position_dodge(width = 0.8)

p_pred <- ggplot(tot_pred, aes(x = month, y = prediction_R, fill = layer)) +
  
  # Violin
  geom_violin(position = pos_dodge, color = NA, width = 0.9, alpha = 0.8) +
  
  # Boxplot
  geom_boxplot(position = pos_dodge, width = 0.25, outlier.shape = NA, alpha = 0.7) +
  
  # Points (jitter)
  geom_jitter(
    aes(color = layer),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.1, size = 0.8, show.legend = FALSE
  ) +
  
  # Median point
  stat_summary(fun = median, geom = "point",
               position = pos_dodge, size = 2.3, color = "black",
               show.legend = FALSE) +
  
  # Custom scale colors
  scale_fill_manual(values = c("Fond" = col_bottom, "Surface" = col_surface)) +
  scale_color_manual(values = c("Fond" = col_points, "Surface" = col_points)) +
  
  labs(
    title = "Prédictions de R par mois et couche",
    x = "Mois",
    y = "R prédit",
    fill = "Couche"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_pred
ggsave(
  p_pred,
  filename = "./figures/T1.33_R_ViolinPlot_Predictions_by_Month_and_Layer.png",
  width = 20, height = 6, dpi = 300
)




# T1.37 ----
bottom <- readRDS("./output/predictions/rf/T1.37/T1.37_RedList_predictions_bottom_2023-av-sept.rds")
surf <- readRDS("./output/predictions/rf/T1.37/T1.37_RedList_predictions_surface_2023-av-sept.rds")


# Check
dt <- bottom$`2023-05-01`
dt <- surf$`2023-05-01`
summary(dt$prediction_R)

# Loop through all months of surf + bottom and print summary of prediction_R
for (date in names(bottom)) {
  cat("--- Bottom - Date:", date, "\n")
  print(summary(bottom[[date]]$prediction_RedList))
  
  cat("--- Surface - Date:", date, "\n")
  print(summary(surf[[date]]$prediction_RedList))
}

# Violin plots of predictions
# 1) Gather bottom + surface predictions into one data.frame
bottom_df <- map_dfr(names(bottom), ~ {
  bottom[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),   # e.g. "mai 2023"
      layer = "Fond",
      prediction_RedList
    )
})

surf_df <- map_dfr(names(surf), ~ {
  surf[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),
      layer = "Surface",
      prediction_RedList
    )
})

tot_pred <- bind_rows(bottom_df, surf_df) %>%
  mutate(
    month = factor(month, levels = unique(month))  # keep chronologic order
  )

# 2) Violin + boxplot + jitter, Surface vs Fond for each month
# Manual colors
col_bottom  <- "#08306B"   # dark blue
col_surface <- "#6BAED6"   # light blue
col_points  <- "grey80"

pos_dodge <- position_dodge(width = 0.8)

p_pred <- ggplot(tot_pred, aes(x = month, y = prediction_RedList, fill = layer)) +
  
  # Violin
  geom_violin(position = pos_dodge, color = NA, width = 0.9, alpha = 0.8) +
  
  # Boxplot
  geom_boxplot(position = pos_dodge, width = 0.25, outlier.shape = NA, alpha = 0.7) +
  
  # Points (jitter)
  geom_jitter(
    aes(color = layer),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.1, size = 0.8, show.legend = FALSE
  ) +
  
  # Median point
  stat_summary(fun = median, geom = "point",
               position = pos_dodge, size = 2.3, color = "black",
               show.legend = FALSE) +
  
  # Custom scale colors
  scale_fill_manual(values = c("Fond" = col_bottom, "Surface" = col_surface)) +
  scale_color_manual(values = c("Fond" = col_points, "Surface" = col_points)) +
  
  labs(
    title = "Prédictions de RedList par mois et couche",
    x = "Mois",
    y = "R prédit",
    fill = "Couche"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_pred

ggsave(
  p_pred,
  filename = "./figures/T1.37_RedList_ViolinPlot_Predictions_by_Month_and_Layer.png",
  width = 20, height = 6, dpi = 300
)










# crypto ------
bottom <- readRDS("./output/predictions/rf/T1.37/T1.37_Crypto_predictions_bottom_2023-av-sept.rds")
surf <- readRDS("./output/predictions/rf/T1.37/T1.37_Crypto_predictions_surface_2023-av-sept.rds")


# Check
dt <- bottom$`2023-05-01`
dt <- surf$`2023-05-01`
summary(dt$prediction_R)

# Loop through all months of surf + bottom and print summary of prediction_R
for (date in names(bottom)) {
  cat("--- Bottom - Date:", date, "\n")
  print(summary(bottom[[date]]$prediction_Crypto))
  
  cat("--- Surface - Date:", date, "\n")
  print(summary(surf[[date]]$prediction_Crypto))
}

# Violin plots of predictions
# 1) Gather bottom + surface predictions into one data.frame
bottom_df <- map_dfr(names(bottom), ~ {
  bottom[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),   # e.g. "mai 2023"
      layer = "Fond",
      prediction_Crypto
    )
})

surf_df <- map_dfr(names(surf), ~ {
  surf[[.x]] %>%
    st_drop_geometry() %>%
    transmute(
      date  = as.Date(.x),
      month = format(date, "%b %Y"),
      layer = "Surface",
      prediction_Crypto
    )
})

tot_pred <- bind_rows(bottom_df, surf_df) %>%
  mutate(
    month = factor(month, levels = unique(month))  # keep chronologic order
  )

# 2) Violin + boxplot + jitter, Surface vs Fond for each month
# Manual colors
col_bottom  <- "#08306B"   # dark blue
col_surface <- "#6BAED6"   # light blue
col_points  <- "grey80"

pos_dodge <- position_dodge(width = 0.8)

p_pred <- ggplot(tot_pred, aes(x = month, y = prediction_Crypto, fill = layer)) +
  
  # Violin
  geom_violin(position = pos_dodge, color = NA, width = 0.9, alpha = 0.8) +
  
  # Boxplot
  geom_boxplot(position = pos_dodge, width = 0.25, outlier.shape = NA, alpha = 0.7) +
  
  # Points (jitter)
  geom_jitter(
    aes(color = layer),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.1, size = 0.8, show.legend = FALSE
  ) +
  
  # Median point
  stat_summary(fun = median, geom = "point",
               position = pos_dodge, size = 2.3, color = "black",
               show.legend = FALSE) +
  
  # Custom scale colors
  scale_fill_manual(values = c("Fond" = col_bottom, "Surface" = col_surface)) +
  scale_color_manual(values = c("Fond" = col_points, "Surface" = col_points)) +
  
  labs(
    title = "Prédictions de Crypto par mois et couche",
    x = "Mois",
    y = "R prédit",
    fill = "Couche"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_pred

ggsave(
  p_pred,
  filename = "./figures/T1.37_Crypto_ViolinPlot_Predictions_by_Month_and_Layer.png",
  width = 20, height = 6, dpi = 300
)



