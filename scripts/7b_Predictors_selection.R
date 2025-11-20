#------------- Description ---------------------
# Purpose: 
# This script aims to compute select predictors. 

# Author: Marieke Schultz

# Date script created: 16/11/2025

#------------- Setting up ------------------
# Remove existing objects
rm(list = ls())

# Set current working directory
setwd("/media/marieke/Shared/Chap-1/Model/Scripts/Chap_1_2018-2024")

# Libraries
library(corrplot)
library(dplyr)
library(exactextractr)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(terra)
library(stringr)
library(pMEM)



# Load functions
source("./utils/Fct_Data-Prep.R")

#------------- Load data ------------------
# predictors_tr_v.1.1 ----
pred <- st_read("./data/processed_data/predictors/predictors_tr_v1.1.gpkg")
# predictors_tr_v.1.2 ----
pred <- st_read("./data/processed_data/predictors/predictors_tr_v1.2.gpkg")

# mtdt_7_sel_v1.0 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.0.gpkg")
# mtdt_7_sel_v1.1 ----
mtdt <- st_read("./data/processed_data/Mtdt/mtdt_7_sel_v1.1.gpkg")

#---------------- REMOVE ROWS BASED ON SITE SELECTION -----

# Remove pred rows not in mtdt ----
pred <- pred %>%
  dplyr::filter(replicates %in% mtdt$replicates)

rm(mtdt)

#---------------- PRELIMINAR PREDICTORS SELECTION / CLEANING ----------------

# Drop geometry ----
sel <- pred %>%
  st_drop_geometry() 


# Remove range ----
# Remove predictors containing "range"
range_cols <- grepl("_range", colnames(sel))

colnames(sel[which(range_cols)])


sel <- sel %>%
  dplyr::select(which(!range_cols))

rm(range_cols)
# Remove non weighted mean distances ----
sel <- sel %>%
  dplyr::select(-c(shore_dist_m_mean_log, mpa_dist_m_mean_log, port_dist_m_mean))
#--------------- PREDICTORS SELECTION v.1.0 ---------------- 

# Remove min and max cols ----
min_cols <- grepl("_min", colnames(sel))
max_cols <- grepl("_max", colnames(sel))
colnames(sel[which(min_cols)])
colnames(sel[which(max_cols)])

drop_cols <- min_cols | max_cols

sel <- sel[, !drop_cols]

rm(min_cols, max_cols, drop_cols)

# Remove day week and year cols ----
# remove cols containing "7j", "24h", "1y", "_year", "5years", "week", "day"
day_cols <- grepl("_7j", colnames(sel)) | grepl("_24h", colnames(sel)) | grepl("_1y", colnames(sel)) | grepl("_year", colnames(sel)) | grepl("5years", colnames(sel)) | grepl("week", colnames(sel)) | grepl("day", colnames(sel))

colnames(sel[which(day_cols)])

sel <- sel[, !day_cols]

rm(day_cols)

# Remove canyon associated cols (except distance) ----

canyon_cols <- grepl("canyon", colnames(sel)) & !grepl("canyon_dist_m_weight_log", colnames(sel))
colnames(sel[which(canyon_cols)])

sel <- sel[, !canyon_cols]


rm(canyon_cols)
# Remove not grouped habitat ----
hab_cols <- c("nb_habitat_per_km2", "main_habitat")

sel <- sel %>%
  dplyr::select(-all_of(hab_cols))

rm(hab_cols)

# Remove single habitat surface -----
# cols containing "clr"
clr_cols <- grepl("_clr", colnames(sel))

colnames(sel[which(clr_cols)])

sel <- sel %>%
  dplyr::select(-which(clr_cols))

rm(clr_cols)
# Remove sampling effort cols ----
samp_cols <- c("dist_seabed_depthsampling_log", "area_km2")

sel <- sel %>%
  dplyr::select(-all_of(samp_cols))

rm(samp_cols)

# Remove mpa_protection and mpa_fully ----
sel <- sel %>%
  dplyr::select(-mpa_protection) %>%
  dplyr::select(-mpa_fully)
# # Correlation matrix -----
# # Compute Spearman correlation matrix
# # Select numerical cols
# df <- sel %>%
#   dplyr::select(where(is.numeric))
# 
# # Compute correlation matrix
# cor_mat <- cor(df, method = "spearman", use = "pairwise.complete.obs")
# 
# # Open PNG device *before* plotting
# png("./figures/Predictors/Selection/corr_matrix_sel_v1.0_pre-selection.png",
#     width = 2000, height = 2000, res = 300)
# 
# # Plot correlation matrix
# corrplot(cor_mat,
#          method = "color",
#          type = "upper",
#          tl.col = "black",
#          addCoef.col = "black", 
#          tl.cex = 0.5,
#          cl.cex = 0.5,
#          number.cex = 0.5)
# 
# # Close device (this actually writes the file)
# dev.off()
# 
# # Results 
# # Very high correlation > 0.9 between : 
# # cop_analysed_sst_month_mean and temp_mean_1m --> we remove cop_analysed_sst_month_mean
# # roughness & tri & slope
# 
# # High correlation > 0.7 between :
# # canyon_dist and roughness, slope and tri
# 
# 
# # VIF ----
# # Explanation : VIF can be used to detect collinearity (Strong correlation between two or more predictor variables). Collinearity causes instability in parameter estimation in regression-type models. The VIF is based on the square of the multiple correlation coefficient resulting from regressing a predictor variable against all other predictor variables. If a variable has a strong linear relationship with at least one other variables, the correlation coefficient would be close to 1, and VIF for that variable would be large. A VIF greater than 10 is a signal that the model has a collinearity problem. (source = usdm documentation)
# 
# # vifcor : first finds a pair of variables which has the maximum linear correlation (greater than the threshold; th), and exclude the one with a greater VIF. The procedure is repeated untill no pair of variables with a high corrrelation coefficient (grater than the threshold) remains.
# 
# # vifstep calculates VIF for all variables, excludes the one with the highest VIF (if it is greater than the threshold), repeat the procedure untill no variables with a VIF greater than th remains.
# 
# 
# # Select numerical cols
# df <- sel %>%
#   dplyr::select(where(is.numeric))
# 
# 
# # Print all VIF values
# vif <- usdm::vif(df)
# vif[order(-vif$VIF), ]
# 
# # Perform stepwise VIF selection
# p_vifstep <- usdm::vifstep(df, th = 10, keep = NULL)
# 
# # Perform VIF correlation selection
# p_vifcor <- usdm::vifcor(df, th = 0.7, keep = NULL, method = 'spearman') # remaining predictors with corr < th # default th is 0.9 # # make again with th = 0.7 method = spearman because we have non-normal data.
# 
# # Filter vif selected predictors
# df_vifcor <- df[, p_vifcor@results$Variable]
# df_vifstep <- df[, p_vifstep@results$Variable]
# 
# # Check removed predictors
# setdiff(names(df), names(df_vifcor))
# # [1] "temp_mean_1m"       "roughness_mean_log" "slope_mean_log"     "tri_mean_log"      
# setdiff(names(df), names(df_vifstep))
# # [1] "slope_mean_log" "tri_mean_log"  
# 
# 
# rm(cor_mat, vif, p_vifcor, p_vifstep, df, df_vifcor, df_vifstep)

# Selection based on VIF and correlation results ----
# Remove roughness slope and tri -----
# slope and tri avec removed with vistep
# roughness is correlated to canyon distance and is removed by vifcor

sel <- sel %>%
  dplyr::select(-c(roughness_mean_log, slope_mean_log, tri_mean_log))


# Remove sst ----
# sst is highly correlated to temp_mean_1m and temp_mean is removed by vifcor. We prefer to keep mars3d temp that takes into account depth. 

sel <- sel %>%
  dplyr::select(-cop_analysed_sst_month_mean)



# # Correlation matrix after selection -----
# 
# # Compute Spearman correlation matrix
# df <- sel %>%
#   dplyr::select(where(is.numeric))
# 
# cor_mat <- cor(df, method = "spearman", use = "pairwise.complete.obs")
# 
# # Open PNG device *before* plotting
# png("./figures/Predictors/Selection/corr_matrix_sel_v1.0_post-selection.png",
#     width = 2000, height = 2000, res = 300)
# 
# # Plot correlation matrix
# corrplot(cor_mat,
#          method = "color",
#          type = "upper",
#          tl.col = "black",
#          addCoef.col = "black", 
#          tl.cex = 0.5,
#          cl.cex = 0.5,
#          number.cex = 0.5)
# 
# # Close device (this actually writes the file)
# dev.off()
# 
# # Clean 
# rm(cor_mat, df)



# Add back geometry from pred ----
sel <- sel %>%
  left_join(pred %>% dplyr::select(replicates, geom),
            by = "replicates") %>%
  st_as_sf()



#--------------- EXPORTS ----------------
# # Export predictors_sel_v1.0.gpkg ----
# # based on mtdt_7_sel_v1.0.gpkg (755 obs) and predictors_tr_v1.1.gpkg (142 var + replicates)
# # selected predictors (16 + + x, y, replicates, geom): 
# colnames(sel)
# # [1] "grouped_main_habitat"       "x"                          "y"                          "replicates"                 "port_dist_m_weight"        
# # [6] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# # [11] "sal_mean_1m"                "northness"                  "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"      "shore_dist_m_weight_log"   
# # [16] "gravity_mean_log"           "tpi_mean_log"               "cop_chl_month_mean_log"     "eastness_log"               "geom"  
# 
# 
# 
# st_write(sel, "./data/processed_data/predictors/predictors_sel_v1.0.gpkg", delete_dsn = TRUE)
# 
# 
# # Export predictors_sel_v1.2.gpkg ----
# # based on mtdt_7_sel_v1.0.gpkg (755 obs) and predictors_tr_v1.2.gpkg (142 var + replicates --> x and y not transformed + negative log)
# # selected predictors (16 + + x, y, replicates, geom): 
# colnames(sel)
# # [1] "grouped_main_habitat"       "x"                          "y"                          "replicates"                 "port_dist_m_weight"        
# # [6] "grouped_nb_habitat_per_km2" "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"                "temp_mean_1m"              
# # [11] "sal_mean_1m"                "northness"                  "canyon_dist_m_weight_log"   "mpa_dist_m_weight_log"      "shore_dist_m_weight_log"   
# # [16] "gravity_mean_log"           "tpi_mean_log"               "cop_chl_month_mean_log"     "eastness_log"               "geom"    
# 
# 
# st_write(sel, "./data/processed_data/predictors/predictors_sel_v1.1.gpkg", delete_dsn = TRUE)
# 

# Export predictors_sel_v1.3.gpkg ----
# based on mtdt_7_sel_v1.1.gpkg (637 obs) and predictors_tr_v1.2.gpkg (142 var + replicates) 
# selected predictors (16 + + x, y, replicates, geom):
# [1] "grouped_main_habitat"       "x"                          "y"                         
# [4] "replicates"                 "northness"                  "eastness"                  
# [7] "tpi_mean_log"               "port_dist_m_weight"         "grouped_nb_habitat_per_km2"
# [10] "bathy_mean"                 "wind_mean_1m"               "vel_mean_1m"               
# [13] "temp_mean_1m"               "sal_mean_1m"                "canyon_dist_m_weight_log"  
# [16] "mpa_dist_m_weight_log"      "shore_dist_m_weight_log"    "gravity_mean_log"          
# [19] "cop_chl_month_mean_log"    "geom" 
st_write(sel, "./data/processed_data/predictors/predictors_sel_v1.3.gpkg", delete_dsn = TRUE)
# Export predictors_raw_sel_v1.3.gpkg ----
# based on predictors_raw_v2.2 and predictors_sel_v1.3 
# selects the same predictors and rows but from RAW predictors (for data explo)

pred_raw <- st_read("./data/processed_data/predictors/predictors_raw_v2.2.gpkg")
predictors_sel_v1.3 <- st_read("./data/processed_data/predictors/predictors_sel_v1.3.gpkg")

sel_raw <- pred_raw %>%
  dplyr::filter(replicates %in% predictors_sel_v1.3$replicates) %>% 
  dplyr::select(c("replicates",
                  "grouped_main_habitat",
                  "x",
                  "y",
                  "aspect_mean",
                  "tpi_mean",
                  "port_dist_m_weight",
                  "grouped_nb_habitat_per_km2",
                  "bathy_mean",
                  "wind_mean_1m",
                  "vel_mean_1m",
                  "temp_mean_1m",
                  "sal_mean_1m",
                  "canyon_dist_m_weight",
                  "mpa_dist_m_weight",
                  "shore_dist_m_weight",
                  "gravity_mean",
                  "cop_chl_month_mean",
                  "geom"))

st_write(sel_raw, "./data/processed_data/predictors/predictors_raw_sel_v1.3.gpkg", delete_dsn = TRUE)

#----------------- AUTRES -----------------
#---- Visualise selected pred ----
method = 'vifcor' # 'vifcor'

if (method == 'vifstep') {
  pred_sel <- df_vifstep
  color <- "hotpink"
} else if (method == 'vifcor') {
  pred_sel <- df_vifcor
  color <- "skyblue"
}


# Correlation matrix
cm <- generate_correlation_matrix(data = df, exclude_column = NULL) # Spearman correlation

# Distance matrix
dm <- as.dist(1 - cm) # Dissimilarity measure (pearson = 1 -> 0 = identical, pearson = -1 -> -2 = maw dist)

# Dimension reduction for mapping/ visualise the data
mds_df <- as.data.frame(cmdscale(dm, k = 2))  # 2D scaling = PCA
mds_df <- cbind(Variable = rownames(mds_df), mds_df)  # Add Variable column
colnames(mds_df) <- c("Variable", "Dim1", "Dim2")  # Rename columns
rownames(mds_df) <- NULL  # Remove row names

# Add a column indicating whether the variable is in pred_sel
mds_df$Selected <- ifelse(mds_df$Variable %in% colnames(pred_sel), "Selected", "Not Selected")

# Create the plot
ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Variable, color = Selected)) +
  geom_point(size = 2.5) +
  geom_text_repel(size = 2.5, max.overlaps = Inf) +  # Ensure all labels appear
  scale_color_manual(values = c("Selected" = color, "Not Selected" = "grey30")) +  # Highlight selected variables
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = "white"),  # White panel background
    plot.background = element_rect(fill = "white", color = "white")   # White plot background
  ) +
  labs(
    title = paste0("MDS with Predictors selected ", method),
    x = "Dimension 1",
    y = "Dimension 2",
    color = paste0("Variable Selection\n(", sum(mds_df$Selected == "Selected"), " Pred selected)")
  )






#---- Multivariate analyses [TO DO] ----


""" From Celia """
"
# see if some covariates are too coorelated and choose 6 to 7 predictors (compare with the mean number of occurences)
# dataframe with all occurences and all covariates + Row.names

id_spygen_values <- c(occ_fitting$id_spygen, occ_validation$id_spygen)
filtered_occ <- occ[occ$id_spygen %in% id_spygen_values, ]
filtered_cov <- data.frame(cov_med$id_spygen)

common_id_spygen <- intersect(filtered_occ$id_spygen, filtered_cov$cov_med.id_spygen)

# filtered_occ[, -1] <- as.data.frame(lapply(filtered_occ[, -1], as.numeric))

filtered_occ_2 <- occ[occ$id_spygen %in% common_id_spygen, ]

combined_data <- dplyr::full_join(filtered_occ_2, cov_med, by = "id_spygen")

## RDA
# utiliser data_MED_final
RDA2=vegan::capscale(combined_data[,-c(1,238:250)] ~ nb_substrats + protection + principal_substrat + bathymetry + dist_fully_protected_MPA + chlorophyll  + temperature + salinity, combined_data, dist="jaccard", na.action = na.omit, add =TRUE)
summary(RDA2)
plot(RDA2)

## AFC 
afc <- FactoMineR::CA(occ[,-1], graph = FALSE, axes = c(1,2))
factoextra::fviz_ca_biplot(afc, axes = c(1, 2), col.row = "blue", col.col = "black", repel= FALSE, ellipse.type="confidence", addELLipse=TRUE)


# ACP
# ici problème avec les NA vérifier qu'on a bien aucun NA dans le fichier final et sauter cette étape 
# Create a dataset with all quantitative covariates
data_quantitative <- cov_med[c("gravity", "dist_fully_protected_MPA", "dist_shore", "fishing_pressure", "chlorophyll", "turbidity","salinity", "temperature", "nb_substrats", "bathymetry")]

pca = FactoMineR::PCA(data_quantitative[,-1], scale.unit = TRUE, ncp = 5, graph = FALSE)
# Get eigenvalues
eig.val <- factoextra::get_eigenvalue(pca)
eig.val
# Draw screeplot
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))

# Variable contributions to axes
var <- factoextra::get_pca_var(pca)
# Corrplot of the cos2 of variables
corrplot::corrplot(var$cos2, is.corr=FALSE)
# Total cos2 of variables on Dim.1 to Dim.3
factoextra::fviz_cos2(pca, choice = "var", axes = 1:3)
# Contributions of variables to PC1 : 
factoextra::fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2 :
factoextra::fviz_contrib(pca, choice = "var", axes = 2, top = 10)
# Variable contributions to axes
factoextra::fviz_pca_var(pca, col.var = "contrib",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
""

