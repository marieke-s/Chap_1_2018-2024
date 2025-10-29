#---------------------------------- Function : buffer_transect --------------------
buffer_transect <- function(df,
                            start_lon_col,
                            start_lat_col,
                            end_lon_col = NULL,
                            end_lat_col = NULL,
                            buffer_dist = 500) # buffer distance is meters
  {
  
  # Build WKT geometry
  df <- df |>
    mutate(
      wkt_geometry = if (is.null(end_lon_col) || is.null(end_lat_col)) {
        paste0("POINT(", .data[[start_lon_col]], " ", .data[[start_lat_col]], ")")
      } else {
        ifelse(
          is.na(.data[[end_lon_col]]) | is.na(.data[[end_lat_col]]),
          paste0("POINT(", .data[[start_lon_col]], " ", .data[[start_lat_col]], ")"),
          paste0("LINESTRING(",
                 .data[[start_lon_col]], " ", .data[[start_lat_col]], ", ",
                 .data[[end_lon_col]], " ", .data[[end_lat_col]], ")")
        )
      }
    )
  
  # Convert to sf object
  sf_obj <- st_as_sf(df, wkt = "wkt_geometry", crs = 4326)
  
  # Apply buffer and return
  st_buffer(sf_obj, dist = buffer_dist)
}







#---------------------------------- Function to compute field replicate count
compute_field_replicates <- function(rep_string) {
  if (is.na(rep_string) || rep_string == "") return(NA_integer_)
  
  # Split on `/` to get number of samples
  parts <- strsplit(rep_string, "/")[[1]]
  
  # Count each pooled group (with `_`) or single sample as 1 unit
  replicate_units <- sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) == 1)) +  # unpooled
    sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) > 1))     # pooled
  
  return(replicate_units)
}














#---------------------------------- Function to compute PCR replicate count
compute_pcr_replicates <- function(rep_string) {
  if (is.na(rep_string) || rep_string == "") return(NA_integer_)
  
  # Split on `/` to get number of samples
  parts <- strsplit(rep_string, "/")[[1]]
  
  # Count each pooled group (with `_`) or single sample as 1 unit
  replicate_units <- sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) == 1)) +  # unpooled
    sum(sapply(parts, function(x) length(strsplit(x, "_")[[1]]) > 1))     # pooled
  
  return(12 * replicate_units)
}










#---------------------------------- Function : spatial_extraction --------------------
# This function extracts pixel values from a raster within a polygon (ie replicates buffer).
# It reprojects the polygon to the raster's CRS if necessary.
# It computes the mean, min, max, and range of the pixel values within the polygon.
# The mean is weighted by the area of the pixels intersected with the.


spatial_extraction <- function(var, rast, poly, stats = c("mean", "min", "max", "range")) {
  
  # Load necessary libraries
  require(terra)
  require(sf)
  require(dplyr)
  
  # Ensure stats contains only allowed values
  stats <- intersect(stats, c("mean", "min", "max", "range"))
  
  # Reproject polygons if CRS differs
  raster_crs <- sf::st_crs(terra::crs(rast))
  poly_crs <- sf::st_crs(poly)
  
  if (!identical(poly_crs, raster_crs)) {
    poly <- sf::st_transform(poly, crs = raster_crs)
  }
  
  # Initialize vectors to store the statistics
  rast_means <- numeric(length = nrow(poly))
  rast_mins <- numeric(length = nrow(poly))
  rast_maxs <- numeric(length = nrow(poly))
  rast_ranges <- numeric(length = nrow(poly))
  
  for (i in 1:nrow(poly)) {
    
    # Extract raster values within the current polygon with weights
    extracted_values <- terra::extract(x = rast, y = poly[i, ], exact = TRUE) 
    
    # Filter out NA values
    extracted_values <- extracted_values[!is.na(extracted_values[, 2]), ]
    
    values <- extracted_values[, 2]
    weights <- extracted_values[, 3]
    
    if (length(values) > 0) {
      if ("mean" %in% stats) {
        rast_means[i] <- sum(values * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
      }
      if ("min" %in% stats) {
        rast_mins[i] <- min(values, na.rm = TRUE)
      }
      if ("max" %in% stats) {
        rast_maxs[i] <- max(values, na.rm = TRUE)
      }
      if ("range" %in% stats) {
        rast_ranges[i] <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
      }
    } else {
      if ("mean" %in% stats) rast_means[i] <- NA
      if ("min" %in% stats)  rast_mins[i]  <- NA
      if ("max" %in% stats)  rast_maxs[i]  <- NA
      if ("range" %in% stats) rast_ranges[i] <- NA
    }
  }
  
  # Append results to poly
  if ("mean" %in% stats) poly[[paste0(var, "_mean")]] <- rast_means
  if ("min" %in% stats)  poly[[paste0(var, "_min")]]  <- rast_mins
  if ("max" %in% stats)  poly[[paste0(var, "_max")]]  <- rast_maxs
  if ("range" %in% stats) poly[[paste0(var, "_range")]] <- rast_ranges
  
  return(poly)
}



#---------------------------------- Function to compute SPECIFIC terra::terrain indices ------------------

compute_selected_terrain_ind <- function(rast, folder_path, neighbors, name, indices = c("slope", "aspect", "flowdir", "TPI", "TRI", "roughness")) {
  
  tile_name <- name  
  
  if ("slope" %in% indices) {
    filename_slope <- paste0(folder_path, "slope", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "slope", unit = "degrees", neighbors = neighbors, filename = filename_slope)
  }
  
  if ("aspect" %in% indices) {
    filename_aspect <- paste0(folder_path, "aspect", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "aspect", unit = "degrees", neighbors = neighbors, filename = filename_aspect)
  }
  
  if ("flowdir" %in% indices) {
    filename_flowdir <- paste0(folder_path, "flowdir", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "flowdir", neighbors = neighbors, filename = filename_flowdir)
  }
  
  if ("TPI" %in% indices) {
    filename_TPI <- paste0(folder_path, "TPI", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "TPI", neighbors = neighbors, filename = filename_TPI)
  }
  
  if ("TRI" %in% indices) {
    filename_TRI <- paste0(folder_path, "TRI", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "TRI", neighbors = neighbors, filename = filename_TRI)
  }
  
  if ("roughness" %in% indices) {
    filename_roughness <- paste0(folder_path, "roughness", neighbors, "_", tile_name, ".tif")
    terra::terrain(x = rast, v = "roughness", neighbors = neighbors, filename = filename_roughness)
  }
}

#---------------------------------- Function : identify_columns   ------------------
 # identify species and predictors columns from total dataframe 
identify_columns <- function(df) {
  # Identify response variables (species columns) based on the pattern "Genus.species"
  response_vars <- grep("^[A-Za-z]+\\.[A-Za-z]+$", colnames(df), value = TRUE)
  
  # Identify predictor variables as the columns that are not species columns
  predictor_vars <- setdiff(colnames(df), response_vars)
  
  # Get the column numbers for response_vars and predictor_vars
  response_columns <- which(colnames(df) %in% response_vars)
  predictor_columns <- which(colnames(df) %in% predictor_vars)
  
  # Return the result as a list containing the column numbers
  return(list(response_columns = response_columns, predictor_columns = predictor_columns))
}
#---------------------------------- Function to list species under a certain occurrence threshold --------------
list_rare_sp <- function(df, sum = NULL, exact = TRUE) {
  # Initialize an empty character vector to store column names
  result <- character()  # Using character() instead of list()
  
  # Loop through each column in df
  for (col_name in colnames(df)) {
    # Calculate the sum of the column
    column_sum <- sum(df[[col_name]], na.rm = TRUE)  # Handling NA values
    
    # Check the condition based on the parameters
    if (exact) {  # If exact is TRUE
      if (column_sum == sum) {
        result <- c(result, col_name)  # Append to character vector
      }
    } else {  # If exact is FALSE
      if (column_sum <= sum) {
        result <- c(result, col_name)  # Append to character vector
      }
    }
  }
  
  # Return the result (character vector of column names)
  return(result)
}

#----------------------------------  Function to list species columns and predictor columns
identify_columns <- function(df) {
  # Identify response variables (species columns) based on the pattern "Genus.species"
  response_vars <- grep("^[A-Za-z]+\\.[A-Za-z]+$", colnames(df), value = TRUE)
  
  # Identify predictor variables as the columns that are not species columns
  predictor_vars <- setdiff(colnames(df), response_vars)
  
  # Get the column numbers for response_vars and predictor_vars
  response_columns <- which(colnames(df) %in% response_vars)
  predictor_columns <- which(colnames(df) %in% predictor_vars)
  
  # Return the result as a list containing the column numbers
  return(list(response_columns = response_columns, predictor_columns = predictor_columns))
}

#---------------------------------- Function: plot_species_rarity --------------------
plot_species_rarity <- function(df,
                                threshold = 5,
                                col_to_del = c("spygen_code", "field_replicates", "replicates", "geom"),
                                round_digits = 0,
                                plot_type = c("bar", "hist", "both")) {
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the packages: dplyr, tidyr, ggplot2")
  }
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Validate plot type
  plot_type <- match.arg(plot_type)
  
  # Defensive: remove only existing columns
  col_to_del <- intersect(col_to_del, names(df))
  
  # ---- Data preparation ----
  plot_data <- df %>%
    dplyr::select(-all_of(col_to_del)) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "species", values_to = "tot") %>%
    mutate(
      rank = rank(-tot),
      above_mean = ifelse(tot > mean(tot, na.rm = TRUE), 1, 0),
      below_mean = ifelse(tot < mean(tot, na.rm = TRUE), 1, 0)
    ) %>%
    arrange(desc(tot))
  
  # ---- Summary statistics ----
  mean_df <- mean(plot_data$tot, na.rm = TRUE)
  percentage_above_mean <- (sum(plot_data$above_mean) / nrow(plot_data)) * 100
  percentage_below_mean <- (sum(plot_data$below_mean) / nrow(plot_data)) * 100
  
  highlighted_species <- plot_data %>% filter(tot < threshold)
  num_highlighted <- nrow(highlighted_species)
  percentage_highlighted <- (num_highlighted / nrow(plot_data)) * 100
  
  # ---- Plotting ----
  plots <- list()
  
  # --- (1) Bar plot ---
  if (plot_type %in% c("bar", "both")) {
    n_species <- nrow(plot_data)
    x_pos_base <- max(seq_len(n_species))
    shift <- min(34, max(0, x_pos_base - 1))
    x_pos <- x_pos_base - shift
    
    bar_plot <- ggplot(plot_data, aes(x = reorder(species, -tot), y = tot)) +
      geom_bar(stat = "identity", aes(fill = ifelse(tot < threshold, "Rare", "Common"))) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "darkslategrey", size = 0.5) +
      geom_hline(yintercept = mean_df, linetype = "dashed", color = "firebrick4", size = 0.5) +
      geom_text(
        x = x_pos,
        y = threshold,
        label = paste0("Rarity threshold: ", round(threshold, round_digits)),
        vjust = -1.5,
        hjust = 0,
        color = "darkslategrey",
        size = 3.8
      ) +
      geom_text(
        x = x_pos,
        y = mean_df,
        label = paste0("Mean occurrence: ", round(mean_df, round_digits)),
        vjust = -1.5,
        hjust = 0,
        color = "firebrick4",
        size = 3.8
      ) +
      scale_fill_manual(values = c("Rare" = "darkslategrey", "Common" = "darkslategray3"), name = "Species Type") +
      theme_test() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(color = "black", size = 10, hjust = 0.5),
        legend.position = "bottom"
      ) +
      labs(
        title = "Total Occurrences and Ranking for Each Species",
        x = "Species",
        y = "Total Occurrences",
        caption = paste0("Rare species: ≤ ", threshold, " occurrences. Total: ",
                         num_highlighted, "/", nrow(plot_data),
                         " species (", round(percentage_highlighted, 1), "%).")
      )
    
    plots$bar <- bar_plot
  }
  
  # --- (2) Histogram plot ---
  if (plot_type %in% c("hist", "both")) {
    hist_plot <- ggplot(plot_data, aes(x = tot)) +
      geom_histogram(bins = 10, fill = "lightblue", color = "grey40") +
      geom_vline(xintercept = threshold, color = "darkslategrey", linetype = "dashed") +
      geom_vline(xintercept = mean_df, color = "firebrick4", linetype = "dashed") +
      labs(
        title = "Distribution of Total Occurrences per Species",
        x = "Total Occurrences",
        y = "Frequency",
        caption = paste0("Mean: ", round(mean_df, 1), " | Threshold: ", threshold)
      ) +
      theme_minimal()
    plots$hist <- hist_plot
  }
  
  # ---- Return output ----
  out <- list(
    plot_data = plot_data,
    stats = list(
      mean_occurrence = mean_df,
      pct_above_mean = percentage_above_mean,
      pct_below_mean = percentage_below_mean,
      num_highlighted = num_highlighted,
      pct_highlighted = percentage_highlighted,
      threshold = threshold,
      n_species = nrow(plot_data)
    ),
    plots = plots
  )
  
  return(out)
}

#---------------------------------- vv --------------------
# Function: ggplot histogram with full summary (vector or data frame column)
gg_hist_summary <- function(x, col_name = NULL, bins = 50, 
                            fill_low = "#80ffdb", fill_high = "#03045e") {
  library(ggplot2)
  
  # If input is a data frame + col_name
  if (is.data.frame(x)) {
    if (is.null(col_name)) stop("Please provide 'col_name' when passing a data frame")
    if (!col_name %in% names(x)) stop("Column not found in data frame")
    vec <- as.numeric(x[[col_name]])
    plot_title <- paste("Distribution of", col_name)
    x_label <- col_name
  } else {
    # If input is a numeric vector
    vec <- as.numeric(x)
    plot_title <- "Distribution"
    x_label <- "Values"
  }
  
  # Remove NAs
  vec <- vec[!is.na(vec)]
  
  # Compute summary statistics
  mean_val <- mean(vec)
  sd_val <- sd(vec)
  median_val <- median(vec)
  min_val <- min(vec)
  max_val <- max(vec)
  n_val <- length(vec)
  
  # Prepare summary text
  summary_text <- paste0(
    "Mean: ", round(mean_val, 2), 
    "\nSD: ", round(sd_val, 2),
    "\nMedian: ", round(median_val, 2),
    "\nMin: ", round(min_val, 2),
    "\nMax: ", round(max_val, 2),
    "\nN: ", n_val
  )
  
  # Determine position for annotation (top-right with margin)
  x_pos <- max(vec, na.rm = TRUE)
  y_pos <- max(hist(vec, plot = FALSE)$density)  # use density scale
  y_pos <- y_pos * 0.95  # slightly below top
  
  # Build histogram
  hist_plot <- ggplot(data.frame(vec = vec), aes(x = vec)) +
    geom_histogram(aes(y = ..density.., fill = ..density..), 
                   bins = bins, color = "black", linewidth = 0.01) +
    scale_fill_gradient(low = fill_low, high = fill_high) +
    ggtitle(plot_title) +
    xlab(x_label) + ylab("Density") +
    theme_minimal() +
    annotate("text", x = x_pos, y = y_pos, label = summary_text,
             hjust = 1, vjust = 1, size = 4, color = "black")
  
  return(hist_plot)
}
gg_hist_summary <- function(x, col_name = NULL, bins = 50, 
                            fill_low = "#80ffdb", fill_high = "#03045e") {
  library(ggplot2)
  
  # Determine if input is vector or data frame column
  if (is.data.frame(x)) {
    if (is.null(col_name)) stop("Please provide 'col_name' when passing a data frame")
    if (!col_name %in% names(x)) stop("Column not found in data frame")
    vec <- as.numeric(x[[col_name]])
    x_label <- col_name
    plot_title <- paste("Distribution of", col_name)
  } else {
    vec <- as.numeric(x)
    if (is.null(col_name)) {
      x_label <- "Values"
      plot_title <- "Distribution"
    } else {
      x_label <- col_name
      plot_title <- paste("Distribution of", col_name)
    }
  }
  
  # Remove NAs
  vec <- vec[!is.na(vec)]
  
  # Skip empty vectors
  if (length(vec) == 0) {
    warning(paste("Skipping", col_name, "- no valid values"))
    return(NULL)
  }
  
  # Compute summary statistics
  mean_val <- mean(vec)
  sd_val <- sd(vec)
  median_val <- median(vec)
  min_val <- min(vec)
  max_val <- max(vec)
  n_val <- length(vec)
  
  summary_text <- paste0(
    "Mean: ", round(mean_val, 2), 
    "\nSD: ", round(sd_val, 2),
    "\nMedian: ", round(median_val, 2),
    "\nMin: ", round(min_val, 2),
    "\nMax: ", round(max_val, 2),
    "\nN: ", n_val
  )
  
  # Determine annotation position
  x_pos <- max(vec, na.rm = TRUE)
  y_pos <- max(hist(vec, plot = FALSE)$density) * 0.95
  
  # Build plot
  hist_plot <- ggplot(data.frame(x = vec), aes(x = x)) +
    geom_histogram(aes(y = ..density.., fill = ..density..), 
                   bins = bins, color = "black", linewidth = 0.01) +
    scale_fill_gradient(low = fill_low, high = fill_high) +
    ggtitle(plot_title) +
    xlab(x_label) + ylab("Density") +
    theme_minimal() +
    annotate("text", x = x_pos, y = y_pos, label = summary_text,
             hjust = 1, vjust = 1, size = 4, color = "black")
  
  return(hist_plot)
}
