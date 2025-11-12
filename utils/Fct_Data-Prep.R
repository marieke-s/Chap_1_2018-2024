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
        axis.text.x = element_text(angle = 45, hjust = 1, size = 3.8),
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

#---------------------------------- Function to plot hist + summary --------------------
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
  y_pos <- max(hist(vec, plot = FALSE)$counts) * 0.95
  
  # Build plot (using counts/frequency)
  hist_plot <- ggplot(data.frame(x = vec), aes(x = x)) +
    geom_histogram(aes(y = ..count.., fill = ..count..), 
                   bins = bins, color = "black", linewidth = 0.01) +
    scale_fill_gradient(low = fill_low, high = fill_high) +
    ggtitle(plot_title) +
    xlab(x_label) + ylab("Frequency") +
    theme_minimal() +
    annotate("text", x = x_pos, y = y_pos, label = summary_text,
             hjust = 1, vjust = 1, size = 4, color = "black")
  
  return(hist_plot)
}




#---------------------------------- Function to order spygen_code increasingly --------------------

reorder_spygen_codes <- function(x) {
  # If there's no underscore, return as is
  if (!grepl("_", x)) return(x)
  
  # Split into codes
  codes <- unlist(strsplit(x, "_"))
  
  # Extract numeric part for sorting
  nums <- as.numeric(sub("SPY", "", codes))
  
  # Reorder and rebuild
  ordered_codes <- codes[order(nums)]
  paste(ordered_codes, collapse = "_")
}


#---------------------------------- Function to reorder SPY within replicates groups ----------------
# Function to reorder SPY within within "_" groups, keeping "/" structure
reorder_replicates <- function(x) {
  # Split the string by "/" to handle each group separately
  groups <- unlist(strsplit(x, "/"))
  
  # Process each group independently
  reordered_groups <- sapply(groups, function(g) {
    # Only reorder if multiple SPY codes separated by "_"
    if (grepl("_", g)) {
      codes <- unlist(strsplit(g, "_"))
      nums <- as.numeric(sub("SPY", "", codes))
      codes <- codes[order(nums)]
      paste(codes, collapse = "_")
    } else {
      g
    }
  })
  
  # Recombine the groups using "/"
  paste(reordered_groups, collapse = "/")
}
#---------------------------------- Map + histogram for numeric vars (factor-aware, tiny legends, full stats) ------------------

map_index_plots <- function(df,
                            cols_to_plot = NULL,
                            output_directory = "./figures/Div_indices/Maps/",
                            version = "v1",
                            bins = 50,
                            x_col = "x",
                            y_col = "y",
                            factor = NULL,               # split into panels by this categorical column (e.g., "season")
                            factor_levels = NULL,        # optional explicit order for factor levels
                            panel_ncol = NULL,           # columns in panel grid (auto if NULL)
                            panel_text_scale = c("auto","xs","s","m"),  # control text sizes in panels
                            stats_in_panel_title = FALSE # also show compact stats in each panel's map title
) {
  panel_text_scale <- match.arg(panel_text_scale)
  
  # ---- libs
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(gridExtra)
  library(rlang)
  library(rnaturalearth)
  library(grid)
  
  # ---- text-size presets
  single_sizes <- list(title=14, axis_title=11, axis_text=9, annot=3.8, legend_title=10, legend_text=9)
  panel_sizes_auto <- function(n) {
    if (n <= 2)   list(title=12, axis_title=9,  axis_text=8,   annot=3.2, legend_title=8,  legend_text=7)
    else if (n<=4)list(title=10, axis_title=8,  axis_text=7,   annot=3.0, legend_title=7,  legend_text=6.5)
    else if (n<=6)list(title=9,  axis_title=7,  axis_text=6.5, annot=2.8, legend_title=6,  legend_text=6)
    else          list(title=8,  axis_title=6.5,axis_text=6,   annot=2.4, legend_title=5,  legend_text=5)
  }
  panel_sizes_named <- function(scale){
    switch(scale,
           xs=list(title=6.5, axis_title=5.5, axis_text=5,   annot=2.0, legend_title=4, legend_text=4),
           s =list(title=8.0, axis_title=6.5, axis_text=6.0, annot=2.1, legend_title=5, legend_text=5),
           m =list(title=9.0, axis_title=7.0, axis_text=6.5, annot=2.8, legend_title=6, legend_text=6),
           list(title=8.0, axis_title=6.5, axis_text=6.0, annot=2.6, legend_title=5, legend_text=5))
  }
  
  # ---- helpers
  .short_stats <- function(v) {
    v <- suppressWarnings(as.numeric(v)); v <- v[!is.na(v)]
    if (!length(v)) return("N=0")
    paste0("N=", length(v),
           " • μ=", round(mean(v),2),
           " • σ=", round(sd(v),2),
           " • med=", round(median(v),2),
           " • min=", round(min(v),2),
           " • max=", round(max(v),2))
  }
  
  gg_hist_summary <- function(x, col_name = NULL, bins = 50,
                              fill_low = "#80ffdb", fill_high = "#03045e",
                              sizes = single_sizes) {
    if (is.data.frame(x)) {
      if (is.null(col_name)) stop("Please provide 'col_name' when passing a data frame")
      if (!col_name %in% names(x)) stop(paste0("Column '", col_name, "' not found"))
      vec <- suppressWarnings(as.numeric(x[[col_name]]))
      xlab_txt <- col_name
      ttl_txt  <- paste("Distribution of", col_name)
    } else {
      vec <- suppressWarnings(as.numeric(x))
      xlab_txt <- if (is.null(col_name)) "Values" else col_name
      ttl_txt  <- if (is.null(col_name)) "Distribution" else paste("Distribution of", col_name)
    }
    vec <- vec[!is.na(vec)]
    if (!length(vec)) return(NULL)
    
    # stats text
    summary_text <- .short_stats(vec)
    
    # headroom to avoid overlap; annotate top-left
    h <- hist(vec, plot = FALSE, breaks = bins)
    y_max <- max(h$counts, 1)
    y_lim <- c(0, y_max * 1.35)
    x_pos <- min(vec, na.rm = TRUE)
    y_pos <- y_max * 1.24
    
    ggplot(data.frame(x = vec), aes(x = x)) +
      geom_histogram(aes(y = after_stat(count), fill = after_stat(count)),
                     bins = bins, color = "black", linewidth = 0.01) +
      scale_fill_gradient(low = fill_low, high = fill_high) +
      ggtitle(ttl_txt) +
      xlab(xlab_txt) + ylab("Frequency") +
      coord_cartesian(ylim = y_lim) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = sizes$title, face = "bold"),
        axis.title = element_text(size = sizes$axis_title),
        axis.text  = element_text(size = sizes$axis_text),
        legend.position = "none"
      ) +
      annotate("text", x = x_pos, y = y_pos, label = summary_text,
               hjust = 0, vjust = 0, size = sizes$annot, color = "black")
  }
  
  # ---- output dir & basemap
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # ---- coordinates
  is_sf <- inherits(df, "sf")
  df_local <- df
  if (is_sf) {
    if (!all(c(x_col, y_col) %in% names(df_local))) {
      coords <- sf::st_coordinates(df_local)
      df_local <- df_local %>%
        mutate(!!x_col := coords[, 1],
               !!y_col := coords[, 2])
    }
  } else {
    if (!all(c(x_col, y_col) %in% names(df_local))) {
      stop(paste0("x_col='", x_col, "' and y_col='", y_col,
                  "' must exist in 'df' when df is not an sf object."))
    }
  }
  
  # ---- columns to plot (numeric)
  if (is.null(cols_to_plot)) {
    num_cols <- df_local %>%
      { if (is_sf) st_drop_geometry(.) else . } %>%
      dplyr::select(where(is.numeric)) %>%
      colnames()
    cols_to_plot <- setdiff(num_cols, c(x_col, y_col))
  } else {
    missing_cols <- setdiff(cols_to_plot, names(df_local))
    if (length(missing_cols)) {
      stop("These cols_to_plot are missing in df: ", paste(missing_cols, collapse = ", "))
    }
  }
  
  # ---- factor handling
  use_factor <- !is.null(factor)
  if (use_factor) {
    if (!factor %in% names(df_local)) stop("factor column '", factor, "' not found in df.")
    if (!is.null(factor_levels)) df_local[[factor]] <- factor(df_local[[factor]], levels = factor_levels)
    else df_local[[factor]] <- as.factor(df_local[[factor]])
    factor_vals <- levels(droplevels(df_local[[factor]]))
    factor_vals <- factor_vals[!is.na(factor_vals)]
    if (!length(factor_vals)) stop("factor column '", factor, "' has no non-NA levels.")
  }
  
  # ---- iterate
  for (col in cols_to_plot) {
    
    if (!use_factor) {
      # single plot sizing
      sizes <- single_sizes
      
      hist_plot <- gg_hist_summary(df_local, col_name = col, bins = bins, sizes = sizes)
      if (is.null(hist_plot)) next
      
      map_plot <- ggplot() +
        geom_sf(data = world, fill = "white", color = "lightblue") +
        geom_point(
          data = df_local,
          aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[col]]),
          size = 1, alpha = 0.8
        ) +
        scale_color_gradient(low = "yellow", high = "red") +
        ggtitle(paste("Map of", col)) +
        xlab("Longitude") + ylab("Latitude") +
        coord_sf(
          xlim = c(min(df_local[[x_col]], na.rm = TRUE) - 0.1,
                   max(df_local[[x_col]], na.rm = TRUE) + 0.0),
          ylim = c(min(df_local[[y_col]], na.rm = TRUE) - 0.0,
                   max(df_local[[y_col]], na.rm = TRUE) + 0.0)
        ) +
        theme_test() +
        theme(
          panel.background = element_rect(fill = "#e0edff"),
          plot.title = element_text(size = sizes$title, face = "bold"),
          axis.title = element_text(size = sizes$axis_title),
          axis.text  = element_text(size = sizes$axis_text),
          legend.title = element_text(size = sizes$legend_title),
          legend.text  = element_text(size = sizes$legend_text)
        )
      
      combined_plot <- grid.arrange(map_plot, hist_plot,
                                    layout_matrix = rbind(c(1), c(2)),
                                    heights = c(2, 1))
      print(combined_plot)
      ggsave(file.path(output_directory, paste0("map_", col, "_", version, ".png")),
             plot = combined_plot, width = 8, height = 10, dpi = 300)
      
    } else {
      # panel branch
      n_levels <- length(factor_vals)
      sizes <- if (panel_text_scale == "auto") panel_sizes_auto(n_levels) else panel_sizes_named(panel_text_scale)
      
      # global title (one per numeric var)
      global_title <- grid::textGrob(
        paste("Map and Distribution of", col),
        gp = gpar(fontface = "bold", fontsize = sizes$title + 2)
      )
      
      level_grobs <- list()
      for (lev in factor_vals) {
        df_sub <- dplyr::filter(df_local, .data[[factor]] == lev)
        
        # histogram for this level
        hist_plot <- gg_hist_summary(df_sub, col_name = col, bins = bins, sizes = sizes)
        if (is.null(hist_plot)) next
        
        # small per-panel map title; optional compact stats appended
        sub_title <- paste(factor, ":", lev)
        if (stats_in_panel_title) sub_title <- paste0(sub_title, "  (", .short_stats(df_sub[[col]]), ")")
        
        map_plot <- ggplot() +
          geom_sf(data = world, fill = "white", color = "lightblue") +
          geom_point(
            data = df_sub,
            aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[col]]),
            size = 0.9, alpha = 0.8
          ) +
          scale_color_gradient(low = "yellow", high = "red") +
          ggtitle(sub_title) +
          xlab("Longitude") + ylab("Latitude") +
          coord_sf(
            xlim = c(min(df_sub[[x_col]], na.rm = TRUE) - 0.1,
                     max(df_sub[[x_col]], na.rm = TRUE) + 0.0),
            ylim = c(min(df_sub[[y_col]], na.rm = TRUE) - 0.0,
                     max(df_sub[[y_col]], na.rm = TRUE) + 0.0)
          ) +
          theme_test() +
          theme(
            panel.background = element_rect(fill = "#e0edff"),
            plot.title = element_text(size = sizes$title, face = "bold"),
            axis.title = element_text(size = sizes$axis_title),
            axis.text  = element_text(size = sizes$axis_text),
            legend.title = element_text(size = sizes$legend_title),  # VERY small
            legend.text  = element_text(size = sizes$legend_text)
          )
        
        combo <- grid.arrange(map_plot, hist_plot,
                              layout_matrix = rbind(c(1), c(2)),
                              heights = c(2, 1))
        level_grobs[[lev]] <- combo
      }
      
      if (!length(level_grobs)) next
      ncol_use <- if (!is.null(panel_ncol)) panel_ncol else ceiling(sqrt(length(level_grobs)))
      panel_grid <- arrangeGrob(grobs = level_grobs, ncol = ncol_use)
      panel_all  <- gridExtra::arrangeGrob(global_title, panel_grid, ncol = 1, heights = c(0.5, 10))
      
      print(panel_all)
      
      panel_width <- max(8, 4 * ncol_use)
      ggsave(file.path(output_directory, paste0("panel_", col, "_by_", factor, "_", version, ".png")),
             plot = panel_all, width = panel_width, height = 10, dpi = 300)
    }
  }
  
  message("Saved to: ", normalizePath(output_directory))
  invisible(NULL)
}

#---------------------------------- Function to map + hist + summary categorical var --------------------
# map_categorical_plots(): map + bar (mode summary) for categorical columns
# Adds: separate_maps = TRUE to facet one map per category level


map_categorical_plots <- function(df,
                                  cols_to_plot = NULL,
                                  output_directory = "./figures/Div_indices/Maps/",
                                  version = "v1",
                                  x_col = "x",
                                  y_col = "y",
                                  top_n_bars = 20,         # show top-N categories in bar chart
                                  legend_limit = 20,       # hide legend on single-map if too many levels
                                  separate_maps = FALSE    # facet one map per category
) {
  # ---- libs
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(gridExtra)
  library(rlang)
  library(rnaturalearth)
  library(forcats)
  
  # ---- helpers ------------------------------------------------
  get_modes <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(character(0))
    tab <- sort(table(x), decreasing = TRUE)
    max_ct <- tab[1]
    names(tab)[tab == max_ct]
  }
  
  gg_bar_mode <- function(x, col_name, top_n_bars = 20,
                          fill_low = "#80ffdb", fill_high = "#03045e") {
    stopifnot(col_name %in% names(x))
    vec <- x[[col_name]]
    # coerce to factor for consistent plotting
    if (!is.factor(vec)) vec <- as.factor(vec)
    vec <- fct_explicit_na(vec, na_level = "(NA)")
    
    counts <- as.data.frame(table(vec), stringsAsFactors = FALSE)
    colnames(counts) <- c("level", "n")
    if (!nrow(counts) || all(counts$n == 0)) {
      warning(paste("Skipping", col_name, "- no valid categories"))
      return(NULL)
    }
    
    # top-N handling for readability
    counts <- counts |>
      arrange(desc(n))
    if (nrow(counts) > top_n_bars) {
      top <- counts[1:top_n_bars, ]
      other_sum <- sum(counts$n[(top_n_bars + 1):nrow(counts)])
      counts <- rbind(top, data.frame(level = "Other", n = other_sum))
    }
    counts$level <- factor(counts$level, levels = counts$level)
    
    # modes (could be multiple)
    modes <- get_modes(x[[col_name]])
    mode_txt <- if (length(modes) == 0) {
      "Mode: (none)"
    } else if (length(modes) == 1) {
      paste0("Mode: ", modes[1])
    } else {
      paste0("Modes: ", paste(modes, collapse = ", "))
    }
    n_non_na <- sum(!is.na(x[[col_name]]))
    summary_text <- paste0(mode_txt, "\nN (non-NA): ", n_non_na)
    
    # y cap near max
    y_max <- max(counts$n, 1)
    y_lim <- c(0, y_max * 1.05)
    y_pos <- y_max * 0.98
    
    ggplot(counts, aes(x = level, y = n, fill = n)) +
      geom_col(color = "black", linewidth = 0.1) +
      scale_fill_gradient(low = fill_low, high = fill_high) +
      ggtitle(paste("Distribution of", col_name)) +
      xlab(col_name) + ylab("Frequency") +
      coord_cartesian(ylim = y_lim) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      annotate("text", x = Inf, y = y_pos, label = summary_text,
               hjust = 1.02, vjust = 1, size = 4, color = "black")
  }
  # -------------------------------------------------------------
  
  # ---- output dir
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
  
  # ---- world basemap
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # ---- ensure coordinates
  is_sf <- inherits(df, "sf")
  df_local <- df
  
  if (is_sf) {
    has_xy <- all(c(x_col, y_col) %in% names(df_local))
    if (!has_xy) {
      coords <- sf::st_coordinates(df_local)
      df_local <- df_local %>%
        mutate(!!x_col := coords[, 1],
               !!y_col := coords[, 2])
    }
  } else {
    if (!all(c(x_col, y_col) %in% names(df_local))) {
      stop(paste0("x_col='", x_col, "' and y_col='", y_col,
                  "' must exist in 'df' when df is not an sf object."))
    }
  }
  
  # ---- columns to plot (categorical)
  if (is.null(cols_to_plot)) {
    base <- if (is_sf) sf::st_drop_geometry(df_local) else df_local
    cat_cols <- names(base)[vapply(base, function(col)
      is.factor(col) || is.character(col) || is.logical(col), logical(1))]
    cols_to_plot <- setdiff(cat_cols, c(x_col, y_col))
  } else {
    missing_cols <- setdiff(cols_to_plot, names(df_local))
    if (length(missing_cols)) {
      stop("These cols_to_plot are missing in df: ", paste(missing_cols, collapse = ", "))
    }
  }
  
  # ---- iterate
  for (col in cols_to_plot) {
    # Bar chart with mode (uses top_n_bars + "Other" if needed)
    bar_plot <- gg_bar_mode(df_local, col_name = col, top_n_bars = top_n_bars)
    if (is.null(bar_plot)) next
    
    # SINGLE MAP (default): color by category
    lev_ct <- nlevels(as.factor(df_local[[col]]))
    hide_legend <- lev_ct > legend_limit
    
    if (!separate_maps) {
      map_plot <- ggplot() +
        geom_sf(data = world, fill = "white", color = "lightblue") +
        geom_point(
          data = df_local,
          aes(x = .data[[x_col]], y = .data[[y_col]], color = as.factor(.data[[col]])),
          size = 1, alpha = 0.8, na.rm = TRUE
        ) +
        scale_color_discrete(name = col) +
        ggtitle(paste("Map of", col)) +
        xlab("Longitude") + ylab("Latitude") +
        coord_sf(
          xlim = c(min(df_local[[x_col]], na.rm = TRUE) - 0.1,
                   max(df_local[[x_col]], na.rm = TRUE) + 0.0),
          ylim = c(min(df_local[[y_col]], na.rm = TRUE) - 0.0,
                   max(df_local[[y_col]], na.rm = TRUE) + 0.0)
        ) +
        theme_test() +
        theme(
          panel.background = element_rect(fill = "#e0edff"),
          legend.position = if (hide_legend) "none" else "right"
        )
      
      combined_plot <- grid.arrange(
        map_plot, bar_plot,
        layout_matrix = rbind(c(1), c(2)),
        heights = c(2, 1)
      )
      print(combined_plot)
      out_path <- file.path(output_directory, paste0("mapcat_", col, "_", version, ".png"))
      ggsave(filename = out_path, plot = combined_plot, width = 8, height = 10, dpi = 300)
      
    } else {
      # SEPARATE MAPS: facet one map per category (limited to top_n_bars most frequent)
      # Determine top levels to facet (same criterion as bar chart, excluding "Other")
      base_vec <- df_local[[col]]
      base_vec <- fct_explicit_na(as.factor(base_vec), na_level = "(NA)")
      counts_all <- as.data.frame(table(base_vec), stringsAsFactors = FALSE) |>
        arrange(desc(Freq))
      keep_levels <- counts_all$base_vec[seq_len(min(nrow(counts_all), top_n_bars))]
      
      df_panel <- df_local %>%
        mutate(.cat = fct_explicit_na(as.factor(.data[[col]]), na_level = "(NA)")) %>%
        filter(.cat %in% keep_levels)
      
      # pick grid size: ~square
      n_levels <- nlevels(factor(df_panel$.cat))
      ncol_facets <- ceiling(sqrt(n_levels))
      
      map_panel <- ggplot() +
        geom_sf(data = world, fill = "white", color = "lightblue") +
        geom_point(
          data = df_panel,
          aes(x = .data[[x_col]], y = .data[[y_col]]),
          size = 1, alpha = 0.8, na.rm = TRUE
        ) +
        ggtitle(paste("Maps by category -", col)) +
        xlab("Longitude") + ylab("Latitude") +
        coord_sf(
          xlim = c(min(df_panel[[x_col]], na.rm = TRUE) - 0.1,
                   max(df_panel[[x_col]], na.rm = TRUE) + 0.0),
          ylim = c(min(df_panel[[y_col]], na.rm = TRUE) - 0.0,
                   max(df_panel[[y_col]], na.rm = TRUE) + 0.0)
        ) +
        facet_wrap(~ .cat, ncol = ncol_facets) +
        theme_test() +
        theme(
          panel.background = element_rect(fill = "#e0edff"),
          legend.position = "none",
          strip.text = element_text(size = 9)
        )
      
      # Combine panel + bar
      combined_panel <- grid.arrange(
        map_panel, bar_plot,
        layout_matrix = rbind(c(1), c(2)),
        heights = c(2, 1)
      )
      print(combined_panel)
      
      # Width scaled to facet columns
      panel_width <- max(8, 3.5 * ncol_facets)
      ggsave(
        filename = file.path(output_directory, paste0("mapcatpanel_", col, "_", version, ".png")),
        plot = combined_panel, width = panel_width, height = 10, dpi = 300
      )
    }
  }
  
  message("Saved to: ", normalizePath(output_directory))
}

#---------------------------------- Function to export count tables --------------------


# export_count_tables(): export counts + % per category as a figure 
export_count_tables <- function(mtdt,
                                cat_cols = c("year", "month", "season"),
                                output_directory = "./figures/Summaries/",
                                filename_stub = "sample_counts",
                                version = "v1",
                                layout = c("vertical", "horizontal"),
                                width = NULL, height = NULL, dpi = 300) {
  layout <- match.arg(layout)
  # ---- libs
  library(dplyr)
  library(sf)
  library(rlang)
  library(grid)
  library(gridExtra)
  library(ggplot2)   # for ggsave
  
  # ---- prepare output directory
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
  
  # ---- drop geometry if sf
  base <- if (inherits(mtdt, "sf")) sf::st_drop_geometry(mtdt) else mtdt
  
  # ---- compute counts + percentage
  count_list <- lapply(cat_cols, function(col) {
    if (!col %in% names(base)) {
      warning(paste0("Skipping '", col, "': not found"))
      return(NULL)
    }
    tbl <- base %>%
      group_by(.data[[col]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(.data[[col]]) %>%
      rename(!!col := .data[[col]]) %>%
      mutate(`%` = round(100 * count / sum(count), 1))
    tbl
  })
  names(count_list) <- cat_cols
  count_list <- Filter(Negate(is.null), count_list)
  if (length(count_list) == 0) stop("No valid columns to tabulate.")
  
  # ---- define table theme
  tab_theme <- gridExtra::ttheme_minimal(
    core = list(fg_params = list(cex = 0.9), bg_params = list(fill = c(NA, "#f7f7f7"))),
    colhead = list(fg_params = list(fontface = 2), bg_params = list(fill = "#eaeaea"))
  )
  
  # ---- build table grobs (no titles)
  make_table_panel <- function(tbl) {
    gridExtra::tableGrob(tbl, rows = NULL, theme = tab_theme)
  }
  panels <- lapply(count_list, make_table_panel)
  
  # ---- arrange panels
  n <- length(panels)
  if (layout == "horizontal") {
    ncol <- n; nrow <- 1
  } else {
    ncol <- 1; nrow <- n
  }
  fig <- gridExtra::arrangeGrob(grobs = panels, ncol = ncol, nrow = nrow)
  
  # ---- auto size if not provided
  if (is.null(width) || is.null(height)) {
    if (layout == "horizontal") {
      width  <- if (is.null(width)) 3.5 * n  else width
      height <- if (is.null(height)) 4.5     else height
    } else {
      width  <- if (is.null(width)) 5.0      else width
      height <- if (is.null(height)) 2.8 * n else height
    }
  }
  
  # ---- save figure
  out_file <- file.path(
    output_directory,
    paste0(filename_stub, "_", version, ".png")
  )
  ggsave(filename = out_file, plot = fig, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", normalizePath(out_file))
  invisible(out_file)
}


#---------------------------------- Function to identify outliers using the 1.5 * IQR rule --------------
# see here for method : https://www.khanacademy.org/math/statistics-probability/summarizing-quantitative-data/box-whisker-plots/a/identifying-outliers-iqr-rule

find_outliers <- function(x) {
  if (!is.numeric(x)) return(FALSE)  # Skip non-numeric columns
  
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  
  any(x < lower_bound | x > upper_bound, na.rm = TRUE)  # Returns TRUE if outliers exist
}
#---------------------------------- Function to plot variables distrbution ----------------
plot_distributions <- function(data, 
                               n_per_page = 8, 
                               hist_color = "lightblue", 
                               output_path = NULL, 
                               version_number = 1) {
  n_predictors <- ncol(data)
  n_pages <- ceiling(n_predictors / n_per_page)
  
  for (page in 1:n_pages) {
    start <- (page - 1) * n_per_page + 1
    end <- min(page * n_per_page, n_predictors)
    predictors <- colnames(data)[start:end]
    
    # Create filename for this page
    if (!is.null(output_path)) {
      file_name <- file.path(output_path, 
                             paste0("distribution_v", version_number, "_page", page, ".png"))
      png(filename = file_name, width = 1600, height = 900)
    }
    
    # Set up a grid layout (2 rows, 4 columns)
    par(mfrow = c(2, 4))
    
    # Plot histograms
    for (predictor in predictors) {
      hist(data[[predictor]], 
           main = predictor, 
           xlab = predictor,
           col = hist_color, 
           border = "white", 
           breaks = 20)
    }
    
    # Close the plotting device if outputting to file
    if (!is.null(output_path)) {
      dev.off()
      message("Saved: ", file_name)
    } else {
      # Pause for manual browsing if not saving
      readline(prompt = "Press [Enter] to see the next set of plots...")
    }
  }
}
