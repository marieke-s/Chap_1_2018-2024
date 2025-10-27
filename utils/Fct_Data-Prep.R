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

# Function to species columns and predictor columns
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
