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
#---------------------------------- Function : combine_rows----------------
combine_rows <- function(data, codes_to_combine, new_code, proceed = TRUE) {
  # Ensure exactly two codes
  if (length(codes_to_combine) != 2) {
    stop("Please provide exactly two spygen_codes to combine.")
  }
  
  # Get the rows to combine
  rows_to_combine <- data[data$spygen_code %in% codes_to_combine, ]
  
  if (nrow(rows_to_combine) != 2) {
    stop("The specified spygen_codes must match exactly two rows in the data.")
  }
  
  # Set up columns to compare
  other_columns <- setdiff(names(rows_to_combine), "spygen_code")
  
  # Compare values
  differing_columns <- which(rows_to_combine[1, other_columns] != rows_to_combine[2, other_columns])
  differing_col_names <- other_columns[differing_columns]
  
  # Get row IDs for message
  code_1 <- rows_to_combine$spygen_code[1]
  code_2 <- rows_to_combine$spygen_code[2]
  
  # Warn if differences exist
  if (length(differing_columns) > 0) {
    warning(sprintf(
      "Warning: Values differ between %s and %s in columns: %s",
      code_1,
      code_2,
      paste(differing_col_names, collapse = ", ")
    ))
    if (!proceed) stop("Operation aborted by the user.")
  }
  
  # Use the first row as the base
  new_row <- rows_to_combine[1, ]
  new_row$spygen_code <- new_code
  
  # Custom handling for differing 'comments'
  if ("comments" %in% differing_col_names) {
    comment_1 <- as.character(rows_to_combine$comments[1])
    comment_2 <- as.character(rows_to_combine$comments[2])
    
    comment_1 <- ifelse(is.na(comment_1), "", comment_1)
    comment_2 <- ifelse(is.na(comment_2), "", comment_2)
    
    if (nzchar(comment_2)) {
      combined_comment <- trimws(paste0(
        comment_1,
        ifelse(nzchar(comment_1), " | ", ""),
        code_2, ": ", comment_2
      ))
      new_row$comments <- combined_comment
    }
  }
  
  # Remove original rows and add new one
  data <- data[!data$spygen_code %in% codes_to_combine, ]
  data <- rbind(data, new_row)
  
  return(data)
}
