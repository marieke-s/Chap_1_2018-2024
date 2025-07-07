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