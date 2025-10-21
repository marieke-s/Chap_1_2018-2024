# Pipeline with exactextract

def get_dates(date, time_step):
    """
    Calculate the range of dates for a given time step relative to the provided date.
    """
    from datetime import datetime, timedelta

    if isinstance(date, str):
        date = datetime.strptime(date, "%Y-%m-%d")

    end_date = date - timedelta(days=1)

    time_deltas = {
        'day': 1,
        'week': 7,
        'month': 30,
        'year': 365,
        '5years': 5 * 365 + 1,
    }

    if time_step not in time_deltas:
        raise ValueError("Unsupported time step. Choose from 'day', 'week', 'month', 'year', or '5years'.")

    start_date = date - timedelta(days=time_deltas[time_step])
    return start_date, end_date






def compute_stats(data_array, shape_geometry):
    """
    Extract values and compute weighted mean - automatically with exactextract- , min and max.
    """
    
    try:
        feature = {"type": "Feature", "geometry": mapping(shape_geometry), "properties": {}}
        res = exact_extract(data_array, [feature], ["mean", "min", "max"])


        if not res or len(res) == 0:
            return None, None, None

        props = res[0]["properties"]

        # Multi-band keys: band_1_mean, band_2_mean, etc.
        mean_vals = [v for k, v in props.items() if k.endswith("_mean")]
        min_vals  = [v for k, v in props.items() if k.endswith("_min")]
        max_vals  = [v for k, v in props.items() if k.endswith("_max")]

         # Single-band keys: mean, min, max
        if not mean_vals:
            if "mean" in props:
                mean_vals = [props["mean"]]
        if not min_vals:
            if "min" in props:
                min_vals = [props["min"]]
        if not max_vals:
            if "max" in props:
                max_vals = [props["max"]]
   

        if not mean_vals or not min_vals or not max_vals:
            return None, None, None

        mean_val = float(np.nanmean(mean_vals))
        min_val  = float(np.nanmin(min_vals))
        max_val  = float(np.nanmax(max_vals))

        return mean_val, min_val, max_val

    except Exception as e:
        print(f"compute_stats ERROR: {e}")
        return None, None, None







def open_nc(shape_geometry, date, netcdf_path, variable="CHL"):
    """
    Compute NCDF statistics for a given geometry and date using a netCDF file.
    """
    results = {}

    try:
        if isinstance(date, str):
            date = datetime.strptime(date, "%Y-%m-%d")

        ds = xr.open_dataset(netcdf_path)
        ds = ds.rio.write_crs("EPSG:4326", inplace=True)


        target_date = date - timedelta(days=1)
        time_steps = ["day", "week", "month", "year", "5years"]
        date_ranges = {label: get_dates(date, label) for label in time_steps}


        for label, (start_date, end_date) in date_ranges.items():
            ds_time_range = ds.sel(time=slice(start_date, end_date))

            if ds_time_range.time.size == 0:
                results[label] = (None, None, None, 0)
                continue

            chl_data = ds_time_range[variable]
            valid_data = chl_data.dropna(dim="time", how="all")


            if valid_data.size > 0:                 
                mean_val, min_val, max_val = compute_stats(valid_data, shape_geometry)                
                results[label] = (mean_val, min_val, max_val)
            else:
                print("valid_data.size == 0")
                results[label] = (None, None, None)

        return results

    except Exception as e:
        print(f"Error processing shape with target date: {date}: {e}")
        return {}





def process_geojson(geojson_path, netcdf_path, output_path, variable="CHL"):
    """
    Process the GeoJSON file and compute statistics for each shape using a netCDF file.
    """
    shapes = gpd.read_file(geojson_path)
    shapes = shapes.set_crs("EPSG:4326", allow_override=True)
    shapes = shapes[0:20]

    results = []

    for _, row in tqdm(shapes.iterrows(), total=shapes.shape[0], desc="Processing shapes"):
        shape_geometry = row.geometry
        date = row["date"]
        polygon_id = row.get("replicates", None)

        nc_stats = open_nc(shape_geometry, date, netcdf_path, variable)

        print("------------replicates--------------")
        print(polygon_id)
        print("------------nc_stats--------------")
        print(nc_stats)


        result_entry = {"replicates": polygon_id}
        for label, (mean, min_val, max_val) in nc_stats.items():
            result_entry[f"Cop_CHL_{label}_mean"] = mean
            result_entry[f"Cop_CHL_{label}_min"] = min_val
            result_entry[f"Cop_CHL_{label}_max"] = max_val

        results.append(result_entry)

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)
