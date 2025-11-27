#---------- SET UP ----------------    

import calendar
import dask
from datetime import datetime, timedelta
import exactextract as ee
from exactextract import exact_extract
import gc
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
import rioxarray as rxr
from shapely.geometry import mapping, shape
from shapely.geometry import mapping, Point
from scipy.spatial import cKDTree
import time 
from tqdm import tqdm
import xarray as xr
import os
import glob
import re

















#---------- FUNCTIONS ----------------    


# Initial extraction optimised for memory

def get_dates(date, time_step):
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
        raise ValueError("Unsupported time step.")
    start_date = date - timedelta(days=time_deltas[time_step])
    return start_date, end_date






def compute_stats_from_2d_raster(raster_2d, shape_geometry, agg_name="mean"):
    """
    raster_2d: xarray.DataArray 2D (lat x lon) ou numpy array avec méta spatiale portée par l'objet
    shape_geometry: shapely geometry
    Returns: (mean, min, max) pour le polygone basé sur exact_extract
    """

    try:
        # exact_extract a besoin d'une géométrie GeoJSON
        feature = {"type": "Feature", "geometry": mapping(shape_geometry), "properties": {}}

        # On laisse raster_2d tel quel, car s'il s'agit d'un DataArray xarray
        # il porte les infos spatiales nécessaires (transform, CRS, etc.)
        res = exact_extract(raster_2d, [feature], ["mean", "min", "max"])

        # Rien renvoyé -> on retourne NaN
        if not res or len(res) == 0:
            return np.nan, np.nan, np.nan

        # Certains bindings renvoient une clé "properties", d'autres directement les stats
        props_container = res[0]
        props = props_container.get("properties", props_container)

        def collect_stats(suffix, base_key):
            """
            Récupère les valeurs associées à *_suffix, ou à base_key,
            en filtrant les None et NaN.
            """
            vals = [v for k, v in props.items()
                    if k.endswith(suffix) and v is not None]

            # fallback: clé simple ("mean", "min", "max")
            if not vals and base_key in props and props[base_key] is not None:
                vals = [props[base_key]]

            if not vals:
                return np.array([])

            # Conversion en float + filtrage des NaN
            arr = np.array(vals, dtype=float)
            arr = arr[~np.isnan(arr)]
            return arr

        mean_vals = collect_stats("_mean", "mean")
        min_vals  = collect_stats("_min", "min")
        max_vals  = collect_stats("_max", "max")

        # Si aucune valeur valide
        if mean_vals.size == 0 and min_vals.size == 0 and max_vals.size == 0:
            return np.nan, np.nan, np.nan

        mean_val = float(mean_vals.mean()) if mean_vals.size > 0 else np.nan
        min_val  = float(min_vals.min())  if min_vals.size  > 0 else np.nan
        max_val  = float(max_vals.max())  if max_vals.size  > 0 else np.nan

        return mean_val, min_val, max_val

    except Exception as e:
        # Log minimal + retour NaN
        print(f"compute_stats_from_2d_raster ERROR: {e}")
        return np.nan, np.nan, np.nan







def process_geojson(geojson_path, netcdf_path, output_path, variable="CHL",
                              time_steps=None, xr_chunks=None):
    """
    Memory-optimized processing:
      - open netCDF once (with dask chunks)
      - for each polygon: for each time window compute temporal aggregates (mean/min/max) -> pass 2D to exact_extract
    Parameters:
      xr_chunks: dictionary passed to xr.open_dataset(..., chunks=xr_chunks)
                 e.g. {'time': 10, 'lat': 512, 'lon': 512}
    """

    if time_steps is None:
        time_steps = ["day", "week", "month", "year", "5years"]

    shapes = gpd.read_file(geojson_path)
    # shapes = shapes[10:11] # a modifier : for test
    shapes = shapes.set_crs("EPSG:4326", allow_override=True)

    # open dataset once, with dask chunking
    if xr_chunks is None:
        xr_chunks = {'time': 10}  # minimal sensible chunking — tune to your data and memory

    ds = xr.open_dataset(netcdf_path, chunks=xr_chunks)
    ds = ds.rio.write_crs("EPSG:4326", inplace=True)

    results = []
    total = shapes.shape[0]
    for _, row in tqdm(shapes.iterrows(), total=total, desc="Processing shapes"):
        shape_geometry = row.geometry
        date = row["date"]
        polygon_id = row.get("id", None) # à modifier : ID col

        # ensure date is datetime
        if isinstance(date, str):
            date = datetime.strptime(date, "%Y-%m-%d")

        # Precompute date ranges dictionary for this polygon's date
        date_ranges = {label: get_dates(date, label) for label in time_steps}

        result_entry = {"id": polygon_id} # à modifier : ID col

        for label, (start_date, end_date) in date_ranges.items():
            try:
                # Select time slice lazily
                ds_slice = ds.sel(time=slice(start_date, end_date))

                # If no times in slice -> short-circuit
                if getattr(ds_slice, "time", None) is None or ds_slice.time.size == 0:
                    result_entry[f"Cop_{variable}_{label}_mean"] = None
                    result_entry[f"Cop_{variable}_{label}_min"]  = None
                    result_entry[f"Cop_{variable}_{label}_max"]  = None
                    # cleanup
                    del ds_slice
                    gc.collect()
                    continue

                da = ds_slice[variable]

                # Compute temporal aggregates lazily (these are still dask-backed if chunks were set)
                da_mean = da.mean(dim="time", skipna=True)
                da_min  = da.min(dim="time", skipna=True)
                da_max  = da.max(dim="time", skipna=True)

                # Evaluate only the aggregated 2D arrays into memory (smaller than the full 3D)
                # Convert to numpy arrays *once* and pass to exact_extract
                # If exact_extract accepts xarray.DataArray directly, you can pass da_mean.compute()
                da_mean_np = da_mean.compute()
                da_min_np = da_min.compute()
                da_max_np = da_max.compute()

                # now compute zonal stats using the 2D arrays
                mean_val, _, _ = compute_stats_from_2d_raster(da_mean_np, shape_geometry, agg_name='mean')
                _, min_val, _ = compute_stats_from_2d_raster(da_min_np, shape_geometry, agg_name='min')
                _, _, max_val = compute_stats_from_2d_raster(da_max_np, shape_geometry, agg_name='max')

                result_entry[f"Cop_{variable}_{label}_mean"] = mean_val
                result_entry[f"Cop_{variable}_{label}_min"]  = min_val
                result_entry[f"Cop_{variable}_{label}_max"]  = max_val

            except Exception as e:
                print(f"Error for polygon {polygon_id} label {label}: {e}")
                result_entry[f"Cop_{variable}_{label}_mean"] = None
                result_entry[f"Cop_{variable}_{label}_min"]  = None
                result_entry[f"Cop_{variable}_{label}_max"]  = None

            finally:
                # ensure we release references
                for v in ("ds_slice", "da", "da_mean", "da_min", "da_max", "da_mean_np", "da_min_np", "da_max_np"):
                    if v in locals():
                        try:
                            del locals()[v]
                        except Exception:
                            pass
                gc.collect()

        results.append(result_entry)

    # close dataset and free resources
    try:
        ds.close()
    except Exception:
        pass
    gc.collect()

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)
    return results_df












#---------- EXECUTION DU CODE & PARAMETRES ----------------    
# Parameters
geojson_dir = "./grids"     
netcdf_path = "/marbec-data/BiodivMed/Copernicus/SST_MED_SST_L4_NRT_OBSERVATIONS_010_004_c_V2_SST_20230331-20230901.nc"   
output_dir = "/marbec-data/BiodivMed/output/Extraction_SST_2018-2024/grid_v1.0/"         









# Find all matching .geojson files
pattern = os.path.join(geojson_dir, "*.geojson")
geojson_files = glob.glob(pattern)


""" OPTIONAL 
# Remove already done 
geojson_files = [
    f for f in geojson_files
    if not (
        f.endswith("grid_v1.0_2023-07-01.geojson") or
        f.endswith("grid_v1.0_2023-09-01.geojson")
    )
]
"""

for geojson_path in geojson_files:
    # get date string from filename
    fname = os.path.basename(geojson_path)
    m = re.search(r"(\d{4}-\d{2}-\d{2})", fname)
    if m:
        date_str = m.group(1)
    else:
        # fallback if no date detected
        date_str = "nodate"

    # build output CSV name, including the date
    output_path = os.path.join(
        output_dir,
        f"Cop_SST_stats_{date_str}.csv" # à modifier avec variable SST/CHL
    )

    print(f"Running extraction for {fname} -> {output_path}")
    process_geojson(
        geojson_path=geojson_path,
        netcdf_path=netcdf_path,
        output_path=output_path,
        variable="analysed_sst",
        time_steps=["month"]
    )