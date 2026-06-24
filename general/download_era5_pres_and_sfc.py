#!/usr/bin/env python
# download_era5_pres_and_sfc.py
import cdsapi
import xarray as xr

c = cdsapi.Client()

# ---------------------------
# 1. Download pressure-level data
# ---------------------------
pres_file = "era5_pres.nc"
sfc_file = "era5_sfc.nc"
combined_file = 'era5_pres_and_sfc_2026012418_2026012418.nc'

c.retrieve(
    "reanalysis-era5-pressure-levels",
    {
        "product_type": ["reanalysis"],
        "variable": [
            "geopotential",
            "ozone_mass_mixing_ratio",
            "potential_vorticity",
            "relative_humidity",
            "specific_humidity",
            "temperature",
            "u_component_of_wind",
            "v_component_of_wind",
            "vertical_velocity",
            "vorticity",
        ],
        "year": ["2026"],
        "month": ["01"],
        "day": ["24"],
        "time": ["18:00"],
        "pressure_level": [
            "1","5","10","20","50","70","100","125","150","175",
            "200","225","250","300","350","400","450","500","550",
            "600","650","700","750","775","800","825","850","875",
            "900","925","950","975","1000"
        ],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [90, -180, 0, 180],
    },
    pres_file,
)

# ---------------------------
# 2. Download single-level data
# ---------------------------

c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": ["reanalysis"],
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_temperature",
            "mean_sea_level_pressure",
			"vertical_integral_of_northward_water_vapour_flux",
			"total_column_water_vapour"
        ],
        "year": ["2026"],
        "month": ["01"],
        "day": ["24"],
        "time": ["18:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [90, -180, 0, 180],
    },
    sfc_file,
)

# ---------------------------
# 3. Merge into one file
# ---------------------------
ds_pres = xr.open_dataset(pres_file)
ds_sfc = xr.open_dataset(sfc_file)

# Merge datasets
ds_merged = xr.merge([ds_pres, ds_sfc])

# Save final combined file
ds_merged.to_netcdf(combined_file)

print("Combined file")