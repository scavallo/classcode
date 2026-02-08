# plot_gfs_fields_archived.py
# 
# Download and plots GFS fields from saved NCEP grb files
#
# Steven Cavallo
# March 2025
import os
import sys, datetime, time
import xarray as xr
import pygrib
#import requests
#from tqdm import tqdm
import urllib.request
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
#from datetime import datetime
from scipy import ndimage


import weather_modules as wm
import utilities_modules as um
from mstats import *

# Note: You may need to install pygrib and/or cfgrib.  If you do, use the following lines below to install them.
# pip install pygrib
# pip install xarray cfgrib

####################################################
# User Options
####################################################
# Define GFS file parameters
base_url = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
cycle = "12"  # GFS cycle time (00, 06, 12, 18 UTC)
forecast_hour = "000"  # Forecast hour (000, 003, 006, etc.)
resolution = "0p25"  # Choose from 0p25, 0p50, 1p00 (higher resolution = larger file)
date = "20260206"  # YYYYMMDD format (update as needed)

# Custom options
overlay_windbarbs = False
plot_region = 'NAMERICA' # CONUS or NAMERICA
plot_pressure_levels = [700,500,900] #[main_level, upper_level, lower_level] 
windbarb_pressure_level = 700 # If overlay_windbarbs = True, then this is the pressure level for wind barbs that are plotted
barb_interval = 10 # interval for wind barbs; only used if overlay_windbarbs = True
num_smoothing_iterations = 5 # Number of smoothings passes to make field more synoptic-scale
title_fontsize = 14
label_fontsize = 14
contour_fontsize = 10

imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
save_dir = '/Users/scavallo/Documents/Work/Classes/Synoptic-Dynamics/spring_2026/homework/data/'
#figname_prefix = f"gfs_700t_700ghgt_" 
#figname_prefix = f"gfs_1000ghgt_500ghgt_"
#figname_prefix = f"gfs_850temperature_850wind_"
figname_prefix = f"gfs_deformation_700mb_"

# Generally do not change the settings below
proj_latlon = [35. , -95.] # Do not change 
analysis_time = date+cycle
hinc = 3 # hour increment in files

plot_option = 14 # 0 to print file content to screen
                # 1 for 500 mb heights, 
                # 2 for 500 and 1000 mb heights, 
                # 3 for 1000:500 mb thickness; 
                # 4 for slp and 500 hPa heights, 
                # 5 for slp and windbarbs
                # 6 for 2-m temperature contours and windbarbs
                # 7 for plot_pressure_levels[0] temperature contours and windbarbs
                # 8 for plot_pressure_levels[0] height contours and windbarbs
                # 9 for plot_pressure_levels[0] mb vorticity fill and windbarbs
                # 10 for 700 mb vorticity fill and 1000:500 mb thickness
                # 11 for plot_pressure_levels[0] temperature fill and wind barbs
                # 12 for 700 mb temperature and geopotential height contours
                # 13 for Deformation
                # 14 for Frontogenesis (needs Metpy to be installed!)

####################################################
# End of User options
####################################################   

# Find the time index
time_index = int(int(forecast_hour) / hinc)
figname = figname_prefix + analysis_time + '_f' + str(forecast_hour)

print(figname)

# Construct file path
gfs_folder = f"gfs.{date}/{cycle}/atmos"
filename = f"gfs.t{cycle}z.pgrb2.{resolution}.f{forecast_hour}"
file_url = f"{base_url}/{gfs_folder}/{filename}"

# Save location
save_path = os.path.join(save_dir, filename)

url = f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_{resolution}.pl?dir=%2Fgfs.{date}%2F{cycle}%2Fatmos&file=gfs.t{cycle}z.pgrb2.{resolution}.f{forecast_hour}&all_var=on&all_lev=onp"  # Replace with your target URL
print(url)
output_file = f"gfs_{date}{cycle}_f{forecast_hour}.grb"  # Desired local file name
output_path = os.path.join(save_dir, output_file)   
description_file = f"gfs_{date}{cycle}_f{forecast_hour}.idx"
description_path = os.path.join(save_dir, description_file)  


if os.path.exists(output_path):
    print(f"{output_path} exists. Skipping download.")
else:
    try:
        # Get the data file
        urllib.request.urlretrieve(url, output_path)
        print(f"File successfully downloaded as {output_path}")
        # Get the description file
        urllib.request.urlretrieve(url, description_path)
        print(f"File successfully downloaded as {description_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


dataset = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
#dataset = cfgrib.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
dataset_slp = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'meanSea'})
dataset_sfc = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})
dataset_trop = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'tropopause'})
dataset_refl = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'atmosphere'})
dataset_rh = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'atmosphereSingleLayer'})

                    
if plot_option == 0:
    dataset = xr.open_dataset(output_path,engine="cfgrib")
    print(dataset)
    print("Available variables in the dataset:")
    for var in dataset.variables:
        print(f"- {var}: {dataset[var].attrs.get('long_name', 'No description')}")

    exit()

if plot_option == 14:
    import metpy.calc as mpcalc
    from metpy.units import units

    # --- Load and parse CF metadata (important for coords/units) ---
    dataset = xr.open_dataset(output_path,engine="cfgrib",filter_by_keys={'typeOfLevel': 'isobaricInhPa'}).metpy.parse_cf()
    #print("Available variables in the dataset:")
    #for var in dataset.variables:
    #    print(f"- {var}: {dataset[var].attrs.get('long_name', 'No description')}")
   
    # --- Example: choose a pressure level and time (adapt variable names to your dataset) ---
    level = plot_pressure_levels[0] * units.hPa

    T = dataset["t"].metpy.sel(vertical=level).squeeze()  # K
    u = dataset["u"].metpy.sel(vertical=level).squeeze()  # m/s
    v = dataset["v"].metpy.sel(vertical=level).squeeze()  # m/s

    # Potential temperature for frontogenesis calculation
    theta = mpcalc.potential_temperature(level, T)
    mstats(theta)

    # If your data are lat/lon, compute grid deltas (meters) from 1D lon/lat arrays:
    lons = T.metpy.x  # often "lon"
    lats = T.metpy.y  # often "lat"
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

    # MetPy expects dx,dy compatible with the data grid; broadcast if needed
    # (Common case: theta/u/v are (..., y, x))
    # fronto units: K / (m * s)
    fronto = mpcalc.frontogenesis(theta, u, v, dx=dx, dy=dy, x_dim=-1, y_dim=-2)
    # Convert to standard units for plotting
    #fronto_unit = fronto.to('K / 100 km / 3 h') # .to doesn't work in my version for some reason
    fronto_unit = fronto*1000*3600*300
    fronto_unit = np.array(fronto_unit)
    for ii in range(0,num_smoothing_iterations): 
        fronto_unit = ndimage.gaussian_filter(fronto_unit,0.75)
    mstats(fronto_unit)


if plot_option == 1:
    plot_pressure_levels = [500,500,500]
    windbarb_pressure_level = 500
    overlay_windbarbs = True
elif plot_option == 2:
    plot_pressure_levels = [500,500,1000]
    windbarb_pressure_level = 500
    overlay_windbarbs = False
elif plot_option == 3:
    plot_pressure_levels = [500,1000,500]
    windbarb_pressure_level = 500
elif plot_option == 4:
    plot_pressure_levels = [500,1000,500]
    windbarb_pressure_level = 500
elif plot_option == 5:
    aa = 0
elif plot_option == 6:  
    aa = 0
elif plot_option == 7:
    #plot_pressure_levels = [850,1000,500]
    #windbarb_pressure_level = 850
    overlay_windbarbs = True
    windbarb_pressure_level = plot_pressure_levels[0]
elif plot_option == 8:
    overlay_windbarbs = True
    windbarb_pressure_level = plot_pressure_levels[0]
elif plot_option == 9:
    #plot_pressure_levels = [500,1000,500]
    #windbarb_pressure_level = 500
    windbarb_pressure_level = plot_pressure_levels[0]
elif plot_option == 10:
    #plot_pressure_levels = [500,1000,500]   
    #windbarb_pressure_level = 500 
    plot_pressure_levels = [700,1000,500]   
    windbarb_pressure_level = 700 
elif plot_option == 11:
    #plot_pressure_levels = [700,1000,500]
    #windbarb_pressure_level = 700
    overlay_windbarbs = True
    windbarb_pressure_level = plot_pressure_levels[0]
elif plot_option == 12:
    plot_pressure_levels = [700,700,500]
    windbarb_pressure_level = 700
elif plot_option == 13:
    plot_pressure_levels = [700,700,500]
    windbarb_pressure_level = 700
elif plot_option == 14:
    overlay_windbarbs = True
else:
    print("Invalid plot_option")
    exit()

pressure_level_options = np.array(plot_pressure_levels).astype(int) 

# Extract the variables needed for plotting
lons = dataset['longitude']
lats = dataset['latitude']
heights_mainlevel = dataset['gh'].sel(isobaricInhPa=pressure_level_options[0])
heights_upperlevel = dataset['gh'].sel(isobaricInhPa=pressure_level_options[1])
heights_lowerlevel = dataset['gh'].sel(isobaricInhPa=pressure_level_options[2])
temps_mainlevel = dataset['t'].sel(isobaricInhPa=pressure_level_options[0])
temps_upperlevel = dataset['t'].sel(isobaricInhPa=pressure_level_options[1])
temps_lowerlevel = dataset['t'].sel(isobaricInhPa=pressure_level_options[2])
absv_mainlevel = dataset['absv'].sel(isobaricInhPa=pressure_level_options[0]) * 1e5
if ( (plot_option == 5) or (plot_option == 6) ):
    grbs = pygrib.open(output_path)
    u_wind = grbs.select(name='10 metre U wind component')[0].values
    v_wind = grbs.select(name='10 metre V wind component')[0].values
    t2m = grbs.select(name="2 metre temperature")[0].values
else:
    u_wind = dataset['u'].sel(isobaricInhPa=windbarb_pressure_level, method="nearest").values
    v_wind = dataset['v'].sel(isobaricInhPa=windbarb_pressure_level, method="nearest").values

slp_in = dataset_slp['prmsl'] / 100.


absv_mainlevel = np.array(absv_mainlevel)
slp = np.array(slp_in)
for ii in range(0,num_smoothing_iterations): 
    absv_mainlevel = ndimage.gaussian_filter(absv_mainlevel,0.75)
mstats(absv_mainlevel)

X,Y = np.meshgrid(lons, lats) 

if ( (plot_option == 13) or (plot_option == 14) ):
    results = wm.compute_deformation(u_wind, v_wind, lats, lons)
    stretching_deformation = np.array(results['A'])*10**5
    shearing_deformation = np.array(results['B'])*10**5
    total_deformation = np.array(results['D'])*10**5
    for ii in range(0,num_smoothing_iterations): 
        stretching_deformation = ndimage.gaussian_filter(stretching_deformation,0.75)
        shearing_deformation = ndimage.gaussian_filter(shearing_deformation,0.75)
        total_deformation = ndimage.gaussian_filter(total_deformation,0.75)



# Convert wind to knots
u_wind = u_wind*1.94384
v_wind = v_wind*1.94384

thickness = heights_lowerlevel-heights_upperlevel

####################################################
# Get the title ready
####################################################
print(analysis_time)
dt = datetime.datetime.strptime(analysis_time, '%Y%m%d%H')
analysis_time_title = 'Analysis time: %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))

valid_time = um.advance_time(analysis_time,int(forecast_hour))

dtt = datetime.datetime.strptime(valid_time, '%Y%m%d%H')
valid_time_title = 'Valid time: %s at %s UTC' % (dtt.strftime('%d %b %Y'), dtt.strftime('%H00'))  


golden = (np.sqrt(5)+1.)/2.
####################################################
# Figure 1
####################################################
#proj = ccrs.NorthPolarStereo(central_longitude=-90)
proj = ccrs.LambertConformal(central_longitude=proj_latlon[1], central_latitude=proj_latlon[0]) 
fig = plt.figure(figsize=(14., 14./golden), dpi=128)  
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection=proj)
ax.add_feature(cfeature.LAND, facecolor='beige', edgecolor='black')
ax.add_feature(cfeature.OCEAN, facecolor='skyblue')
ax.add_feature(cfeature.COASTLINE,linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=1)


# Contour plots
if ( (plot_option == 1) or (plot_option == 2) ):
    cntrs1 = np.arange(5040,6000,60)
    cntrs2 = np.arange(0,2000,60)    

    contour = ax.contour(lons, lats, heights_mainlevel, levels=cntrs1,
                     colors='black', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
    

if plot_option == 2:
    contour = ax.contour(lons, lats, heights_lowerlevel, levels=cntrs2,
                     colors='black', linestyles='dashed', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs2, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)

if plot_option == 3:
    cntrs1 = np.arange(4020,5400,30)
    cntrs2 = np.arange(5430,5800,30)    


    contour = ax.contour(lons, lats, thickness, levels=cntrs1,
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=[5400],
                     colors='blue', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, [5400], inline=False, inline_spacing = 0.01, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=cntrs2,
                     colors='red', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs2, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize) 

if plot_option == 4:
    overlay_windbarbs = False
    cntrs1 = np.arange(892,1060,4)
    cntrs2 = np.arange(5040,6000,60)
    contour = ax.contour(lons, lats, slp, levels=cntrs1,
                     colors='black', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, heights_mainlevel, levels=cntrs2,
                     colors='blue', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs2, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)

if plot_option == 5:
    overlay_windbarbs = True
    cntrs1 = np.arange(892,1060,4)
    contour = ax.contour(lons, lats, slp, levels=cntrs1,
                     colors='black', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
if plot_option == 6:
    overlay_windbarbs = True
    t2m_degC = t2m - 273.15
    cntrs1 = np.arange(-60,65,5)
    contour = ax.contour(lons, lats, t2m_degC, levels=cntrs1,
                     colors='red', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, colors='black', fmt = '%i', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
if plot_option == 7:
    overlay_windbarbs = True
    temps850_degC = temps_mainlevel - 273.15
    cntrs1 = np.arange(-80,85,2)
    contour = ax.contour(lons, lats, temps850_degC, levels=cntrs1,
                     colors='black', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, colors='black', fmt = '%i', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
if plot_option == 8:
    overlay_windbarbs = True
    cntrs1 = np.arange(1000,2000,50)
    cntrs2 = np.arange(0,2000,60)    

    contour = ax.contour(lons, lats, heights_mainlevel, levels=cntrs1, colors='black', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
if plot_option == 9:
    overlay_windbarbs = True
   
    cf = ax.contourf(lons, lats, absv_mainlevel, levels=np.arange(5, 50, 5), cmap="gist_heat_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("500 mb Vorticity (1/s)")
if plot_option == 10:
    overlay_windbarbs = True
   
    cf = ax.contourf(lons, lats, absv_mainlevel, levels=np.arange(5, 50, 5), cmap="gist_heat_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("700 mb Vorticity (1/s)")

    cntrs1 = np.arange(4020,5400,30)
    cntrs2 = np.arange(5430,5800,30)    

    contour = ax.contour(lons, lats, thickness, levels=cntrs1,
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=[5400],
                     colors='blue', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, [5400], inline=False, inline_spacing = 0.01, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=cntrs2,
                     colors='red', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs2, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize) 
if plot_option == 11:
    overlay_windbarbs = True
    temps700_degC = temps_mainlevel - 273.15
    cntrs1 = np.arange(-20,20,1)
    
    cf = ax.contourf(lons, lats, temps700_degC, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label(r"700 hPa temperature ($^{\circ}$C)")
if plot_option == 12:
    overlay_windbarbs = True
   
    cntrs1 = np.arange(2700,3330,30)

    cf = ax.contourf(lons, lats, heights_mainlevel, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("700 mb heights (m)")

    cntrs2 = np.arange(-60,40,2)   
    temps700_degC = temps_mainlevel - 273.15 

    contour = ax.contour(lons, lats, temps700_degC, levels=cntrs2,
                     colors='red', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)
if plot_option == 13:
    overlay_windbarbs = True   
    um.filter_numeric_nans(total_deformation,100,100,'high')

    cntrs1 = np.arange(10,71,5)

    cf = ax.contourf(lons, lats, total_deformation, levels=cntrs1, cmap="viridis", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("Deformation (1/s)")
   
    contour = ax.contour(lons, lats, stretching_deformation, levels=cntrs1,
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

    contour = ax.contour(lons, lats, shearing_deformation, levels=cntrs1,
                     colors='red', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)
if plot_option == 14:
    overlay_windbarbs = True   
    cmax = 12
    um.filter_numeric_nans(fronto_unit,cmax,cmax,'high')


    cntrs1 = np.arange(0.5,cmax+0.1,0.25)
    cntrs2 = np.arange(10,71,5)

    cf = ax.contourf(lons, lats, fronto_unit, levels=cntrs1, cmap="hot_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("K / (100 km * 3 h)")
   
    contour = ax.contour(lons, lats, stretching_deformation, levels=cntrs2,
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

    contour = ax.contour(lons, lats, shearing_deformation, levels=cntrs2,
                     colors='red', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.3, linestyle='-',xlocs=np.arange(-180,181,45), ylocs=np.arange(0,90,10))
if plot_region == 'CONUS':
    ax.set_extent([-125, -65, 25, 50], crs=ccrs.PlateCarree()) # CONUS
elif plot_region == 'NAMERICA':
    ax.set_extent([-145, -35, 25, 70], crs=ccrs.PlateCarree()) # North America
um.bold_labels(ax,label_fontsize)

if overlay_windbarbs == True:
# Add wind barbs
# Subsample data for better visualization (e.g., every 5th grid point)
    ax.barbs(
        lons[::barb_interval], lats[::barb_interval],
        u_wind[::barb_interval, ::barb_interval],
        v_wind[::barb_interval, ::barb_interval],
        length=6, transform=ccrs.PlateCarree()
)
    
# Add title and labels
if plot_option == 1:
    ax.set_title(f"500-hPa Geopotential Heights (GFS) (60 m interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 2:
    ax.set_title(f"500-hPa Geopotential Heights (solid) (60 m interval)\n 1000-hPa Geopotential Heights (dashed) (60 m interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 3:
    ax.set_title(f"1000:500 mb Thickness (GFS) (30 m interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 4:
    ax.set_title(f"Mean sea level pressure (GFS) (black) (4 hPa interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 5:
    ax.set_title(f"Mean sea level pressure (GFS) (black) (4 hPa interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 6:
    ax.set_title(f"2-m temperature (GFS) (red) (5 degC interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 7:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa temperature (GFS) (contours) (2 degC interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 8:
    ax.set_title(f"850-hPa Geopotential Heights (solid) (50 m interval) and wind barbs\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 9:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa Vorticity (GFS) (*10-5 s-1) and Wind (knots)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 10:
    ax.set_title(f"700-hPa Vorticity (GFS) (*10-5 s-1), Wind (knots), and 1000:500 hPa thickness (30 m interval)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 11:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa temperature (GFS) (5 degC interval) and Wind (knots)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 12:
    ax.set_title(f"700-hPa temperature (GFS) (2 degC interval), Heights (30 m interval), and Wind (knots)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 13:
    ax.set_title(f"700-hPa Total Deformation (color fill), Stretching Deformation (Dashed blue), and Shearing Deformation (Solid red) (s-1)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 14:
    ax.set_title(f"700-hPa Frontogenesis (color fill) (K / (100 km * 3 h)), Stretching Deformation (Dashed blue), and Shearing Deformation (Solid red) (s-1)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

save_name = imagedir + figname + ".png"
plt.savefig(save_name, bbox_inches='tight')   

plt.show()
