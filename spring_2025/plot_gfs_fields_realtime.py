# plot_gfs_fields_realtime.py
# 
# Download and plots real-time GFS fields
#
# Steven Cavallo
# January 2025

###################################################
# Imports
####################################################
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from datetime import datetime
import utilities_modules as um

####################################################
# User Options
####################################################
plot_option = 7 # 1 for 500 mb heights, 
                #2 for 500 and 1000 mb heights, 
                # 3 for 1000:500 mb thickness; 
                # 4 for slp and 500 hPa heights, 
                # 5 for slp and windbarbs
                # 6 for 2-m temperature and windbarbs
                # 7 for 850 mb temperature and windbarbs
overlay_windbarbs = False
windbarb_pressure_level = 850 # pressure level for wind barbs; only used if overlay_windbarbs = True; If plot_option = 5, this is overwritten with 10-m wind
barb_interval = 10 # interval for wind barbs; only used if overlay_windbarbs = True
title_fontsize = 14
label_fontsize = 14
contour_fontsize = 10
#proj_latlon = [60. , 270.]
proj_latlon = [35. , -95.]

analysis_time = '2025012812'
forecast_hour = 15
hinc = 3 # hour increment in files
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
figname_prefix = 'gfs_analysis_' 
####################################################
# Be very careful when changing below
####################################################
time_index = int(forecast_hour / hinc)
analysis_dir = analysis_time[0:8]

if forecast_hour < 10:
    figname = figname_prefix + analysis_time + '_f00' + str(forecast_hour)
elif forecast_hour >= 10 and forecast_hour < 100:
    figname = figname_prefix + analysis_time + '_f0' + str(forecast_hour)
else:
    figname = figname_prefix + analysis_time + '_f' + str(forecast_hour)

if plot_option == 7:
    windbarb_pressure_level = 850

# Step 1: Download GFS Data
# Replace this URL with the appropriate GFS data URL
base_url = "https://nomads.ncep.noaa.gov/dods/gfs_0p25/"
file_subdir = "gfs_0p25_"
analysis_hour = analysis_time[8:10]
url = base_url + 'gfs' + analysis_dir + "/" + file_subdir + analysis_hour + 'z'
#url = "https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs20250115/gfs_0p25_12z"
dataset = xr.open_dataset(url)

# Step 2: Extract 500 mb Geopotential Heights
# Select the variable "hgt" (geopotential height) and the 500 mb level
heights_500mb = dataset['hgtprs'].sel(lev=500, method="nearest").isel(time=time_index)
heights_850mb = dataset['hgtprs'].sel(lev=850, method="nearest").isel(time=time_index)
heights_1000mb = dataset['hgtprs'].sel(lev=1000, method="nearest").isel(time=time_index)
temps_850mb = dataset['tmpprs'].sel(lev=850, method="nearest").isel(time=time_index)
t2m_in = dataset['tmp2m'].isel(time=time_index)
slp_in = dataset['prmslmsl'].isel(time=time_index)/100.
if ( (plot_option == 5) or (plot_option == 6) ):
    u_wind = dataset['ugrd10m'].isel(time=time_index)  
    v_wind = dataset['vgrd10m'].isel(time=time_index)
else:
    u_wind = dataset['ugrdprs'].sel(lev=windbarb_pressure_level, method="nearest").isel(time=time_index)
    v_wind = dataset['vgrdprs'].sel(lev=windbarb_pressure_level, method="nearest").isel(time=time_index)

# Get latitude, longitude, and height data
lats = heights_500mb['lat']
lons = heights_500mb['lon']
heights500 = heights_500mb.values
heights850 = heights_850mb.values
heights1000 = heights_1000mb.values
temps850 = temps_850mb.values
t2m = t2m_in.values
slp = slp_in.values 

# Convert wind to knots
u_wind = u_wind.values*1.94384
v_wind = v_wind.values *1.94384

thickness = heights500-heights1000

# Extract analysis time for title
analysis_time = str(dataset.time.isel(time=0).values)
analysis_time_formatted = datetime.strptime(analysis_time, "%Y-%m-%dT%H:%M:%S.%f000").strftime("%H UTC %Y-%m-%d")
valid_time = str(heights_500mb['time'].values)
valid_time_formatted = datetime.strptime(valid_time, "%Y-%m-%dT%H:%M:%S.%f000").strftime("%H UTC %Y-%m-%d")

analysis_dt = datetime.strptime(str(analysis_time), "%Y-%m-%dT%H:%M:%S.%f000")
forecast_dt = datetime.strptime(str(valid_time), "%Y-%m-%dT%H:%M:%S.%f000")
forecast_hour = int((forecast_dt - analysis_dt).total_seconds() / 3600)

# Step 3: Plot 500 mb Geopotential Heights
golden = (np.sqrt(5)+1.)/2.
####################################################
# Figure 1
####################################################
#proj = ccrs.NorthPolarStereo(central_longitude=-90)
proj = ccrs.LambertConformal(central_longitude=proj_latlon[1], central_latitude=proj_latlon[0]) 
fig = plt.figure(figsize=(14., 14./golden), dpi=128)  
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection=proj)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE,linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=1)

# Contour plots
if ( (plot_option == 1) or (plot_option == 2) ):
    cntrs1 = np.arange(5040,6000,60)
    cntrs2 = np.arange(0,2000,60)    

    contour = ax.contour(lons, lats, heights500, levels=cntrs1,
                     colors='black', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, inline=False, inline_spacing=0.01,fontsize=contour_fontsize)
    

if plot_option == 2:
    contour = ax.contour(lons, lats, heights1000, levels=cntrs2,
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
    contour = ax.contour(lons, lats, heights500, levels=cntrs2,
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
    temps850_degC = temps850 - 273.15
    cntrs1 = np.arange(-80,85,5)
    contour = ax.contour(lons, lats, temps850_degC, levels=cntrs1,
                     colors='red', linestyles='solid', transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs1, colors='black', fmt = '%i', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)

gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.3, linestyle='-',xlocs=np.arange(-180,181,45), ylocs=np.arange(0,90,10))
#ax.set_extent([-180, 180, proj_latlon[0], 90], crs=ccrs.PlateCarree())
ax.set_extent([-125, -65, 25, 50], crs=ccrs.PlateCarree())
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
    ax.set_title(f"500 mb Geopotential Heights (GFS) (60 m interval)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 2:
    ax.set_title(f"500 mb Geopotential Heights (solid) (60 m interval)\n 1000 mb Geopotential Heights (dashed) (60 m interval)\n Valid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 3:
    ax.set_title(f"1000:500 mb Thickness (GFS) (30 m interval)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 4:
    ax.set_title(f"Mean sea level pressure (GFS) (black) (4 hPa interval)\n 500 mb Geopotential Heights (blue) (60 m interval)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 5:
    ax.set_title(f"Mean sea level pressure (GFS) (black) (4 hPa interval)\n 10-m wind (barbs)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 6:
    ax.set_title(f"2-m temperature (GFS) (red) (5 degC interval)\n 10-m wind (barbs)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
if plot_option == 7:
    ax.set_title(f"850 hPa temperature (GFS) (red) (5 degC interval)\n 850 hPa wind (barbs)\nValid Time: {valid_time_formatted} (Forecast hour {forecast_hour})\nAnalysis Time: {analysis_time_formatted}", fontsize=14)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

save_name = imagedir + figname + ".png"
plt.savefig(save_name, bbox_inches='tight')   

plt.show()


