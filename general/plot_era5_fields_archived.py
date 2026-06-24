# plot_era5_fields_archived.py
# 
# Download and plots ERA5 fields from netcdf files
#
# Steven Cavallo
# May 2026
import os
import sys, datetime, time
import xarray as xr
import urllib.request
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from scipy import ndimage

import weather_modules as wm
import utilities_modules as um
from mstats import *
####################################################
# User Options
####################################################
# Define GFS file parameters
era5_path = '/Users/scavallo/Documents/Work/Research/Students/Will_Durbin/ERA5/era5_pres_and_sfc_2026012506_2026012506.nc'
analysis_time = '2026012506'
time_index = 0
forecast_hour = '000'

# Custom options
overlay_windbarbs = True
plot_region = 'CONUS' # CONUS or NAMERICA
plot_pressure_levels = [850,500,900] #[main_level, upper_level, lower_level] 
windbarb_pressure_level = 850 # If overlay_windbarbs = True, then this is the pressure level for wind barbs that are plotted
barb_interval = 10 # interval for wind barbs; only used if overlay_windbarbs = True
num_smoothing_iterations = 5 # Number of smoothings passes to make field more synoptic-scale
title_fontsize = 14
label_fontsize = 14
contour_fontsize = 10

imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
#imagedir = '/Users/scavallo/Documents/Work/Research/Students/Will_Durbin/images/'
#figname_prefix = f"gfs_700t_700ghgt_" 
#figname_prefix = f"gfs_1000ghgt_500ghgt_"
#figname_prefix = f"gfs_850temperature_850wind_"
#figname_prefix = f"gfs_deformation_850mb_"
#figname_prefix = f"gfs_frontogenesis_850mb_"
#figname_prefix = f"gfs_1000:500thickness_"
#figname_prefix = f"gfs_850hPaTemp_"
#figname_prefix = f"gfs_700hPa_ghgt_temp_wind"
figname_prefix = f"gfs_qomegaTraditional_700mb_"
#figname_prefix = f"gfs_qomegaSutcliffeTrenberth_650mb_"
#figname_prefix = f"gfs_qomegaQvec_700mb_"
#figname_prefix = f"gfs_qheighttendency_500mb_"
#figname_prefix = f"gfs_500ghgt_slp_"

# Generally do not change the setting below
proj_latlon = [35. , -95.] # Do not change 

plot_option = 15 # 0 to print file content to screen
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
                # 15 QG Omega forcings from traditional form of QG Omega equation
                # 16 QG Omega forcings from Sutcliffe-Trenberth form of QG Omega equation
                # 17 QG Omega forcings from Q-vector form of QG Omega equation
                # 18 QG Height Tendency forcings

####################################################
# End of User options
####################################################   
figname = figname_prefix + analysis_time + '_f' + str(forecast_hour)
print(figname)

# ---- ERA5 FILE ----
dataset = xr.open_dataset(era5_path)
# If multiple times exist, select one:
#dataset = dataset.isel(time=time_index)
                    

if plot_option == 14:
    import metpy.calc as mpcalc
    from metpy.units import units

    level = plot_pressure_levels[0] * units.hPa

    T = dataset["t"].metpy.sel(vertical=level).squeeze()  # K
    u = dataset["u"].metpy.sel(vertical=level).squeeze()  # m/s
    v = dataset["v"].metpy.sel(vertical=level).squeeze()  # m/s

    # Potential temperature for frontogenesis calculation
    theta = mpcalc.potential_temperature(level, T)

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
    plot_pressure_levels = plot_pressure_levels #[700,700,500]
    #windbarb_pressure_level = 700
elif plot_option == 14:
    overlay_windbarbs = True
elif plot_option == 15:
    plot_pressure_levels = [700,500,500]
    windbarb_pressure_level = 700
elif plot_option == 16:
    plot_pressure_levels = [700,1000,500]
    windbarb_pressure_level = 700
elif plot_option == 17:
    plot_pressure_levels = [700,500,500]
    windbarb_pressure_level = 700
elif plot_option == 18:
    plot_pressure_levels = [500,500,500]
    windbarb_pressure_level = 500
else:
    print("Invalid plot_option")
    exit()

pressure_level_options = np.array(plot_pressure_levels).astype(int) 

# Extract the variables needed for plotting
lons = dataset['longitude']
lats = dataset['latitude']
if lats[0] > lats[-1]:
    dataset = dataset.sortby('latitude')
    lats = dataset['latitude']
X,Y = np.meshgrid(lons, lats) 
# Pressure levels (ERA5 uses "level" in hPa)
heights_mainlevel = dataset['z'].sel(pressure_level=pressure_level_options[0]).squeeze() / 9.81
heights_upperlevel = dataset['z'].sel(pressure_level=pressure_level_options[1]).squeeze() / 9.81
heights_lowerlevel = dataset['z'].sel(pressure_level=pressure_level_options[2]).squeeze() / 9.81

temps_mainlevel = dataset['t'].sel(pressure_level=pressure_level_options[0]).squeeze()
temps_upperlevel = dataset['t'].sel(pressure_level=pressure_level_options[1]).squeeze()
temps_lowerlevel = dataset['t'].sel(pressure_level=pressure_level_options[2]).squeeze()

u_wind = dataset['u'].sel(pressure_level=windbarb_pressure_level).squeeze()
v_wind = dataset['v'].sel(pressure_level=windbarb_pressure_level).squeeze()

t2m = dataset['t2m'].squeeze()     
slp_in = dataset['msl'].squeeze() / 100.
if( (plot_option == 5) or (plot_option == 6) ):
    del u_wind, v_wind
    u_wind = dataset['u10'].squeeze()
    v_wind = dataset['v10'].squeeze()


#absv_mainlevel = wm.vertical_vorticity_latlon(u, v, lats, lons, abs_opt)

relv_mainlevel = dataset['vo'].sel(pressure_level=pressure_level_options[0]) 
omeg_e = (2*np.pi) / (24*3600)
f = 2.*omeg_e*np.sin(np.deg2rad(Y))
absv_mainlevel = (relv_mainlevel + f)*1e5

thickness = heights_lowerlevel-heights_upperlevel

if ( (plot_option == 15) or (plot_option == 16) or (plot_option == 17) ):
    f0 = 2.*omeg_e*np.sin(np.deg2rad(45.))
    nlat_tup = np.shape(np.array(lats))
    nlon_tup = np.shape(np.array(lons))
    phPa = dataset['pressure_level']
    p = phPa*100.
    parr = np.array(p).astype('f')
    nz_tup = np.shape(parr)
    nz = int(nz_tup[0])
    nlat = int(nlat_tup[0])
    nlon = int(nlon_tup[0])
    print(nz,nlat,nlon)

    p3d = np.zeros((nz,nlat,nlon))

    for kk in range(0,nz):
        p3d[kk,:,:] = parr[kk]

    dp = np.gradient(parr)

    dp3d = np.zeros_like(p3d).astype('f')
    for kk in range(0,nz):
        dp3d[kk,:,:] = dp[kk]

    uin = np.array(dataset['u']).squeeze()
    vin = np.array(dataset['v']).squeeze()
    tin = np.array(dataset['t']).squeeze()

    for ii in range(0,5): 
        uin = ndimage.gaussian_filter(uin,3.75)
        vin = ndimage.gaussian_filter(vin,3.75)
        tin = ndimage.gaussian_filter(tin,3.75)

    forcing, vort_adv_term, temp_adv_term, sutcliffe_trenberth_term, Qvec_term, Qx, Qy, sigma = wm.compute_qg_forcing_latlon(uin, vin, tin, thickness, np.array(lats), np.array(lons), parr, f0, dp, R=287.0, a=6.371e6)
    zindex = np.where(phPa == pressure_level_options[0])

    forcing_level = forcing[zindex,:,:].squeeze()
    vortadv_term_level = vort_adv_term[zindex,:,:].squeeze()
    tempadv_term_level = temp_adv_term[zindex,:,:].squeeze()
    sutcliffe_trenberth_term_level = sutcliffe_trenberth_term[zindex,:,:].squeeze()
    Qvec_term_level = Qvec_term[zindex,:,:].squeeze()
    sigma_level = sigma[zindex,:,:].squeeze()
    Qx_level = Qx[zindex,:,:].squeeze()
    Qy_level = Qy[zindex,:,:].squeeze()

    # cut off data at the poles
    forcing_level[0,:] = float('NaN')
    forcing_level[-1,:] = float('NaN')
    vortadv_term_level[0,:] = float('NaN')
    vortadv_term_level[-1,:] = float('NaN')
    tempadv_term_level[0,:] = float('NaN')
    tempadv_term_level[-1,:] = float('NaN')
    sutcliffe_trenberth_term_level[0,:] = float('NaN')
    sutcliffe_trenberth_term_level[-1,:] = float('NaN')
    Qvec_term_level[0,:] = float('NaN')
    Qvec_term_level[-1,:] = float('NaN')
    Qx_level[0,:] = float('NaN')
    Qx_level[-1,:] = float('NaN')
    Qy_level[0,:] = float('NaN')
    Qy_level[-1,:] = float('NaN')

    
    mstats(forcing_level)
    mstats(Qvec_term_level)
    aa = Qvec_term_level - forcing_level

    mstats(Qx_level)

if ( (plot_option == 18) ):
    omeg_e = (2*np.pi) / (24*3600)
    f0 = 2.*omeg_e*np.sin(np.deg2rad(45.))
    nlat_tup = np.shape(np.array(lats))
    nlon_tup = np.shape(np.array(lons))
    phPa = dataset['pressure_level']
    p = phPa*100.
    parr = np.array(p).astype('f')
    nz_tup = np.shape(parr)
    nz = int(nz_tup[0])
    nlat = int(nlat_tup[0])
    nlon = int(nlon_tup[0])
    print(nz,nlat,nlon)

    p3d = np.zeros((nz,nlat,nlon))

    for kk in range(0,nz):
        p3d[kk,:,:] = parr[kk]

    dp = -1*np.gradient(parr)


    uin = np.array(dataset['u']).squeeze()
    vin = np.array(dataset['v']).squeeze()
    tin = np.array(dataset['t']).squeeze()

    for ii in range(0,5): 
        uin = ndimage.gaussian_filter(uin,3.75)
        vin = ndimage.gaussian_filter(vin,3.75)
        tin = ndimage.gaussian_filter(tin,3.75)

    forcing, vort_adv_term, diff_temp_adv_term, sigma = wm.compute_qg_HeightTendency_latlon(uin, vin, tin, np.array(lats), np.array(lons), parr, f0, dp, R=287.0, a=6.371e6)

    zindex = np.where(phPa == pressure_level_options[0])

    forcing_level = forcing[zindex,:,:].squeeze()
    vortadv_term_level = vort_adv_term[zindex,:,:].squeeze()
    diff_tempadv_term_level = diff_temp_adv_term[zindex,:,:].squeeze()
    sigma_level = sigma[zindex,:,:].squeeze()

    # cut off data at the poles
    forcing_level[0,:] = float('NaN')
    forcing_level[-1,:] = float('NaN')
    vortadv_term_level[0,:] = float('NaN')
    vortadv_term_level[-1,:] = float('NaN')
    diff_tempadv_term_level[0,:] = float('NaN')
    diff_tempadv_term_level[-1,:] = float('NaN')
    
    mstats(forcing_level)
    mstats(diff_tempadv_term_level)
    mstats(vortadv_term_level)

absv_mainlevel = np.array(absv_mainlevel)
slp = np.array(slp_in)
for ii in range(0,num_smoothing_iterations): 
    absv_mainlevel = ndimage.gaussian_filter(absv_mainlevel,0.75)


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

    cntrs3 = np.arange(5040,6000,60)


    contour = ax.contour(lons, lats, thickness, levels=cntrs1,
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=[5400],
                     colors='blue', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, [5400], inline=False, inline_spacing = 0.01, fontsize=contour_fontsize)
    contour = ax.contour(lons, lats, thickness, levels=cntrs2,
                     colors='red', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, cntrs2, colors='black', inline=False, inline_spacing=0.01,fontsize=contour_fontsize) 
    

    #contour = ax.contour(lons, lats, heights_mainlevel, levels=cntrs3,
    #                 colors='black', linestyles='solid', transform=ccrs.PlateCarree())
    #ax.clabel(contour, cntrs3, colors='green', inline=False, inline_spacing=0.01,fontsize=contour_fontsize)

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
                     colors='red', linestyles='solid', transform=ccrs.PlateCarree())
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

if plot_option == 15:
    overlay_windbarbs = True
    scaling_factor = 10**12
    forcing_level = forcing_level*scaling_factor
    vortadv_term_level = vortadv_term_level*scaling_factor
    tempadv_term_level = tempadv_term_level*scaling_factor

    cmin = -6
    cmax = 6
    cint = 0.2
    #forcing_level = um.filter_numeric_nans(forcing_level,cmax,cmax,'high')
    #forcing_level = um.filter_numeric_nans(forcing_level,cmin,cmin,'low')
  
    for ii in range(0,num_smoothing_iterations): 
        forcing_level = ndimage.gaussian_filter(forcing_level,0.75)
        vortadv_term_level = ndimage.gaussian_filter(vortadv_term_level,0.75)
        tempadv_term_level = ndimage.gaussian_filter(tempadv_term_level,0.75)
  

    cntrs1 = np.arange(cmin,cmax+cint,cint)
    cint2 = 3*cint
    cntrs2 = np.arange(cmin,cmax+cint2,cint2)

    cf = ax.contourf(lons, lats, forcing_level, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("x 10^12 Pa s-1 m-1")
   
    contour = ax.contour(lons, lats, tempadv_term_level, levels=cntrs2[np.where(cntrs2>0)],
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

    contour = ax.contour(lons, lats, vortadv_term_level, levels=cntrs2[np.where(cntrs2>0)],
                     colors='red', linestyles='solid',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

if plot_option == 16:
    overlay_windbarbs = True
    scaling_factor = 10**12
    
    sutcliffe_trenberth_term_level = sutcliffe_trenberth_term_level*scaling_factor

    cmin = -6
    cmax = 6
    cint = 0.2
    sutcliffe_trenberth_term_level = um.filter_numeric_nans(sutcliffe_trenberth_term_level,cmax,cmax,'high')
    sutcliffe_trenberth_term_level = um.filter_numeric_nans(sutcliffe_trenberth_term_level,cmin,cmin,'low')
  
    for ii in range(0,num_smoothing_iterations): 
       sutcliffe_trenberth_term_level = ndimage.gaussian_filter(sutcliffe_trenberth_term_level,0.75)
  

    cntrs1 = np.arange(cmin,cmax+cint,cint)
    cint2 = 3*cint
    cntrs2 = np.arange(cmin,cmax+cint2,cint2)

    cf = ax.contourf(lons, lats, sutcliffe_trenberth_term_level, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("x 10^12 Pa s-1 m-1")
   
if plot_option == 17:
    overlay_windbarbs = False
    scaling_factor = 10**12
    
    Qvec_term_level = Qvec_term_level*scaling_factor

    cmin = -6
    cmax = 6
    cint = 0.2
    Qvec_term_level = um.filter_numeric_nans(Qvec_term_level,cmax,cmax,'high')
    Qvec_term_level = um.filter_numeric_nans(Qvec_term_level,cmin,cmin,'low')

    qx_ceil = np.nanmean(Qx_level) + (2.*np.nanstd(Qx_level))
    qx_floor = np.nanmean(Qx_level) - (2.*np.nanstd(Qx_level))
    Qx_level = um.filter_numeric_nans(Qx_level,qx_ceil,qx_ceil,'high')
    Qx_level = um.filter_numeric_nans(Qx_level,qx_floor,qx_floor,'low')
    
    qy_ceil = np.nanmean(Qy_level) + (2.*np.nanstd(Qy_level))
    qy_floor = np.nanmean(Qy_level) - (2.*np.nanstd(Qy_level))   
    Qy_level = um.filter_numeric_nans(Qy_level,qy_ceil,qy_ceil,'high')
    Qy_level = um.filter_numeric_nans(Qy_level,qy_floor,qy_floor,'low')
  
    for ii in range(0,num_smoothing_iterations): 
       Qvec_term_level = ndimage.gaussian_filter(Qvec_term_level,0.75)
  

    cntrs1 = np.arange(cmin,cmax+cint,cint)
    cint2 = 3*cint
    cntrs2 = np.arange(cmin,cmax+cint2,cint2)

    cf = ax.contourf(lons, lats, Qvec_term_level, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("x 10^12 Pa s-1 m-1")

    skip = 6
    cff = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
           Qx_level[::skip, ::skip], Qy_level[::skip, ::skip],
           color='black', transform=ccrs.PlateCarree(),scale=2e-11)
if plot_option == 18:
    overlay_windbarbs = True
    scaling_factor = 10**14
    forcing_level = forcing_level*scaling_factor
    vortadv_term_level = vortadv_term_level*scaling_factor
    diff_tempadv_term_level = diff_tempadv_term_level*scaling_factor

    cmin = -16
    cmax = 16
    cint = 1
  
    for ii in range(0,num_smoothing_iterations): 
        forcing_level = ndimage.gaussian_filter(forcing_level,0.75)
        vortadv_term_level = ndimage.gaussian_filter(vortadv_term_level,0.75)
        diff_tempadv_term_level = ndimage.gaussian_filter(diff_tempadv_term_level,0.75)
  
    cntrs1 = np.arange(cmin,cmax+cint,cint)
    cint2 = 3*cint
    cntrs2 = np.arange(cmin,cmax+cint2,cint2)

    cf = ax.contourf(lons, lats, forcing_level, levels=cntrs1, cmap="RdBu_r", transform=ccrs.PlateCarree())
    # Add colorbar
    cbar = plt.colorbar(cf, orientation="horizontal", pad=0.05)
    cbar.set_label("x 10^14 s-3")
   
    contour = ax.contour(lons, lats, diff_tempadv_term_level, levels=cntrs2[np.where(cntrs2>0)],
                     colors='blue', linestyles='dashed',transform=ccrs.PlateCarree())
    ax.clabel(contour, colors='black',inline=False, inline_spacing = 0.001, fontsize=contour_fontsize)

    contour = ax.contour(lons, lats, vortadv_term_level, levels=cntrs2[np.where(cntrs2>0)],
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
    u_wind = np.array(u_wind)
    v_wind = np.array(v_wind)
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
    #ax.set_title(f"1000:500 hPa Thickness (GFS) (dashed; 60 m interval)\n 500-hPa heights (solid black; 60 m interval; labels in green)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
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
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa Total Deformation (color fill), Stretching Deformation (Dashed blue), and Shearing Deformation (Solid red) (s-1)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 14:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa Frontogenesis (color fill) (K / (100 km * 3 h)), Stretching Deformation (Dashed blue), and Shearing Deformation (Solid red)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 15:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa QG Ascent from total forcings (color fill) (*10^12 Pa s-1 m-1), temperature adv. (Dashed blue), and diff. vorticity adv. (Solid red)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 16:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa QG Ascent from Sutcliffe-Trenberth Thermal Wind Advection Term (color fill) (*10^12 Pa s-1 m-1)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 17:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa QG Ascent from Q-vector Form (color fill) (*10^12 Pa s-1 m-1)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
if plot_option == 18:
    ax.set_title(f"{str(plot_pressure_levels[0])}-hPa QG-$\chi$ RHS (color fill; reds = $\chi<0$) (*10^14 s-3), diff. TADV (Dashed blue = $\chi < 0$), and $\zeta_a$ adv. (Solid red = $\chi < 0$)\n {analysis_time_title}\n  {valid_time_title} (Forecast hour {forecast_hour})", fontsize=14)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

save_name = imagedir + figname + ".png"
plt.savefig(save_name, bbox_inches='tight')   

plt.show()
