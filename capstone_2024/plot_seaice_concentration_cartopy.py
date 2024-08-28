##################################################################
# plot_seaice_concentration_cartopy                        
# Created: 07/31/2024                                              
#                 
# Sea ice concentration data is from:
# https://noaadata.apps.nsidc.org/NOAA/G02202_V4/north/aggregate/ 
#                                              
# Steven Cavallo                                                  
##################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import ticker
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4
from scipy import ndimage
from netCDF4 import Dataset
from scipy.stats import gaussian_kde
import utilities_modules as um
from mstats import *

import warnings
warnings.filterwarnings("ignore")
##################################################################
# User Settings
##################################################################
datain = '/Users/scavallo/Documents/Work/Research/Students/Capstone_2024/data/seaice_conc_monthly_nh_197811_202312_v04r00.nc'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
save = True

nyear_anomaly_plot = 2023
proj_latlon = [60.,270.]
label_fontsize = 16
##################################################################
# End user Settings
##################################################################
nyears_data = np.arange(1979,2024,1)

index_plot = nyear_anomaly_plot - nyears_data[0] 
print(index_plot)


f = netCDF4.Dataset(datain,'r')
xloni = f.variables['longitude'][:]
xlati = f.variables['latitude'][:]          
conc = f.variables['cdr_seaice_conc_monthly'][11::12,:,:].squeeze() # Get October's only
f.close()

# Mean from 1979-1996
mean_conc_early = np.nanmean(conc[1:17,:,:],0)
mean_conc_late = np.nanmean(conc[-17:,:,:],0)

mstats(mean_conc_early)
mstats(mean_conc_late)

# concentration from the year set in nyear_anomaly_plot
conc_now = conc[index_plot,:,:].squeeze()

# Anomaly
conc_anom = conc_now - mean_conc_early
conc_anom_late  = mean_conc_late - mean_conc_early


min_val = 0
max_val = 1
cint_val = 0.05             
#cmap_opt = plt.cm.RdBu_r
cmap_opt = plt.cm.Blues_r
cbarlabel = 'Concentration (fraction)'


cflevs =  np.arange(min_val, max_val+(cint_val/2), cint_val)
cflevs_ticks = cflevs[0::2]
cflevs_ticks = np.array(cflevs_ticks).astype('i')

min_val_anom = -1
max_val_anom = 1
cint_val_anom = 0.05             
cmap_opt_anom = plt.cm.RdBu_r
cbarlabel_anom = 'Concentration anomaly (fraction)'


cflevs_anom =  np.arange(min_val_anom, max_val+(cint_val_anom/2), cint_val_anom)
cflevs_anom_ticks = cflevs_anom[0::2]
cflevs_anom_ticks = np.array(cflevs_anom_ticks).astype('i')



####################################################
# Figure 1
####################################################
plotfield = mean_conc_early

proj = ccrs.NorthPolarStereo(central_longitude=-90)
fig = plt.figure(figsize=(8.5,11))
ax = plt.subplot(111, projection=proj)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE,linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.3)
gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='-')
#ax.set_extent([-145, 75, 50, 80], crs=ccrs.PlateCarree())
ax.set_extent([-180, 180, proj_latlon[0], 90], crs=ccrs.PlateCarree())
um.bold_labels(ax,label_fontsize)

out = proj.transform_points(ccrs.PlateCarree(),xloni,xlati,plotfield)
x = out[...,0]
y = out[...,1]
z = out[...,2]
z = um.filter_numeric_nans(z,min_val,float('NaN'),'low') #values below lower density are not plotted
CS1 = plt.contourf(x,y,z,cflevs,cmap=cmap_opt,extend='both',zorder=1)
cbar = plt.colorbar(CS1,shrink=0.75,orientation='horizontal',pad=0.1)
labels2 = [item.get_text() for item in cbar.ax.get_xticklabels()]
labels = labels2[0::2]

cbar.set_ticks(cflevs_ticks)
cbar.set_ticklabels(cflevs_ticks,size=label_fontsize,fontweight='bold')
cbar.set_label(cbarlabel,size=label_fontsize,fontweight='bold')

if (save == True):
   imagename = 'mean_conc_early'
   plt.savefig(imagedir + imagename, bbox_inches='tight')


####################################################
# Figure 2
####################################################
plotfield = conc_anom

proj = ccrs.NorthPolarStereo(central_longitude=-90)
fig = plt.figure(figsize=(8.5,11))
ax = plt.subplot(111, projection=proj)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE,linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.3)
gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='-')
#ax.set_extent([-145, 75, 50, 80], crs=ccrs.PlateCarree())
ax.set_extent([-180, 180, proj_latlon[0], 90], crs=ccrs.PlateCarree())
um.bold_labels(ax,label_fontsize)

out = proj.transform_points(ccrs.PlateCarree(),xloni,xlati,plotfield)
x = out[...,0]
y = out[...,1]
z = out[...,2]
z = um.filter_numeric_nans(z,min_val_anom,float('NaN'),'low') #values lower are not plotted
z = um.filter_numeric_nans(z,max_val_anom,float('NaN'),'high') #values higher are not plotted
CS1 = plt.contourf(x,y,z,cflevs_anom,cmap=cmap_opt_anom,extend='both',zorder=1)
cbar = plt.colorbar(CS1,shrink=0.75,orientation='horizontal',pad=0.1)
labels2 = [item.get_text() for item in cbar.ax.get_xticklabels()]
labels = labels2[0::2]

cbar.set_ticks(cflevs_anom_ticks)
cbar.set_ticklabels(cflevs_anom_ticks,size=label_fontsize,fontweight='bold')
cbar.set_label(cbarlabel_anom,size=label_fontsize,fontweight='bold')


if (save == True):
   imagename = 'conc_anom'
   plt.savefig(imagedir + imagename, bbox_inches='tight')


####################################################
# Figure 3
####################################################
plotfield = conc_anom_late

proj = ccrs.NorthPolarStereo(central_longitude=-90)
fig = plt.figure(figsize=(8.5,11))
ax = plt.subplot(111, projection=proj)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE,linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.3)
gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='-')
#ax.set_extent([-145, 75, 50, 80], crs=ccrs.PlateCarree())
ax.set_extent([-180, 180, proj_latlon[0], 90], crs=ccrs.PlateCarree())
um.bold_labels(ax,label_fontsize)

out = proj.transform_points(ccrs.PlateCarree(),xloni,xlati,plotfield)
x = out[...,0]
y = out[...,1]
z = out[...,2]
z = um.filter_numeric_nans(z,min_val_anom,float('NaN'),'low') #values lower are not plotted
z = um.filter_numeric_nans(z,max_val_anom,float('NaN'),'high') #values higher are not plotted
CS1 = plt.contourf(x,y,z,cflevs_anom,cmap=cmap_opt_anom,extend='both',zorder=1)
cbar = plt.colorbar(CS1,shrink=0.75,orientation='horizontal',pad=0.1)
labels2 = [item.get_text() for item in cbar.ax.get_xticklabels()]
labels = labels2[0::2]

cbar.set_ticks(cflevs_anom_ticks)
cbar.set_ticklabels(cflevs_anom_ticks,size=label_fontsize,fontweight='bold')
cbar.set_label(cbarlabel_anom,size=label_fontsize,fontweight='bold')

if (save == True):
   imagename = 'conc_anom_late'
   plt.savefig(imagedir + imagename, bbox_inches='tight')


plt.show()
