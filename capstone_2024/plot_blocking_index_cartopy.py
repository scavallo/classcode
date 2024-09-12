# plot_blocking_index.py
# 
# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from netCDF4 import Dataset

import warnings
warnings.filterwarnings("ignore")

#####User Inputs################################################

data_file = '/Users/scavallo/Documents/Work/Research/Students/Capstone_2024/data/break_block_index_1979-2021.nc'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
figname = 'myfigure'

delta_phi = 30
start_time_file = '1900010100'
plot_start_time = '2007100100'
plot_end_time = '2007103118'
standardized_option = 0 # 0 for not standardized, 1 for standardized

map_projection = 'npstere' # 'npstere' for northern hemisphere polar stereorgraphic, 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
proj_latlon = [60. , 270.]
zoom = 'False'
label_fontsize = 16
#####Make Grid/Time Variable####################################
f = netCDF4.Dataset(data_file, 'r')
times = f.variables['time'][:]
f.close()
hrs_start = times[0].astype('float')-9

dtimes = []
dtimes[0:len(times)-1] = times[1:]-times[0:-1]
hinc = dtimes[0].astype('float')


date_start = um.advance_time(start_time_file,hrs_start)
print(date_start)  

sindnow = 0
datenow = date_start
while datenow < plot_start_time:
    sindnow+=1
    datenow = um.advance_time(datenow,hinc)  

eindnow = 0
datenow = date_start
while datenow < plot_end_time:
    eindnow+=1
    datenow = um.advance_time(datenow,hinc) 

print(sindnow,eindnow,datenow)


#####Load Data##################################################
data = netCDF4.Dataset(data_file, 'r')
lon = data.variables['longitude'][:]
lat = data.variables['latitude'][:]
time = data.variables['time'][sindnow]
if sindnow == eindnow:
    if standardized_option == 0:    
        db_index = data.variables['DB_index'][sindnow,:,:]
        ri_index = data.variables['RI_index'][sindnow,:,:]
        blocking = data.variables['blocking'][sindnow,:,:]
     
    if standardized_option == 1:
        db_index = data.variables['DB_Std'][sindnow,:,:]
        ri_index = data.variables['RI_Std'][sindnow,:,:]
        blocking = data.variables['B_Std'][sindnow,:,:]
else:
    if standardized_option == 0:    
        db_index = np.nanmean(data.variables['DB_index'][sindnow:eindnow,:,:],0)
        ri_index = np.nanmean(data.variables['RI_index'][sindnow:eindnow,:,:],0)
        blocking = np.nanmean(data.variables['blocking'][sindnow:eindnow,:,:],0)
        
    if standardized_option == 1:
        db_index = np.nanmean(data.variables['DB_Std'][sindnow:eindnow,:,:],0)
        ri_index = np.nanmean(data.variables['RI_Std'][sindnow:eindnow,:,:],0)
        blocking = np.nanmean(data.variables['B_Std'][sindnow:eindnow,:,:],0)

data.close()

X,Y = np.meshgrid(lon, lat)  

mstats(X)

mstats(db_index)
mstats(ri_index)
mstats(blocking)
   
   
plotvar = blocking

if standardized_option == 0:
    cbar_min_plot = -40
    cbar_max_plot = 40
    cint_plot = 2
    cbarlabel = 'Blocking Index'
if standardized_option == 1:
    cbar_min_plot = -2
    cbar_max_plot = 2
    cint_plot = 0.1
    cbarlabel = 'Standardized Blocking Index'


cflevs_plot =  np.arange(cbar_min_plot, cbar_max_plot, cint_plot)
cflevs_plot_cntrs = cflevs_plot[-10:]
cflevs_plot_ticks = np.arange(cbar_min_plot,cbar_max_plot,4*cint_plot)
cmap_opt = plt.cm.RdBu_r

golden = (np.sqrt(5)+1.)/2.
####################################################
# Figure 1
####################################################
proj = ccrs.NorthPolarStereo(central_longitude=-90)
fig = plt.figure(figsize=(7., 14./golden), dpi=128)  
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection=proj)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE,linewidth=1)
ax.add_feature(cfeature.STATES, linewidth=1)
gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.3, linestyle='-',xlocs=np.arange(-180,181,45), ylocs=np.arange(0,90,10))
ax.set_extent([-180, 180, proj_latlon[0], 90], crs=ccrs.PlateCarree())
um.bold_labels(ax,label_fontsize)

lonp,latp,plotvar_plot = um.cartopy_plot_prep(proj,X,Y,plotvar)

CS1 = plt.contourf(lonp,latp,plotvar_plot,cflevs_plot,cmap=cmap_opt,extend='both',zorder=1)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)   
labels2 = [item.get_text() for item in cbar.ax.get_xticklabels()]
labels = labels2[0::2]
    
cbar.set_ticks(cflevs_plot_ticks)
cbar.set_ticklabels(cflevs_plot_ticks,size=label_fontsize,fontweight='bold')
cbar.set_label(cbarlabel,size=label_fontsize,fontweight='bold')

save_name = imagedir + figname + ".png"
plt.savefig(save_name, bbox_inches='tight')   

plt.show()

