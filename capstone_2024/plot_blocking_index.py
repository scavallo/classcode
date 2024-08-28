# plot_blocking_index.py
# 
# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy import ndimage

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from netCDF4 import Dataset

from scipy.stats import ttest_ind, ttest_rel, ttest_ind_from_stats

import warnings
warnings.filterwarnings("ignore")

#####User Inputs################################################

#data_file = '/data2/scavallo/blocking/break_block_index_1979-2021.nc'
data_file = '/data2/scavallo/blocking/break_block_index_1979-2021_std.nc'

imagedir = '/home/scavallo/scripts/python_scripts/images/'
figname = 'myfigure'

delta_phi = 30
start_time_file = '1900010100'
plot_time = '1985010100'
standardized_option = 1 # 0 for not standardized, 1 for standardized

map_projection = 'npstere' # 'npstere' for northern hemisphere polar stereorgraphic, 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
proj_latlon = [60. , 270.]
zoom = 'False'
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

indnow = 0
datenow = date_start
while datenow < plot_time:
    indnow+=1
    datenow = um.advance_time(datenow,hinc)  

print(indnow,datenow)


#####Load Data##################################################
data = netCDF4.Dataset(data_file, 'r')
time = data.variables['time'][indnow]
lon = data.variables['longitude'][:]
lat = data.variables['latitude'][:]
if standardized_option == 0:    
    db_index = data.variables['DB_index'][indnow,:,:]
    ri_index = data.variables['RI_index'][indnow,:,:]
    blocking = data.variables['blocking'][indnow,:,:]
    time = data.variables['time'][indnow]
    lon = data.variables['longitude'][:]
    lat = data.variables['latitude'][:]
    data.close()

if standardized_option == 1:
    db_index = data.variables['DB_Std'][indnow,:,:]
    ri_index = data.variables['RI_Std'][indnow,:,:]
    blocking = data.variables['B_Std'][indnow,:,:]


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
    cbar_min_plot = -3
    cbar_max_plot = 3
    cint_plot = 0.2
    cbarlabel = 'Standardized Blocking Index'


cflevs_plot =  np.arange(cbar_min_plot, cbar_max_plot, cint_plot)
cflevs_plot_cntrs = cflevs_plot[-10:]
cflevs_plot_ticks = np.arange(cbar_min_plot,cbar_max_plot,4*cint_plot)
#cmap_opt = plt.cm.jet
cmap_opt = plt.cm.RdBu
   
golden = (np.sqrt(5)+1.)/2.     
fig = plt.figure(figsize=(12.,12.), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
if map_projection == 'ortho':
    if zoom == 'False':   
        m = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
              resolution = 'l', area_thresh = 1000.,ax=ax1)
    else:
        m1 = Basemap(projection='ortho',lon_0=proj_latlon[1],lat_0=proj_latlon[0],resolution=None)

        width = m1.urcrnrx - m1.llcrnrx
        height = m1.urcrnry - m1.llcrnry
    
        width = width*coef
        height = height*coef
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        m = Basemap(projection='ortho',lon_0=proj_latlon[1],lat_0=proj_latlon[0],resolution='l',llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)      
          
elif map_projection == 'lcc':
    m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
               rsphere=(6378137.00,6356752.3142),\
               resolution='l',area_thresh=1000.,projection='lcc',\
               lat_1=50.,lon_0=-107.,ax=ax1)           
elif map_projection == 'npstere':
    if zoom == 'false':
        m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
    else:
        m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
elif map_projection == 'spstere':
    if zoom == 'false':
        m = Basemap(projection='spstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
    else:
        m = Basemap(projection='spstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')

# draw countries, states, and differentiate land from water areas.
m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
m.drawcountries(linewidth=1, color='#444444', zorder=5)
m.drawstates(linewidth=0.66, color='#444444', zorder=4)
m.drawmapboundary
#m.fillcontinents(color='Wheat',lake_color='lightblue', zorder=1) 

# draw lat/lon grid lines every 30 degrees.
#m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 10))
m.drawmeridians(np.arange(0, 360, 30),labels=[False,False,False,False])

x, y = m(X, Y)

CS1 = m.contourf(x,y,plotvar,cmap=cmap_opt,levels=cflevs_plot, extend='both',zorder=1)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
#clabs = ['%i' % f for f in cflevs_trth]
   
labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
cbar.ax.set_xticklabels(labels, size=20)
cbar.set_label(cbarlabel,size=20)

save_name = imagedir + figname + ".png"
plt.savefig(save_name, bbox_inches='tight')   

plt.show()

