#!/usr/bin/python 

# make_tpv_densities_FromNetcdf
#
# Computes tpv track densities and saves to numpy array
#
# Steven Cavallo
# December 2020

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import datetime as dt
import numpy.ma as ma


from mstats import *
import tpv_tracker_plotting_modules as tpvmod

r2d = 180./np.pi
###### File options ################################################
fTracks = '/arctic3/datasets/tpv_SH/tracks/tracks_low_horizPlusVert.nc'
fSave = '/home/scavallo/scripts/python_scripts/images/'
aSave = '/data2/scavallo/era_interim/sh_tracks/density_array_sh_60N65.npy'

gridpath = '/data2/scavallo/era_interim/track_files/wrf_arcticgrid.nc'

###### Track options ###############################################
min_lifetime_days = 2
lat_origin_range = [-30., -90.]
metric_min = -9999.
radius = 555000. #in meters

###### Filter options ##############################################
hemisphere = 'southern'
area_weight = False
make_tpv = True
genesis = False #filter for genesis
lysis = False #filter for lysis

######## Constants #################################################
hinc = 6.0
R_earth = 6371200.
pid = np.pi/180.

#########################################################
# End user options
#########################################################
latboxsize = 2
lonboxsize = 5

if hemisphere == 'northern':
    xlat = np.arange(30,90+(latboxsize/2),latboxsize)
    xlon = np.arange(0,365,lonboxsize)
else:
    xlat = np.arange(-90,-30+(latboxsize/2),latboxsize)
    xlon = np.arange(0,365,lonboxsize)


xloni,xlati = np.meshgrid(xlon,xlat)
den_arr = np.zeros([len(xlat),len(xlon)])
count_arr = np.zeros([len(xlat),len(xlon)])
hit_arr = np.zeros([len(xlat),len(xlon)])    

ntracks_min = (24./hinc)*min_lifetime_days


# Load data in
data = netCDF4.Dataset(fTracks,'r') # you will need to change the path to point to your data
trackLen = data.variables['lenTrack'][:]
Dates = data.variables['timeStamp'][:]
mons = np.array([int(str(x).replace('-','')[4:6]) for x in Dates])
years = np.array([int(str(x).replace('-','')[:4]) for x in Dates])
Dates = np.array([int(str(x).replace('-','')) for x in Dates])
trackStartDate = data.variables['iTimeStart'][:]
trackStartMon = mons[trackStartDate]
trackStartYear = years[trackStartDate]
trackStartDate = Dates[trackStartDate]

print(len(trackLen))
for x in range(0,np.int(0.01*len(trackLen))):
#for x in range(0,len(trackLen)):
    print("On point %d of %d" %(x,np.int(0.01*len(trackLen))))
    if x%10000 == 0:
        print ("On track {0}/{1}".format(x,len(trackLen)))
    if trackLen[x] < ntracks_min: # checking to make sure TPV track was longer than two days
        continue

    lat = data.variables['latExtr'][x,:]
    lon = data.variables['lonExtr'][x,:]    
    if not ma.is_masked(lat):
        per_life_in_polar = float(np.where(lat<=-60)[0].shape[0])/float(lat.shape[0]) # checking if TPV spent 60% of lifetime in Antarctic
    else:
        per_life_in_polar = float(np.where((lat.data<=-60)&(lat.mask!=True))[0].shape[0])/float(np.where((lat.mask!=True))[0].shape[0])
    if per_life_in_polar < 0.6:
        istpv = False
    else:
        istpv = True

    if (make_tpv == True):
        if istpv == True:
	    perc = 0.0
	    perc_thresh = -1.0
	elif istpv == False:
	    perc = 0.0
	    perc_thresh = 1.0
    else:
	perc = 0.0
	perc_thresh = -1.0            

    if ( (np.max(lat) <= lat_origin_range[0]) and (np.min(lat) >= lat_origin_range[1]) and (trackLen[x] >= ntracks_min) and (perc >= perc_thresh) ):
        print(perc,perc_thresh,istpv)
	if (genesis == True):
            vortlat = lat[0]
            vortlon = lon[0]
            fulldata = False 
		
            latnowin = vortlat	    
        elif (lysis == True):
            vortlat = lat[trackLen[x]]
            vortlon = lon[trackLen[x]]
            fulldata = False
		
            latnowin = vortlat 	    
        else:	 
	    npoints = ma.count(lat)	    
	    vortlat = lat[0:npoints]
	    vortlon = lon[0:npoints]
            fulldata = True 
		
            latnowin = np.nanmean(vortlat)    

	    for ii in xrange(0,len(xlat)):
		for jj in xrange(0,len(xlon)):	       
        	   vortices = 0		   
        	   if (fulldata == True):
        	      for kk in xrange(0,len(vortlat)):
        		 # If using NETCDF, subtract 360 from any value greater than 180 in order to
        		 # use a longitude scale of -180-180 degrees (possibly not needed now due to radian conversion)
        		 #if (False):
                	 if vortlon[kk] < 0:
                	     vortlon[kk] = vortlon[kk] + 360        	     
			 try:
			     distance = abs(tpvmod.rEarth*tpvmod.distance_on_unit_sphere(vortlat[kk],vortlon[kk],xlat[ii],xlon[jj]))
                	 except:
			     distance = 9999999999.

			 if distance <= radius:
	        	    vortices = vortices + 1
        	   else:
        	      #distance = abs(tpvmod.rEarth*helpers.distance_on_unit_sphere(vortlat,vortlon,xlat[ii],xlon[jj]))
		      try:		      
			  distance = abs(tpvmod.rEarth*tpvmod.distance_on_unit_sphere(vortlat,vortlon,xlat[ii],xlon[jj]))
        	      except:
			  distance = 9999999999.
		      if distance <= radius:
        		 vortices = vortices + 1
		   if area_weight == True:
		       den_arr[ii,jj] = den_arr[ii,jj] + (vortices/np.cos(latnowin*np.pi/180.))
		   else:	     
        	       den_arr[ii,jj] = den_arr[ii,jj] + vortices
		   if vortices > 0:
		       hit_arr[ii,jj] = hit_arr[ii,jj] + 1
		   count_arr[ii,jj] = count_arr[ii,jj] + 1
    
    continue


data.close()
file = open(aSave, 'w+')
np.savez(file,den_arr=den_arr,hit_arr=hit_arr,count_arr=count_arr)
