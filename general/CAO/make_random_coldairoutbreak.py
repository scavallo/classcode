#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

import tpv_tracker_plotting_modules as tpvmod

from scipy import ndimage
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

###################################
# Set user options
###################################
#event_fdir = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/'
#eventfile_in_full = 'CAOinfo.txt'
#eventfile_out_random = 'test_CAO_60Ngen_randomTPV_with95confidence.txt'
event_fdir = '/data1/scavallo/data/seaice/VRILEs/1979_2019/'
eventfile_in_full = 'VRILE_masks_jja_bwfilter_3d_5percentile_dtUnique01_1979_2019.dat'
eventfile_out_random = 'VRILE_60N65_randomTPV_with95confidence.txt'

track_fdir = '/data2/scavallo/era_interim/track_files/'
#track_file = 'tpvs_60Ngen_djf_CAOboxfull_1979010100_2018123118.dat'
track_file = 'trth_tracks_60N65_jja_1979010100_2018123118.dat'

date_format = 'yyyymmddhh'
numiters = 10000
percentiles = [2.5,97.5]
#confidence_interval_limits = [5,95]

random_year_range = [1979,2018]
random_day_range = [1,31]

month_choices = [6,7,8]
hour_choices = [0,6,12,18]
###################################
# END user options
###################################

event_outfile = open(event_fdir + eventfile_out_random,'w')
event_outfile.write('Nposs N100 N250 N500 N750 N1000 N1250 N1500 N1750 N2000 N2500 N3000 N3500 N4000 N4500 N5000 N6000')

# Read in event dates
bb = np.loadtxt(event_fdir + eventfile_in_full, skiprows=1)       
dateevents = bb[:,0]
event_lats = bb[:,1]
event_lons = bb[:,2]
event_mags = bb[:,3]
Nposs = np.int(len(dateevents)-2) # Subtract 2 because of 2019   

# Read in full track file
aa = np.loadtxt(track_fdir + track_file, skiprows=1)       
datelist = aa[:,0]
lat = aa[:,1]	
lon = aa[:,2]
thetamin = aa[:,3]	    
thetaamp = aa[:,4]
vort_circ = aa[:,5]	    
vort_radius = aa[:,6]

nrows, ncols = np.shape(aa)
print nrows, ncols

[inds] = np.where(datelist > 190000000)
dates_tpvs = datelist[inds]
lats_tpvs = lat[inds]
lons_tpvs = lon[inds]

dates_tpvs_nohour = []
for iii in range(0,len(inds)):
    datenow = dates_tpvs[iii]
    datestrnow = str(datenow)
    yyyymmdd_tpvnow = datestrnow[0:8]
    dates_tpvs_nohour.append(yyyymmdd_tpvnow)

dates_tpvs_nohour = np.array(dates_tpvs_nohour).astype('float')


N100_reps = []
N250_reps = []
N500_reps = []
N750_reps = []
N1000_reps = []
N1250_reps = []
N1500_reps = []
N1750_reps = []
N2000_reps = []
N2500_reps = []
N3000_reps = []
N3500_reps = []
N4000_reps = []
N4500_reps = []
N5000_reps = []
N6000_reps = []
for iternow in range(0,numiters):

    rand_year = np.random.randint(random_year_range[0],high=random_year_range[1]+1,size=Nposs)
    rand_month = np.random.choice(month_choices,size=Nposs)
    rand_hour = np.random.choice(hour_choices,size=Nposs)


    distarr = []
    for jjj in range(0,Nposs):

	date_record = dateevents[jjj]
	date_record_x100 = date_record*100
	event_latnow = event_lats[jjj]
	event_lonnow = event_lons[jjj]
	event_magnow = event_mags[jjj]

	date_record_str = str(date_record)	
	yyyy_event = date_record_str[0:4]
	mm_event = date_record_str[4:6]
	dd_event = date_record_str[6:8]
	hh_event = date_record_str[8:10]

	rand_yearnow = rand_year[jjj]
	rand_monthnow = rand_month[jjj]
	if date_format == 'yyyymmddhh':
	    rand_hournow = rand_hour[jjj]
	else:
	    rand_hournow = 00
	if ( (rand_monthnow == 1) or (rand_monthnow == 3) or (rand_monthnow == 5) or (rand_monthnow == 7) or (rand_monthnow == 8) or (rand_monthnow == 10) or (rand_monthnow == 12)):
	    rand_daynow = np.random.randint(0,32)
	elif ( (rand_monthnow == 4) or (rand_monthnow == 6) or (rand_monthnow == 9) or (rand_monthnow == 11) ):
	    rand_daynow = np.random.randint(0,31)
	elif ( rand_monthnow == 2):
	    rand_daynow = np.random.randint(0,29)    

	#print(yyyy_event,mm_event,dd_event,hh_event)
	#print(rand_yearnow,rand_monthnow,rand_daynow,rand_hournow)	
	
	if rand_monthnow < 10:
	    rand_monthnow_str = '0' + str(rand_monthnow)
	else:
	    rand_monthnow_str = str(rand_monthnow)
	if rand_daynow < 10:
	    rand_daynow_str = '0' + str(rand_daynow)
	else:
	    rand_daynow_str = str(rand_daynow)	
	    
	if date_format == 'yyyymmdd':
	    date_record_random = str(rand_yearnow) + rand_monthnow_str + rand_daynow_str + '00'
        else:
	    if rand_hournow < 10:
	        rand_hournow_str = '0' + str(rand_hournow)
	    else:
	        rand_hournow_str = str(rand_hournow)
	    date_record_random = str(rand_yearnow) + rand_monthnow_str + rand_daynow_str + rand_hournow_str
	#print(date_record_random)
        [indnow] = np.where(int(date_record_random) == dates_tpvs)
	#[indnow] = np.where(date_record_x100 == dates_tpvs)
	#print(len(indnow))

	dates_tpvs_now = dates_tpvs[indnow]
	lats_tpvs_now = lats_tpvs[indnow]
	lons_tpvs_now = lons_tpvs[indnow]
	distsave = 9999.     
	for iii in range(0,len(indnow)):
            distance = abs(tpvmod.rEarth*tpvmod.distance_on_unit_sphere(event_latnow,event_lonnow,lats_tpvs_now[iii],lons_tpvs_now[iii]))
            distnow = distance/1000.

	    if iii == 0:
		distsave = distnow
		dist_prev = distnow
	    else:
		if distnow < dist_prev:
	            distsave = distnow	
	distarr.append(distsave)
    #print(distarr)

    distarr = np.array(distarr)
    mstats(distarr)

    del indnow
    [indnow] = np.where(distarr<100.)
    N100 = len(indnow)
    
    del indnow
    [indnow] = np.where(distarr<250.)
    N250 = len(indnow)
    
    del indnow
    [indnow] = np.where(distarr<500.)
    N500 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<750.)
    N750 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<1000.)
    N1000 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<1250.)
    N1250 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<1500.)
    N1500 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<1750.)
    N1750 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<2000.)
    N2000 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<2500.)
    N2500 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<3000.)
    N3000 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<3500.)
    N3500 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<4000.)
    N4000 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<4500.)
    N4500 = len(indnow)

    del indnow
    [indnow] = np.where(distarr<5000.)
    N5000 = len(indnow)
    
    del indnow
    [indnow] = np.where(distarr<6000.)
    N6000 = len(indnow)


    N100_reps.append(N100)
    N250_reps.append(N250)
    N500_reps.append(N500)
    N750_reps.append(N750)
    N1000_reps.append(N1000)
    N1250_reps.append(N1250)
    N1500_reps.append(N1500)
    N1750_reps.append(N1750)
    N2000_reps.append(N2000)
    N2500_reps.append(N2500)
    N3000_reps.append(N3000)
    N3500_reps.append(N3500)
    N4000_reps.append(N4000)
    N4500_reps.append(N4500)
    N5000_reps.append(N5000)
    N6000_reps.append(N6000)


N100_avg = np.nanmean(N100_reps)
N100_5 = np.percentile(N100_reps, percentiles[0])
N100_95 = np.percentile(N100_reps, percentiles[1])

N250_avg = np.nanmean(N250_reps)
N250_5 = np.percentile(N250_reps, percentiles[0])
N250_95 = np.percentile(N250_reps, percentiles[1])

N500_avg = np.nanmean(N500_reps)
N500_5 = np.percentile(N500_reps, percentiles[0])
N500_95 = np.percentile(N500_reps, percentiles[1])

N750_avg = np.nanmean(N750_reps)
N750_5 = np.percentile(N750_reps, percentiles[0])
N750_95 = np.percentile(N750_reps, percentiles[1])

N1000_avg = np.nanmean(N1000_reps)
N1000_5 = np.percentile(N1000_reps, percentiles[0])
N1000_95 = np.percentile(N1000_reps, percentiles[1])

N1250_avg = np.nanmean(N1250_reps)
N1250_5 = np.percentile(N1250_reps, percentiles[0])
N1250_95 = np.percentile(N1250_reps, percentiles[1])

N1500_avg = np.nanmean(N1500_reps)
N1500_5 = np.percentile(N1500_reps, percentiles[0])
N1500_95 = np.percentile(N1500_reps, percentiles[1])

N1750_avg = np.nanmean(N1750_reps)
N1750_5 = np.percentile(N1750_reps, percentiles[0])
N1750_95 = np.percentile(N1750_reps, percentiles[1])

N2000_avg = np.nanmean(N2000_reps)
N2000_5 = np.percentile(N2000_reps, percentiles[0])
N2000_95 = np.percentile(N2000_reps, percentiles[1])

N2500_avg = np.nanmean(N2500_reps)
N2500_5 = np.percentile(N2500_reps, percentiles[0])
N2500_95 = np.percentile(N2500_reps, percentiles[1])

N3000_avg = np.nanmean(N3000_reps)
N3000_5 = np.percentile(N3000_reps, percentiles[0])
N3000_95 = np.percentile(N3000_reps, percentiles[1])

N3500_avg = np.nanmean(N3500_reps)
N3500_5 = np.percentile(N3500_reps, percentiles[0])
N3500_95 = np.percentile(N3500_reps, percentiles[1])

N4000_avg = np.nanmean(N4000_reps)
N4000_5 = np.percentile(N4000_reps, percentiles[0])
N4000_95 = np.percentile(N4000_reps, percentiles[1])

N4500_avg = np.nanmean(N4500_reps)
N4500_5 = np.percentile(N4500_reps, percentiles[0])
N4500_95 = np.percentile(N4500_reps, percentiles[1])

N5000_avg = np.nanmean(N5000_reps)
N5000_5 = np.percentile(N5000_reps, percentiles[0])
N5000_95 = np.percentile(N5000_reps, percentiles[1])

N6000_avg = np.nanmean(N6000_reps)
N6000_5 = np.percentile(N6000_reps, percentiles[0])
N6000_95 = np.percentile(N6000_reps, percentiles[1])

#alpha = (confidence_interval_limits[0]*2.0)/100.
#boot_means,ci_limits_out = wm.bootstrap(N1000_reps,10000,alpha)
#print(ci_limits_out)
#exit()

event_outfile.write('\n')
print(Nposs, N100_avg, N250_avg, N500_avg, N750_avg, N1000_avg, N1250_avg, N1500_avg, N1750_avg, N2000_avg, N2500_avg, N3000_avg, N3500_avg, N4000_avg, N4500_avg, N5000_avg, N6000_avg)
print(Nposs, N100_5, N250_5, N500_5, N750_5, N1000_5, N1250_5, N1500_5, N1750_5, N2000_5, N2500_5, N3000_5, N3500_5, N4000_5, N4500_5, N5000_5, N6000_5)

#for jjj in range(0,Nposs):
event_outfile.write('%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' % (Nposs, N100_avg, N250_avg, N500_avg, N750_avg, N1000_avg, N1250_avg, N1500_avg, N1750_avg, N2000_avg, N2500_avg, N3000_avg, N3500_avg, N4000_avg, N4500_avg, N5000_avg, N6000_avg))  
event_outfile.write('%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' % (Nposs, N100_5, N250_5, N500_5, N750_5, N1000_5, N1250_5, N1500_5, N1750_5, N2000_5, N2500_5, N3000_5, N3500_5, N4000_5, N4500_5, N5000_5, N6000_5)) 
event_outfile.write('%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' % (Nposs, N100_95, N250_95, N500_95, N750_95, N1000_95, N1250_95, N1500_95, N1750_95, N2000_95, N2500_95, N3000_95, N3500_95, N4000_95, N4500_95, N5000_95, N6000_95))  
