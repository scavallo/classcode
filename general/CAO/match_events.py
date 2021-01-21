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
#eventpath = '/data2/scavallo/era_interim/track_files/okc_snowdates_50percentile.dat'
#eventpath = '/data1/scavallo/data/seaice/seaice_loss_events/paper_data/rapid_seaice_loss_events_int_withboth_annual_bwfilter_3d_5percentile.dat'
#eventpath = '/data1/scavallo/data/seaice/seaice_loss_events/5percentile_2017/rapid_seaice_loss_events_int_withboth_jja_3d.dat'
#trackpath_in = '/data2/scavallo/era_interim/track_files/trth_tracks_jja_1979010100_2014123118.dat'
#trackpath_in = '/data2/scavallo/era_interim/track_files/trth_tracks_all_1979010100_2016123118.dat'
#trackpath_out = '/data1/scavallo/data/seaice/seaice_loss_events/5percentile_2017/trth_tracks_seaicelossevents_5percentile_jja_1979010100_2016123118.dat'
#trackpath_in = '/data2/scavallo/nnrp/cyclones_nsidc/ncepstorms_allarctic_1958_2016.txt'
#trackpath_out = '/data1/scavallo/data/seaice/seaice_loss_events/5percentile_2017/sfccyclone_tracks_seaicelossevents_5percentile_jja_1979010100_2016123118.dat'
#trackpath_out = '/data1/scavallo/data/seaice/seaice_loss_events/paper_data/sfccyclone_tracks_bwfilter_5percentile_son_1979010100_2016123118.dat'
#trackpath_out = '/data2/scavallo/era_interim/track_files/trth_tracks_all_okc_snowdates_50percentile.dat'

#trackpath_in = '/data2/scavallo/era_interim/track_files/trth_tracks_60N65_1979010100_2016123118.dat'
#trackpath_out = '/data1/scavallo/data/seaice/seaice_loss_events/paper_data/trth_tracks_60N65_bwfilter_5percentile_son_1979010100_2016123118.dat'
#trackpath_out = '/data1/scavallo/data/seaice/seaice_loss_events/paper_data/trth_tracks_60N65_bwfilter_dtUnique01_5percentile_jja_1979010100_2016123118.dat'


#trackpath_in = '/data2/scavallo/era_interim/track_files/tpvs_60N65_djf_CAOrange_500km_1979010100_2018123118.dat'
#eventpath_in = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/CAOinfo.txt'
#eventpath_out1 = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/CAOinfo_hits_60N65_500km.txt'
#eventpath_out2 = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/CAOinfo_nulls_60N65_500km.txt'
#eventpath_out3 = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/CAOinfo_all_60N65_500km.txt'

#trackpath_in = '/data1/scavallo/data/seaice/VRILEs/1979_2019/tpvs_60N65_jja_VRILErange_500km_1979010100_2018123118.dat'
#eventpath_in = '/data1/scavallo/data/seaice/VRILEs/1979_2019/VRILE_masks_jja_bwfilter_3d_5percentile_dtUnique01_1979_2019.dat'
#eventpath_out1 = '/data1/scavallo/data/seaice/VRILEs/1979_2019/TPVinfo_hits_60N65_500km.txt'
#eventpath_out2 = '/data1/scavallo/data/seaice/VRILEs/1979_2019/TPVinfo_nulls_60N65n_500km.txt'
#eventpath_out3 = '/data1/scavallo/data/seaice/VRILEs/1979_2019/TPVinfo_all_60N65_500km.txt'


if 1 == 0:
    date_format = 'yyyymmddhh'
    #trackpath_in = '/data1/scavallo/data/seaice/VRILEs/1979_2019/tpvs_60N65_jja_VRILErange_500km_1979010100_2018123118.dat'
    eventpath_in = '/data1/scavallo/data/seaice/VRILEs/1979_2019/VRILE_masks_jja_bwfilter_3d_5percentile_dtUnique01_1979_2019.dat'
    trackdir = '/data1/scavallo/data/seaice/VRILEs/1979_2019/'
    eventdir = '/data1/scavallo/data/seaice/VRILEs/1979_2019/'
    #trackfile_prefix = 'tpvs_60Ngen_jja_VRILErange_'
    trackfile_desc = '60Ngen_jja_VRILErange_'
    eventfile_prefix = 'TPVinfo_'
    eventfile_desc = '60Ngen_'
else:
    date_format = 'yyyymmdd'
    #trackpath_in = '/data2/scavallo/era_interim/track_files/tpvs_60N65_djf_CAOrange_500km_1979010100_2018123118.dat'
    eventpath_in = '/data1/scavallo/data/cold_air_outbreaks/CAOinfo.txt'
    trackdir = '/data1/scavallo/data/cold_air_outbreaks/'
    eventdir = '/data1/scavallo/data/cold_air_outbreaks/test/'    
    trackfile_desc = '60N65_djf_CAOrange_'
    eventfile_prefix = 'CAOinfo_'
    eventfile_desc = '60N65_'


suffs = [500,1000,1500,1750,2000,2500,3000,3500,4000,5000]


# Read in event dates
bb = np.loadtxt(eventpath_in, skiprows=1)       
#dateevents = bb[:,0]
dateevents = bb[:,0]
event_lats = bb[:,1]
event_lons = bb[:,2]
event_mags = bb[:,3]
Ntimeranges = len(dateevents)    

print Ntimeranges

dateevents = np.array(dateevents)
event_lats = np.array(event_lats)
event_lons = np.array(event_lons,dtype=float)
event_mags = np.array(event_mags)

[inds] = np.where(event_lons<0)
event_lons[inds] = event_lons[inds] + 360.

eventpath_out_stats = eventdir + eventfile_prefix + 'stats_' + eventfile_desc[:-1] + '.txt'
event_outfile_stats = open(eventpath_out_stats,'w')
event_outfile_stats.write('Threshold (km) hits nulls total')
event_outfile_stats.write('\n')

for tt in range(0,len(suffs)):

    trackpath_in_now = trackdir + 'tpvs_' + trackfile_desc + str(suffs[tt]) + 'km_1979010100_2018123118.dat'
    print(trackpath_in_now)
    
    eventpath_out1_now = eventdir + eventfile_prefix + 'hits_' + eventfile_desc + str(suffs[tt]) + 'km.txt'
    eventpath_out2_now = eventdir + eventfile_prefix + 'nulls_' + eventfile_desc + str(suffs[tt]) + 'km.txt'
    eventpath_out3_now = eventdir + eventfile_prefix + 'all_' + eventfile_desc + str(suffs[tt]) + 'km.txt'
        
    
    # Read in full track file
    aa = np.loadtxt(trackpath_in_now, skiprows=1)       
    datelist = aa[:,0]
    lat = aa[:,1]	
    lon = aa[:,2]
    thetamin = aa[:,3]	    
    thetaamp = aa[:,4]
    vort_circ = aa[:,5]	    
    vort_radius = aa[:,6]

    nrows, ncols = np.shape(aa)
    print nrows, ncols

    if date_format == 'yyyymmddhh':
        [inds] = np.where(datelist > 1900000000)
    else:
        [inds] = np.where(datelist > 19000000)

    dates_tpvs = datelist[inds]
    lats_tpvs = lat[inds]
    lons_tpvs = lon[inds]
    dates_tpvs_nohour = []
    for iii in range(0,len(inds)):
	datenow = dates_tpvs[iii]
	datestrnow = str(datenow)
	
	yyyymmddhh_tpvnow = datestrnow[0:10]
	yyyymmdd_tpvnow = datestrnow[0:8]
	#dates_tpvs.append(yyyymmddhh_tpvnow)
	dates_tpvs_nohour.append(yyyymmdd_tpvnow)
	
	

    dates_tpvs = np.array(dates_tpvs).astype('float')
    dates_tpvs_nohour = np.array(dates_tpvs_nohour).astype('float')
    lats_tpvs = np.array(lats_tpvs)
    lons_tpvs = np.array(lons_tpvs)

    event_outfile = open(eventpath_out1_now,'w')
    event_outfile.write('Date Lat Lon Mag')
    event_outfile.write('\n')

    event_outfile_null = open(eventpath_out2_now,'w')
    event_outfile_null.write('Date Lat Lon Mag')
    event_outfile_null.write('\n')

    event_outfile_all = open(eventpath_out3_now,'w')
    event_outfile_all.write('Date Lat Lon Mag')
    event_outfile_all.write('\n')

    nullcount = 0
    hitcount = 0
    for jjj in range(0,Ntimeranges):

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

	#[indnow] = np.where(date_record_x100 == dates_tpvs)
	if date_format == 'yyyymmdd':
	    dates_tpvs_now = dates_tpvs_nohour
	    [indnow] = np.where(date_record == dates_tpvs_now)
	else:
	    dates_tpvs_now = dates_tpvs
	    [indnow] = np.where(date_record == dates_tpvs_now)
	
	#print(jjj,len(indnow))
	if ( len(indnow) > 0 ):
	
	    for iii in range(0,len(indnow)):
	        print(date_record,dates_tpvs_now[indnow[iii]])	    
		print(lats_tpvs[indnow[iii]],lons_tpvs[indnow[iii]])
	        track_latnow = lats_tpvs[indnow[iii]]
		track_lonnow = lons_tpvs[indnow[iii]]
		#[track_ind] = np.where(dates_tpvs_now == date_record)
	        #print(indnow,track_ind)
	        #track_latnow = lats_tpvs[track_ind[0]]
	        #track_lonnow = lons_tpvs[track_ind[0]]
	
                #event_outfile.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, event_latnow, event_lonnow, event_magnow))  
	        event_outfile.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, track_latnow, track_lonnow, event_magnow))
		event_outfile_all.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, track_latnow, track_lonnow, event_magnow))   
	    hitcount += 1
	else:
            event_outfile_null.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, event_latnow, event_lonnow, event_magnow))
	    event_outfile_all.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, event_latnow, event_lonnow, event_magnow))  
	    #event_outfile.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, track_latnow, track_lonnow, event_magnow)) 
	    nullcount += 1

	#event_outfile_all.write('%-8s %7.2f %7.2f %7.2f\n' % (date_record, event_latnow, event_lonnow, event_magnow))  
	

    print(hitcount, nullcount)
    event_outfile_stats.write('%5d %4d %4d %4d\n' % (suffs[tt],hitcount,nullcount,hitcount+nullcount))
    
    event_outfile.close()
    event_outfile_null.close()
    event_outfile_all.close()
event_outfile_stats.close()
