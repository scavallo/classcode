#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from scipy import ndimage
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")



###################################
# Set user options
###################################


figname_description = 'test'
#figname_description = 'allarctic_rankedevents_19792016_08'

hinc = 6.
min_lifetime_days = 2
minimum_latitude = 30.
maximum_latitude = 60.
label_fontsize = 18
month_filt = [6,8]
years_filt = [2007,2016]
num_ranks = 15
pres_ref = 980.
dumval = -999.00

#fdir1 = '/data2/scavallo/nnrp/cyclones_nsidc/'
#fdir2 = '/data2/scavallo/nnrp/cyclones_nsidc/'
fdir1 = '/Users/scavallo/Documents/data/nnrp/cyclones_nsidc/'
fdir2 = '/Users/scavallo/Documents/data/nnrp/cyclones_nsidc/'


textfile1 = 'ncepstorms_allarctic_1958_2016_jja.txt'
#textfile2 = 'ncepstorms_allarctic_1958_2016_jja.txt'
#textfile1 = 'ncepstorms_allarctic_highseaiceIQR_1958_2016_jja.txt'
#textfile2 = 'ncepstorms_allarctic_lowseaiceIQR_1958_2016_jja.txt'
#textfile1 = 'ncepstorms_allarctic_1981_2010_son.txt'
#textfile1 = 'ncepstorms_allarctic_1958_2016.txt'
#textfile1 = 'ncepstorms_1958_2016.txt.tracksFmt.txt'
#textfile1 = 'ncepstorms_60N65_1958_2016.txt'
#textfile2 = 'ncepstorms_allarctic_2007_2016.txt'
#textfile1 = 'sfccyclone_tracks_seaicelossevents_5percentile_jja_1979010100_2014123118.dat'
#textfile2 = 'ncepstorms_allarctic_1958_2016_jja.txt'

write_outfile = 'True'
textfile2 = 'allmidlatcyclones_minS60_genS60_jja_2007_2016_top15.txt'
#textfile2 = 'allcyclones_jja_2007_2016_top15.txt'

legend2text = 'Top Arctic cyclones'
legend1text = 'Arctic cyclones, 2007-2016'
#legend1text = 'Arctic cyclones, DJF'
#legend2text = 'Arctic cyclones, JJA'
#legend1text = 'All Arctic cyclones, JJA'
#legend2text = '60N65 Arctic cyclones, JJA'
#legend1text = '1981-2010'
#legend2text = '2007-2016'
#legend1text = 'High sea ice IQR, JJA'
#legend2text = 'Low sea ice IQR, JJA'
#legend1text = 'Event JJA'
#legend2text = 'Climate JJA'


imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
###################################
# END user options
###################################

fpath = fdir1 + textfile1
aa = np.loadtxt(fpath, skiprows=1)
datelist = aa[:,0]
lat = aa[:,1]	
lon = aa[:,2]
slpmin = aa[:,3]	    
slpamp = aa[:,4]
vort_radius = aa[:,5]	    
vort_ptend = aa[:,6]

nrows, ncols = np.shape(aa)
print(nrows, ncols)

tinds = np.where(datelist>10000)
lats_all = lat[tinds]
lons_all = lon[tinds]
Y = lats_all[0:-1]
X = lons_all[0:-1]
lat2 = lats_all[1:]
lon2 = lons_all[1:]

slp_min = []
year_min = []
month_min = []
day_min = []
yyyymmdd_min = []
lat_subset = []
lon_subset = []
radius_subset = []
ptend_subset = []
amps_subset = []
startind_subset = []
endind_subset = []

ntracks_min = (24./hinc)*min_lifetime_days


tt = 0
tracknum = -1
while tt < nrows:
	if datelist[tt] < 10000:
		ntracks = datelist[tt]
		nhours = ntracks*hinc
		sind = tt+1
		eind = np.int(sind+ntracks)

		datestart = str(datelist[tt+1])
		slpnow = slpmin[sind:eind]
		[indnow] = np.where(slpnow == np.nanmin(slpnow))
		latnow = lat[sind:eind]
		#maxlatnow = np.max(latnow)	
		maxlatnow = np.max(latnow[indnow][0])
		genlatnow = lat[sind]
		yearstart = datestart[0:4] 
	
	
		#if ( (lat[sind] >= lat_origin_range[0]) and (lat[sind] <= lat_origin_range[1]) and (ntracks >= ntracks_min) ):
		if ( (maxlatnow >= minimum_latitude) and (maxlatnow <= maximum_latitude) and (genlatnow < maximum_latitude) and (genlatnow >= minimum_latitude) and (ntracks >= ntracks_min) and (int(yearstart)>=years_filt[0]) and (int(yearstart)<=years_filt[1]) ):
			tracknum += 1
	    	#print sind,eind,tracknum
	    	#print lat[sind],lat[eind]
			
			lonnow = lon[sind:eind]	    
			ampsnow = slpamp[sind:eind]	    
			radsnow = vort_radius[sind:eind]	    
			datesnow = datelist[sind:eind]
			ptendsnow = vort_ptend[sind:eind]	
	    
	    	    
	    
			[minind] = np.where(slpnow==np.nanmin(slpnow))		
			datemin = str(datesnow[minind[0]])		
			yearminnow = datemin[0:4]
			monthminnow = datemin[4:6]		
			dayminnow = datemin[6:8]	
			yyyymmddnow = datemin[0:8]		    
	    
	    
			if ( (int(monthminnow) >= month_filt[0]) and (int(monthminnow) <= month_filt[1]) ): 
			#if ( (int(monthminnow) >= month_filt[0]) and (int(monthminnow) <= month_filt[1]) and (int(latnow[minind[0]]) >= 60) ): 	    
			#if ( (int(monthminnow) >= month_filt[0]) and (int(monthminnow) <= month_filt[1]) and (int(latnow[minind[0]]) < 60) ): 	 
			
				slp_min.append(slpnow[minind[0]])
				year_min.append(yearminnow)
				month_min.append(monthminnow)	
				day_min.append(dayminnow)
				yyyymmdd_min.append(yyyymmddnow)   
				lat_subset.append(latnow[minind[0]])
				lon_subset.append(lonnow[minind[0]])
				radius_subset.append(radsnow[minind[0]])
				ptend_subset.append(ptendsnow[minind[0]])
				amps_subset.append(ampsnow[minind[0]])
				startind_subset.append(sind)
				endind_subset.append(eind)
	     
	    
	tt += 1

	
slp_arr = np.array(slp_min).astype('f')	
years_arr = np.array(year_min).astype(int)
month_arr = np.array(month_min).astype(int)
day_arr = np.array(day_min).astype(int)
yyyymmdd_arr = np.array(yyyymmdd_min)



ranked_slp = np.sort(slp_arr)
ranked_slp_indices = np.argsort(slp_arr)
yyyymmdd_ranked = yyyymmdd_arr[ranked_slp_indices]
ranked_slp_values = slp_arr[ranked_slp_indices]

amps_subset_arr = np.array(amps_subset)
radius_subset_arr = np.array(radius_subset)
ptend_subset_arr = np.array(ptend_subset)
lat_subset_arr = np.array(lat_subset)
lon_subset_arr = np.array(lon_subset)
startind_subset_arr = np.array(startind_subset)
endind_subset_arr = np.array(endind_subset)
ranked_amps = amps_subset_arr[ranked_slp_indices]
ranked_radii = radius_subset_arr[ranked_slp_indices]
ranked_ptends = ptend_subset_arr[ranked_slp_indices]
ranked_lats = lat_subset_arr[ranked_slp_indices]
ranked_lons = lon_subset_arr[ranked_slp_indices]
ranked_startinds = startind_subset_arr[ranked_slp_indices]
ranked_endinds = endind_subset_arr[ranked_slp_indices]

#print(ranked_startinds)
#print(ranked_endinds)

if write_outfile == 'True':
    fpath_out = fdir2 + textfile2
    outfile = open(fpath_out,'w')
    outfile.write('date lat lon slpmin amplitude radius ptendency')
    outfile.write('\n')

    #for ii in range(0,len(ranked_startinds)):
    for ii in range(0,num_ranks):
		
        sindnow = ranked_startinds[ii]
        eindnow = ranked_endinds[ii]
        ntracks = eindnow - sindnow 
        wcount = sindnow
        outfile.write('%-10d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n' % (ntracks, dumval, dumval, dumval, dumval, dumval, dumval))
        while wcount < eindnow:
            outfile.write('%-10s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n' % (str(int(datelist[wcount])), lat[wcount], lon[wcount], slpmin[wcount], slpamp[wcount], vort_radius[wcount], vort_ptend[wcount]))
			
            wcount += 1
############################################
# Figure 1
############################################
numinds = num_ranks
values_plot = ranked_slp_values[0:numinds]
times_plot = yyyymmdd_ranked[0:numinds]

n_groups = np.size(values_plot)
index = np.arange(n_groups)
bar_width = 0.5
opacity = 0.8
y_pos = np.arange(len(values_plot))

ylims = [np.min(values_plot)-1,np.max(values_plot)+1]

x0 = np.zeros(len(values_plot))
fig, ax = plt.subplots()

rects1 = plt.bar(index,values_plot,bar_width,alpha=opacity,align='center',edgecolor='k',linewidth=1.0,color='b',label=legend1text)
#plt.axvline(x=0.22058956)
ax.grid(True,linestyle='-')
plt.ylim([ylims[0],ylims[1]])
plt.xlim([-0.5,n_groups-0.5])
plt.xticks(y_pos,times_plot,rotation=90)
plt.ylabel('Minimum SLP (hPa)',fontsize=label_fontsize)
save_name = figname_description + '.png'
plt.savefig(imagedir + save_name, bbox_inches='tight')


############################################
# Figure 2
############################################
[indsbelow] = np.where(ranked_slp_values<=pres_ref)
numinds = indsbelow[-1]
values_plot2 = ranked_slp_values[0:numinds+1]
times_plot2 = yyyymmdd_ranked[0:numinds+1]

ref_percent = ((numinds+2.)/len(ranked_slp_values))*100.

n_groups2 = np.size(values_plot2)
index2 = np.arange(n_groups2)
bar_width = 0.5
opacity = 0.8
y_pos2 = np.arange(len(values_plot2))

ylims = [np.min(values_plot2)-1,np.max(values_plot2)+1]
#legend2text = 'Rank = ' + str(numinds+1) + ' (' + str(ref_percent) + '%)' 
legend2text = "Rank = %d (%5.2f %%)" % (numinds+2, ref_percent)

x0 = np.zeros(len(values_plot2))
fig, ax = plt.subplots()

rects1 = plt.bar(index2,values_plot2,bar_width,alpha=opacity,align='center',edgecolor='k',linewidth=1.0,color='b',label=legend2text)
#plt.axvline(x=0.22058956)
ax.grid(True,linestyle='-')
plt.ylim([ylims[0],ylims[1]])
plt.xlim([-0.5,n_groups-0.5])
plt.xticks(y_pos2,times_plot2,rotation=90)
legend = ax.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
plt.ylabel('Minimum SLP (hPa)',fontsize=label_fontsize)
save_name = figname_description + '_presref_' + str(int(pres_ref)) + 'hPa.png'
plt.savefig(imagedir + save_name, bbox_inches='tight')
#plt.show()

############################################
# Figure 3
############################################
numinds = num_ranks
values_plot = ranked_lats[0:numinds]
times_plot = yyyymmdd_ranked[0:numinds]

n_groups = np.size(values_plot)
index = np.arange(n_groups)
bar_width = 0.5
opacity = 0.8
y_pos = np.arange(len(values_plot))

ylims = [np.min(values_plot)-1,np.max(values_plot)+1]

x0 = np.zeros(len(values_plot))
fig, ax = plt.subplots()

rects1 = plt.bar(index,values_plot,bar_width,alpha=opacity,align='center',edgecolor='k',linewidth=1.0,color='b',label=legend1text)
#plt.axvline(x=0.22058956)
ax.grid(True,linestyle='-')
plt.ylim([ylims[0],ylims[1]])
plt.xlim([-0.5,n_groups-0.5])
plt.xticks(y_pos,times_plot,rotation=90)
plt.ylabel(r'Latitude ($^{\circ}$N)',fontsize=label_fontsize)
#save_name = figname_description + '.png'
#plt.savefig(imagedir + save_name, bbox_inches='tight')

############################################
# Figure 4
############################################
numinds = num_ranks
values_plot = ranked_amps[0:numinds]
times_plot = yyyymmdd_ranked[0:numinds]

n_groups = np.size(values_plot)
index = np.arange(n_groups)
bar_width = 0.5
opacity = 0.8
y_pos = np.arange(len(values_plot))

ylims = [np.min(values_plot)-1,np.max(values_plot)+1]

x0 = np.zeros(len(values_plot))
fig, ax = plt.subplots()

rects1 = plt.bar(index,values_plot,bar_width,alpha=opacity,align='center',edgecolor='k',linewidth=1.0,color='b',label=legend1text)
#plt.axvline(x=0.22058956)
ax.grid(True,linestyle='-')
plt.ylim([ylims[0],ylims[1]])
plt.xlim([-0.5,n_groups-0.5])
plt.xticks(y_pos,times_plot,rotation=90)
plt.ylabel('Amplitude (K)',fontsize=label_fontsize)
#save_name = figname_description + '.png'
#plt.savefig(imagedir + save_name, bbox_inches='tight')
plt.show()
