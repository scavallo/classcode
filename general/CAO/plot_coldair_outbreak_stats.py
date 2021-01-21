# imports
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
import math
import matplotlib.colors as col
from matplotlib.patches import Polygon

import warnings
warnings.filterwarnings("ignore")

###################################
# BEGIN user options
###################################
#fdir_in = '/data2/scavallo/era_interim/track_files/'
fdir_in = '/data1/scavallo/data/cold_air_outbreaks/'
#fdir_in = '/data1/scavallo/data/seaice/VRILEs/1979_2019/'
#fdir_in_random = '/data1/scavallo/data/seaice/VRILEs/1979_2019/'
fdir_in_random = '/data1/scavallo/data/cold_air_outbreaks/'
imagedir = '/home/scavallo/scripts/python_scripts/images/'
#imagedir = '/home/scavallo/chestnut_files/papers/TPV_jan2019_lillo/images/'

plot_pref1 = True
plot_pref2 = True
prefix_main = 'TPVinfo'
pref1 = '60Ngen'
pref2 = '60N65'
labeltext_pref1 = pref1
labeltext_pref2 = pref2

eventfile_in1 = 'CAOinfo_hits_60N65_500km.txt'
eventfile_in2 = 'CAOinfo_hits_60N65_1000km.txt'
eventfile_in3 = 'CAOinfo_hits_60N65_1500km.txt'
eventfile_in4 = 'CAOinfo_hits_60N65_1750km.txt'
eventfile_in5 = 'CAOinfo_hits_60N65_2000km.txt'
eventfile_in6 = 'CAOinfo_hits_60N65_2500km.txt'
eventfile_in7 = 'CAOinfo_hits_60N65_3000km.txt'
eventfile_in_full = 'CAOinfo.txt'
#eventfile_in_full = 'VRILE_masks_jja_bwfilter_3d_5percentile_dtUnique01_1979_2019.dat'

random_60Ngen = 'test_CAO_60Ngen_randomTPV_with95confidence.txt'
random_60N65 = 'test_CAO_60N65_randomTPV_with95confidence.txt'

#random_60Ngen = 'VRILE_60Ngen_randomTPV_with95confidence.txt'
#random_60N65 = 'VRILE_60N65_randomTPV_with95confidence.txt'

hinc = 6
min_lifetime_days = 2

highlight_state = 'False'
proj_latlon = [30. , 270.]

label_fontsize = 18
track_thickness = 4
num_passes_smooth = 2
legendtext = ['Test']
###################################
# END user options
###################################
bb = np.loadtxt(fdir_in_random + random_60Ngen, skiprows=1)       
Nposs = bb[0,0]
N500 = bb[0,3]
N1000 = bb[0,5]
N1500 = bb[0,7]
N1750 = bb[0,8]
N2000 = bb[0,9]
N2500 = bb[0,10]
N3000 = bb[0,11]
N3500 = bb[0,12]
N4000 = bb[0,13]
N5000 = bb[0,15]
percent_climo_60Ngen = [(N500/Nposs)*100., (N1000/Nposs)*100.,(N1500/Nposs)*100.,(N1750/Nposs)*100.,(N2000/Nposs)*100.,(N2500/Nposs)*100.,(N3000/Nposs)*100.,(N3500/Nposs)*100.,(N4000/Nposs)*100.,(N5000/Nposs)*100.]

N500a = bb[1,3]
N1000a = bb[1,5]
N1500a = bb[1,7]
N1750a = bb[1,8]
N2000a = bb[1,9]
N2500a = bb[1,10]
N3000a = bb[1,11]
N3500a = bb[1,12]
N4000a = bb[1,13]
N5000a = bb[1,15]
percent_climo_60Ngen_5per = [(N500a/Nposs)*100., (N1000a/Nposs)*100.,(N1500a/Nposs)*100.,(N1750a/Nposs)*100.,(N2000a/Nposs)*100.,(N2500a/Nposs)*100.,(N3000a/Nposs)*100.,(N3500a/Nposs)*100.,(N4000a/Nposs)*100.,(N5000a/Nposs)*100.]

N500b = bb[2,3]
N1000b = bb[2,5]
N1500b = bb[2,7]
N1750b = bb[2,8]
N2000b = bb[2,9]
N2500b = bb[2,10]
N3000b = bb[2,11]
N3500b = bb[2,12]
N4000b = bb[2,13]
N5000b = bb[2,15]
percent_climo_60Ngen_95per = [(N500b/Nposs)*100., (N1000b/Nposs)*100.,(N1500b/Nposs)*100.,(N1750b/Nposs)*100.,(N2000b/Nposs)*100.,(N2500b/Nposs)*100.,(N3000b/Nposs)*100.,(N3500b/Nposs)*100.,(N4000b/Nposs)*100.,(N5000b/Nposs)*100.]

bb = np.loadtxt(fdir_in_random + random_60N65, skiprows=1)       
Nposs = bb[0,0]
N500 = bb[0,3]
N1000 = bb[0,5]
N1500 = bb[0,7]
N1750 = bb[0,8]
N2000 = bb[0,9]
N2500 = bb[0,10]
N3000 = bb[0,11]
N3500 = bb[0,12]
N4000 = bb[0,13]
N5000 = bb[0,15]
percent_climo_60N65 = [(N500/Nposs)*100., (N1000/Nposs)*100.,(N1500/Nposs)*100.,(N1750/Nposs)*100.,(N2000/Nposs)*100.,(N2500/Nposs)*100.,(N3000/Nposs)*100.,(N3500/Nposs)*100.,(N4000/Nposs)*100.,(N5000/Nposs)*100.]

N500a = bb[1,3]
N1000a = bb[1,5]
N1500a = bb[1,7]
N1750a = bb[1,8]
N2000a = bb[1,9]
N2500a = bb[1,10]
N3000a = bb[1,11]
N3500a = bb[1,12]
N4000a = bb[1,13]
N5000a = bb[1,15]
percent_climo_60N65_5per = [(N500a/Nposs)*100., (N1000a/Nposs)*100.,(N1500a/Nposs)*100.,(N1750a/Nposs)*100.,(N2000a/Nposs)*100.,(N2500a/Nposs)*100.,(N3000a/Nposs)*100.,(N3500a/Nposs)*100.,(N4000a/Nposs)*100.,(N5000a/Nposs)*100.]

N500b = bb[2,3]
N1000b = bb[2,5]
N1500b = bb[2,7]
N1750b = bb[2,8]
N2000b = bb[2,9]
N2500b = bb[2,10]
N3000b = bb[2,11]
N3500b = bb[2,12]
N4000b = bb[2,13]
N5000b = bb[2,15]
percent_climo_60N65_95per = [(N500b/Nposs)*100., (N1000b/Nposs)*100.,(N1500b/Nposs)*100.,(N1750b/Nposs)*100.,(N2000b/Nposs)*100.,(N2500b/Nposs)*100.,(N3000b/Nposs)*100.,(N3500b/Nposs)*100.,(N4000b/Nposs)*100.,(N5000b/Nposs)*100.]

if ( (plot_pref1 == True) and (plot_pref2 == True) ):
    numiter = 2
else:
    numiter = 1

# Read in event dates
bb = np.loadtxt(fdir_in + eventfile_in_full, skiprows=1)       
dateevents = bb[:,0]
event_lats = bb[:,1]
event_lons = bb[:,2]
event_mags = bb[:,3]
Nposs = np.float(len(dateevents)-2) # Subtract 2 because of 2019      

for ii in range(1,numiter+1):    
    print(ii)
    if ii == 1:
        prefnow = pref1
    else:
        prefnow = pref2
    
    eventfile_in_now = prefix_main + '_hits_' + prefnow + '_500km.txt'    
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N1 = np.float(len(dateevents))

    eventfile_in_now = prefix_main + '_hits_' + prefnow + '_1000km.txt'  
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N2 = np.float(len(dateevents))

    eventfile_in_now = prefix_main + '_hits_' + prefnow + '_1500km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N3 = np.float(len(dateevents))

    eventfile_in_now = prefix_main + '_hits_' + prefnow + '_1750km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N4 = np.float(len(dateevents))

    eventfile_in_now = prefix_main + '_hits_' + prefnow + '_2000km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N5 = np.float(len(dateevents))

    eventfile_in_now =  prefix_main + '_hits_' + prefnow + '_2500km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N6 = np.float(len(dateevents))

    eventfile_in_now =  prefix_main + '_hits_' + prefnow + '_3000km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N7 = np.float(len(dateevents))

    eventfile_in_now =  prefix_main + '_hits_' + prefnow + '_3500km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N8 = np.float(len(dateevents))
    
    eventfile_in_now =  prefix_main + '_hits_' + prefnow + '_4000km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N9 = np.float(len(dateevents))

    eventfile_in_now =  prefix_main + '_hits_' + prefnow + '_5000km.txt' 
    del dateevents, event_lats, event_lons, event_mags
    bb = np.loadtxt(fdir_in + eventfile_in_now, skiprows=1)       
    dateevents = bb[:,0]
    event_lats = bb[:,1]
    event_lons = bb[:,2]
    event_mags = bb[:,3]
    N10 = np.float(len(dateevents))

    
    print(Nposs)
    print(N1,N2,N3,N4,N5,N6,N7,N8,N9,N10)
    
    if ii == 1:
	percent_500km_pref1  = (N1/Nposs)*100.
	percent_1000km_pref1  = (N2/Nposs)*100.
	percent_1500km_pref1  = (N3/Nposs)*100.
	percent_1750km_pref1  = (N4/Nposs)*100.
	percent_2000km_pref1  = (N5/Nposs)*100.
	percent_2500km_pref1  = (N6/Nposs)*100.
	percent_3000km_pref1  = (N7/Nposs)*100.
	percent_3500km_pref1  = (N8/Nposs)*100.
	percent_4000km_pref1  = (N9/Nposs)*100.
	percent_5000km_pref1  = (N10/Nposs)*100.
    elif ii == 2:
	percent_500km_pref2  = (N1/Nposs)*100.
	percent_1000km_pref2  = (N2/Nposs)*100.
	percent_1500km_pref2  = (N3/Nposs)*100.
	percent_1750km_pref2  = (N4/Nposs)*100.
	percent_2000km_pref2  = (N5/Nposs)*100.
	percent_2500km_pref2  = (N6/Nposs)*100.
	percent_3000km_pref2  = (N7/Nposs)*100.        
        percent_3500km_pref2  = (N8/Nposs)*100. 
	percent_4000km_pref2  = (N9/Nposs)*100. 
	percent_5000km_pref2  = (N10/Nposs)*100. 


print(percent_1000km_pref1,percent_1500km_pref1,percent_1750km_pref1,percent_2000km_pref1,percent_2500km_pref1,percent_3000km_pref1,percent_3500km_pref1,percent_4000km_pref1,percent_5000km_pref1)
if numiter == 2: 
    print(percent_1000km_pref2,percent_1500km_pref2,percent_1750km_pref2,percent_2000km_pref2,percent_2500km_pref2,percent_3000km_pref2,percent_3500km_pref2,percent_4000km_pref2,percent_5000km_pref2) 

x = [500,1000,1500,1750,2000,2500,3000,3500,4000,5000]
y = [percent_500km_pref1,percent_1000km_pref1,percent_1500km_pref1,percent_1750km_pref1,percent_2000km_pref1,percent_2500km_pref1,percent_3000km_pref1,percent_3500km_pref1,percent_4000km_pref1,percent_5000km_pref1]
if numiter == 2:
    y2 = [percent_500km_pref2,percent_1000km_pref2,percent_1500km_pref2,percent_1750km_pref2,percent_2000km_pref2,percent_2500km_pref2,percent_3000km_pref2,percent_3500km_pref2,percent_4000km_pref2,percent_5000km_pref2]


percent_climo_60Ngen = np.array(percent_climo_60Ngen)
percent_climo_60Ngen_5per = np.array(percent_climo_60Ngen_5per)
percent_climo_60Ngen_95per = np.array(percent_climo_60Ngen_5per)

percent_climo_60N65 = np.array(percent_climo_60N65)
percent_climo_60N65_5per = np.array(percent_climo_60N65_5per)
percent_climo_60N65_95per = np.array(percent_climo_60N65_5per)

yticks = np.arange(0,101,10)

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
plt.plot(x,y,linestyle='-',linewidth='3',color='b',label=labeltext_pref1)
plt.plot(x,percent_climo_60Ngen,linestyle='--',linewidth='3',color='b')
if numiter==2:
    plt.plot(x,y2,linestyle='-',linewidth='3',color='r',label=labeltext_pref2)
    plt.plot(x,percent_climo_60N65,linestyle='--',linewidth='3',color='r')
    legend = ax1.legend(loc='center right', shadow=True,fontsize=label_fontsize)
    
    err1 = percent_climo_60Ngen + percent_climo_60Ngen_5per
    [einds] = np.where(err1>100)
    err1[einds] = 100.
    err2 = percent_climo_60Ngen - percent_climo_60Ngen_95per
    ax1.fill_between(x, err1, err2, alpha=0.7, edgecolor='0.6', facecolor='b')   
    err1 = percent_climo_60N65 + percent_climo_60N65_5per
    [einds] = np.where(err1>100)
    err1[einds] = 100.
    err2 = percent_climo_60N65 - percent_climo_60N65_95per
    ax1.fill_between(x, err1, err2, alpha=0.7, edgecolor='0.3', facecolor='r')    
ax1.set_ylim(-5,105)
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
plt.ylabel('Percentage of CAOs with nearby TPV',fontsize=label_fontsize)
plt.xlabel('Distance threshold from CAO centroid (km)',fontsize=label_fontsize) 
save_name = 'tpv_cao_NumberSensitivity.png'
plt.savefig(imagedir + save_name, bbox_inches='tight')  
plt.show()
