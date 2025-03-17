#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt


# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from scipy import ndimage
import scipy.io as sio
from scipy import stats
from scipy.stats import t
import statsmodels.api as sm
import warnings
warnings.filterwarnings("ignore")



###################################
# Set user options
###################################
date_firstrecord = '1979010100' # only if record_num > -1
date_lastrecord = '2018123118' # only if record_num > -1
map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
zoom = 'true'
proj_latlon = [90. , 270.]
hinc = 6 # number of hours between records;  only used if record_num > -1

plot_option = 1 # 1 for histogram, 2 for probability density function
num_files = 4
overlay_genesis_lysis = 'False'
make_tpv = 'False'
min_lifetime_days = 2
lat_origin_range = [30., 90.]
#lat_origin_range = [-90.,-30.]
label_fontsize = 12
npasses_smooth = 0

#fdir2 = '/Users/scavallo/Documents/data/track_files/southern_hemisphere/'
#fdir1 = '/Users/scavallo/Documents/data/track_files/'
#textfile = '2013_2014_trth_tracks_all.txt'
#textfile = 'trth_tracks_all_2006060100_2006083118.dat'
#textfile = 'trth_tracks_all_1979010100_2014123118.dat'
#textfile = 'trth_tracks_03_1979010100_2014123118.dat'
#textfile1 = 'trth_tracks_all_1979010100_2014123118.dat'
#textfile1 = 'trth_tracks_djf_1979010100_2014123118.dat'
#textfile2 = 'trth_tracks_jja_1979010100_2014123118.dat'

#textfile2 = 'trth_tracks_SH_60S65_1979010100_2018123118.dat'
#textfile2 = 'trth_tracks_SH_60S65_son_1979010100_2018123118.dat'
#textfile1 = 'trth_tracks_60N65_1979010100_2018123118.dat'
#textfile2 = 'trth_tracks_60N65_jja_1979010100_2018123118.dat'

fdir1 = '/Users/scavallo/Documents/data/sfc_cyclones/ERA5/sprenger/'
fdir2 = '/Users/scavallo/Documents/Work/Research/Papers/Arctic_TPV_interactions/data/2024_60Ngen/radius_933km/'#intensity_10percentile/'
fdir3 = '/Users/scavallo/Documents/Work/Research/Papers/Arctic_TPV_interactions/data/2024_60Ngen/radius_933km/'#intensity_10percentile/'
fdir4 = '/Users/scavallo/Documents/Work/Research/Papers/Arctic_TPV_interactions/data/2024_60Ngen/radius_933km/'#intensity_10percentile/'
textfile1 = 'ERA5_sfccyclones_arctictracks_2days_synoptic_1979010100_2019123100.dat' # Enter name of track file 1
textfile2 = 'ERA5_sfccyclones_arctictracks_1ormoreinteractions_synoptic_60N65_1979010100_2021123100.dat'
textfile3 = 'ERA5_sfccyclones_arctictracks_2ormoreinteractions_synoptic_60N65_1979010100_2021123100.dat'
textfile4 = 'ERA5_sfccyclones_arctictracks_3ormoreinteractions_synoptic_60N65_1979010100_2021123100.dat'
#textfile2 = 'ERA5_sfccyclone_NH_tpvlink_1500km_1ormoreinteractions_1979010100_2018123118.dat'
#textfile3 = 'ERA5_sfccyclone_NH_tpvlink_1500km_2ormoreinteractions_1979010100_2018123118.dat' # Enter name of track file 2
#textfile4 = 'ERA5_sfccyclone_NH_tpvlink_1500km_5ormoreinteractions_1979010100_2018123118.dat'

labeltext = ['All Arctic Cyclones', '1 or more', '2 or more', '3 or more']
hemispheres = ['northern','northern','northern','northern','northern']

#imagedir = '/data1/scavallo/data/cases/decjan_2013/for_colucci/'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'

###################################
# END user options
###################################
R_earth = 6371200
R_earth = R_earth / 1000
pid = np.pi/180


if plot_option == 1:
    density_opt = False
elif plot_option == 2:
    density_opt = True

for iter in range(0,num_files):
    hemnow = hemispheres[iter]
    if hemnow == 'northern':
        lat_origin_range = [30., 90.]
    elif hemnow == 'southern':
        lat_origin_range = [-90.,-30.]    
    if iter==0:
        textfile = textfile1
        fdir = fdir1
    elif iter == 1:
        textfile = textfile2
        fdir = fdir2
    elif iter==2:
        textfile = textfile3
        fdir = fdir3
    elif iter == 3:
        textfile = textfile4
        fdir = fdir4


        	
    fpath = fdir + textfile
    print(fpath)
    aa = np.loadtxt(fpath, skiprows=1)       
    datelist = aa[:,0]
    lat = aa[:,1]	
    lon = aa[:,2]
    thetamin = aa[:,3]	    
    thetaamp = aa[:,4]
    vort_circ = aa[:,5]	    
    vort_radius = aa[:,6]

    nrows, ncols = np.shape(aa)
    print(nrows, ncols)

    tinds = np.where(datelist>10000)
    lats_all = lat[tinds]
    lons_all = lon[tinds]
    Y = lats_all[0:-1]
    X = lons_all[0:-1]
    lat2 = lats_all[1:]
    lon2 = lons_all[1:]

    lifetime = []

    ntracks_min = (24./hinc)*min_lifetime_days
    tt = 0
    tracknum = -1
    while tt < nrows:
        if (datelist[tt] < 10000) :
		
            ntracks = datelist[tt]
            nhours = ntracks*hinc
            sind = tt+1
            eind = int(sind+ntracks)

            if make_tpv == 'True':
                latnow = lat[sind:eind]
                nlats_polar = np.where(latnow>=65.0)
                perc = (np.float(np.size(nlats_polar)) / np.float(np.size(latnow)))*100.
                lat_origin_range = [0.,90.]
            if min_lifetime_days < 2:
                min_lifetime_days = 2.0
                ntracks_min = (24./hinc)*min_lifetime_days
                perc_thresh = 60.
                #print  perc
            else:
                perc = 0.0    	
                perc_thresh = -1.0        

			#if ( (lat[sind] >= lat_origin_thresh) and (ntracks >= ntracks_min) ):
            if ( (lat[sind] >= lat_origin_range[0]) and (lat[sind] <= lat_origin_range[1]) and (ntracks >= ntracks_min) and (perc >= perc_thresh) ):
                tracknum += 1

                latnow = lat[sind:eind]
                lonnow = lon[sind:eind]
                trthnow = thetamin[sind:eind]
                ampsnow = thetaamp[sind:eind]	    
                radsnow = vort_radius[sind:eind]	    
                circsnow = vort_circ[sind:eind]	    

                lifetime.append(nhours)

                #distnow = um.earth_distm(latnow[0],lonnow[0],np.array(latnow[0:-2]),np.array(lonnow[1:-1])) 
                Y = latnow[0:-1]
                X = lonnow[0:-1]
                lat2 = latnow[1:]
                lon2 = lonnow[1:]
	

        tt += 1
		
    if iter == 0:
        lifetime1 = lifetime	    
    elif iter == 1:
        lifetime2 = lifetime	 
    elif iter == 2:
        lifetime3 = lifetime
    elif iter == 3:
        lifetime4 = lifetime        	
    elif iter == 4:
        lifetime5 = lifetime		

titletext1 = 'Vortex Latitude'
titletext2 = 'Vortex lifetime'
titletext3 = 'Vortex amplitude'
titletext4 = 'Median vortex amplitude'

#lifetime_days = np.array(lifetime2) / 24.0
lifetime_days1 = np.array(lifetime1) / 24.0
if num_files >= 2:
    lifetime_days2 = np.array(lifetime2) / 24.0
if num_files >= 3:
    lifetime_days3 = np.array(lifetime3) / 24.0
if num_files >= 4:
    lifetime_days4 = np.array(lifetime4) / 24.0

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

#binvals = np.arange(0,96,1)
binvals = np.arange(0,26,1)


lifetime_days = lifetime_days1
mstats(lifetime_days)
#lq = np.percentile(lifetime_days, 25)
uq = np.percentile(lifetime_days, 75)
p90 = np.percentile(lifetime_days, 90)
p95 = np.percentile(lifetime_days, 95)
p99 = np.percentile(lifetime_days, 99)

fit_percentage = 99
p99_1 = np.percentile(lifetime_days1, fit_percentage)
if num_files >= 2:
    p99_2 = np.percentile(lifetime_days2, fit_percentage)
if num_files >= 3:
    p99_3 = np.percentile(lifetime_days3, fit_percentage)
if num_files >= 4:
    p99_4 = np.percentile(lifetime_days4, fit_percentage)
#p99p9 = np.percentile(lifetime_days, 99.9)
#print "The upper and lower lifetime quartiles are %.2f days and %.2f days" % (uq, lq)
print("The 99th, 95th, 90th, and 75th lifetime percentiles are %.2f days, %.2f days, %.2f days, and %.2f days" % (p99, p95, p90, uq))

slopes_arr = []
std_err_arr = []
df_arr = []

#############################
fig = plt.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

if plot_option == 1:
   alphaval = 0.75
   normval = 0
elif plot_option == 2:
   alphaval = 0.0
   normval = 1


n1, bins, patches = ax1.hist(lifetime_days1, bins=binvals, density=density_opt, histtype='bar',color='b', alpha=0, rwidth = 1.00) 
bincenters1 = 0.5*(bins[1:]+bins[:-1])


[p99inds] = np.where(bincenters1<=p99_1)
nn = n1[p99inds]
bincenters1b = bincenters1[p99inds]
nnn = np.log(nn)
inds = np.isfinite(nnn)    
m, b = np.polyfit(bincenters1b[inds],nnn[inds],1)               
plotvar_fit = m*bincenters1b + b    
exp_fit1 = np.exp(b)*np.exp(m*bincenters1b)
pearR = np.corrcoef(nnn[inds],exp_fit1[inds]); # Correlation coefficients
coefprint1 = pearR[0,1]

n1[1:] = um.smooth_onedim(n1[1:],npasses_smooth)

xin = bincenters1b[inds]
yin = nnn[inds]
df = len(yin)-1
slope, intercept, r_value, p_value, std_err = stats.linregress(xin,yin)
slopes_arr.append(slope)
std_err_arr.append(std_err)
df_arr.append(df)    

if num_files >= 2:
    n2, bins, patches = ax1.hist(lifetime_days2, bins=binvals, density=density_opt, histtype='bar',color='b', alpha=0, rwidth = 1.00) 
    bincenters2 = 0.5*(bins[1:]+bins[:-1])

    del p99inds
    [p99inds] = np.where(bincenters2<=p99_2)
    nn = n2[p99inds]
    bincenters2b = bincenters2[p99inds]
    nnn = np.log(nn)
    inds = np.isfinite(nnn)    
    m, b = np.polyfit(bincenters2b[inds],nnn[inds],1)               
    plotvar_fit = m*bincenters2b + b    
    exp_fit2 = np.exp(b)*np.exp(m*bincenters2b)
    pearR = np.corrcoef(nnn[inds],exp_fit2[inds]); # Correlation coefficients
    coefprint2 = pearR[0,1]

    n2[1:] = um.smooth_onedim(n2[1:],npasses_smooth)
    
    xin = bincenters2b[inds]
    yin = nnn[inds]
    df = len(yin)-1
    slope, intercept, r_value, p_value, std_err = stats.linregress(xin,yin)
    slopes_arr.append(slope)
    std_err_arr.append(std_err)
    df_arr.append(df)    
if num_files >= 3:
    n3, bins, patches = ax1.hist(lifetime_days3, bins=binvals, density=density_opt, histtype='bar',color='b', alpha=0, rwidth = 1.00) 
    bincenters3 = 0.5*(bins[1:]+bins[:-1])

    del p99inds
    [p99inds] = np.where(bincenters3<=p99_3)
    nn = n3[p99inds]
    bincenters3b = bincenters3[p99inds]
    nnn = np.log(nn)
    inds = np.isfinite(nnn)    
    m, b = np.polyfit(bincenters3b[inds],nnn[inds],1)               
    plotvar_fit = m*bincenters3b + b    
    exp_fit3 = np.exp(b)*np.exp(m*bincenters3b)
    pearR = np.corrcoef(nnn[inds],exp_fit3[inds]); # Correlation coefficients
    coefprint3 = pearR[0,1]

    n3[1:] = um.smooth_onedim(n3[1:],npasses_smooth)    

    xin = bincenters3b[inds]
    yin = nnn[inds]
    df = len(yin)-1
    slope, intercept, r_value, p_value, std_err = stats.linregress(xin,yin)
    slopes_arr.append(slope)
    std_err_arr.append(std_err)
    df_arr.append(df)    
if num_files >= 4:
    n4, bins, patches = ax1.hist(lifetime_days4, bins=binvals, density=density_opt, histtype='bar',color='b', alpha=0, rwidth = 1.00) 
    bincenters4 = 0.5*(bins[1:]+bins[:-1])

    del p99inds
    [p99inds] = np.where(bincenters4<=p99_4)
    nn = n4[p99inds]
    bincenters4b = bincenters4[p99inds]
    nnn = np.log(nn)
    inds = np.isfinite(nnn)    
    m, b = np.polyfit(bincenters4b[inds],nnn[inds],1)               
    plotvar_fit = m*bincenters4b + b    
    exp_fit4 = np.exp(b)*np.exp(m*bincenters4b)
    pearR = np.corrcoef(nnn[inds],exp_fit4[inds]); # Correlation coefficients
    coefprint4 = pearR[0,1]

    n4[1:] = um.smooth_onedim(n4[1:],npasses_smooth)

    xin = bincenters4b[inds]
    yin = nnn[inds]
    df = len(yin)-1
    slope, intercept, r_value, p_value, std_err = stats.linregress(xin,yin)
    slopes_arr.append(slope)
    std_err_arr.append(std_err)
    df_arr.append(df)    

p95line = np.zeros_like(n1).astype('float')
p95line[:] = p95

# t-test for statistical significance between slopes.  The statistical significance test goes as follows:
# (1) H0 (null hypothesis): slope 1 = slope 2
#     H1 (alternative hypothesis): slope 1 does not equal slope 2
# (2) Use the t-statistic
# (3) Choose desired confidence interval.  Usually 95 percent, so let alpha = 0.025.
# (4) Find the t_{0.025} value in a t-table using 'degfree' below as the number of degrees of freedom
#     https://www.sjsu.edu/faculty/gerstman/StatPrimer/t-table.pdf
# (5) To accept the null hypothesis, the following must be satisfied:
#              -t_{0.025} < t < t_{0.025}
#     where t is computed below.
# (6) If the null hypothesis is rejected, then the slopes are significantly different for the chosen confidence interval.
slopes_arr = np.array(slopes_arr)
std_err_arr = np.array(std_err_arr)
tstats = (slopes_arr[0::-1] - slopes_arr[1::])/(np.sqrt(std_err_arr[0::-1]**2 +std_err_arr[1::]**2));
#dfstat = np.min(df_arr)
print("tstats are ", tstats)
print("degrees of freedom are ", df_arr)

    
    
#n2[1:] = um.smooth_onedim(n2[1:],npasses_smooth)
#n = um.zero2nan(n)   
#p1, = ax1.plot(bincenters1b,n1,'0.3',linewidth=3.0,label='DJF')     
#p2, = ax1.plot(bincenters2b,n2,'0.7',linewidth=3.0,label='JJA') 
plt.scatter(bincenters1,n1,color='b',s=20, marker='o', alpha=0.5, linewidths=None,label=labeltext[0])
if num_files >= 2:
    plt.scatter(bincenters2,n2,color='r',s=50, marker='o', alpha=0.5, linewidths=None,label=labeltext[1])
if num_files >= 3:
    plt.scatter(bincenters3,n3,color='g',s=80, marker='o', alpha=0.5, linewidths=None,label=labeltext[2])
if num_files >= 4:
    plt.scatter(bincenters4,n4,color='m',s=110, marker='o', alpha=0.5, linewidths=None,label=labeltext[3])

colornow = 'b'
plt.plot(bincenters1b, exp_fit1, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint1))
if num_files >= 2:
    colornow = 'r'
    plt.plot(bincenters2b, exp_fit2, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint2))  
if num_files >= 3:
    colornow = 'g'
    plt.plot(bincenters3b, exp_fit3, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint3))
if num_files >= 4:
    colornow = 'm'
    plt.plot(bincenters4b, exp_fit4, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint3))

#plt.axvline(x=p95,color='k',linewidth=8,label='95th percentile')
plt.yscale('log', nonpositive='clip')
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True, fontsize=label_fontsize)
#plt.title(titletext1,fontsize=label_fontsize)
plt.xlim([binvals[0],binvals[-1]])
ylims = ax1.get_ylim()
if plot_option == 1:
    plt.ylim([0.5,ylims[1]])

if plot_option == 1:
    plt.ylabel('Number',fontsize=label_fontsize)
elif plot_option == 2:
    plt.ylabel('Probability',fontsize=label_fontsize)
figname_suffix = 'histogram'
#elif plot_option == 2:
    #plt.ylabel('Probability',fontsize=label_fontsize)
    #figname_suffix = 'pdf'

plt.xlabel('Lifetime (days)',fontsize=label_fontsize)

um.bold_labels(ax1,label_fontsize)
save_name = 'compare_tpv_lifetimes_'  + figname_suffix + '_' + date_firstrecord +  '_' + date_lastrecord + '.png'
plt.savefig(imagedir + save_name, bbox_inches='tight')

plt.show()
#exit()
