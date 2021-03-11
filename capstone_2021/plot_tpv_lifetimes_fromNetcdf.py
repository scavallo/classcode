# plot_tpv_lifetimes_fromNetcdf
#
# Plots histograms and probability distribution functions for TPV track lifetimes.
# 
# Steven Cavallo
# March 2021
############################
# Imports
############################
import numpy.ma as ma
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

import warnings
warnings.filterwarnings("ignore")

import utilities_modules as um
import weather_modules as wm
from mstats import *
############################
# User directories
############################
fpath_sh = '/arctic3/datasets/tpv_SH/tracks/tracks_low_horizPlusVert.nc'
fpath_nh = '/arctic3/datasets/tpv/0.5x0.5/tracks_low.nc' # These are the northern hemisphere tracks 
imagedir = '/home/scavallo/scripts/python_scripts/images/' # Set this to a path on your local computer where you want to save your images, and uncomment out lines 276-277 below

###################################
# Other user options
###################################
hinc = 6 # number of hours between records
plot_option = 1 # 1 for histogram, 2 for probability density function
plot_variable = 1 # 1 to compare annual Northern vs. Southern hemisphere 
                  # 2 to compare winter vs. summer Southern hemisphere
		  # 3 to compare winter vs. summer Northern hemisphere 
		  # 4 to compare summer Northern vs. Southern Hemisphere
		  # 5 to compare winter Northern vs. Southern Hemisphere

min_lifetime_days = 2
label_fontsize = 18
npasses_smooth = 0

#############################
# Load data in
#############################
data_sh = netCDF4.Dataset(fpath_sh,'r') # you will need to change the path to point to your data
trackLen_sh = data_sh.variables['lenTrack'][:]
Dates_sh = data_sh.variables['timeStamp'][:]
mons_sh = np.array([int(str(x).replace('-','')[4:6]) for x in Dates_sh])
years_sh = np.array([int(str(x).replace('-','')[:4]) for x in Dates_sh])
Dates_sh = np.array([int(str(x).replace('-','')) for x in Dates_sh])
trackStartDate_sh = data_sh.variables['iTimeStart'][:]
trackStartMon_sh = mons_sh[trackStartDate_sh]
trackStartYear_sh = years_sh[trackStartDate_sh]
trackStartDate_sh = Dates_sh[trackStartDate_sh]


data_nh = netCDF4.Dataset(fpath_nh,'r') # you will need to change the path to point to your data
trackLen_nh = data_nh.variables['lenTrack'][:]
Dates_nh = data_nh.variables['timeStamp'][:]
mons_nh = np.array([int(str(x).replace('-','')[4:6]) for x in Dates_nh])
years_nh = np.array([int(str(x).replace('-','')[:4]) for x in Dates_nh])
Dates_nh = np.array([int(str(x).replace('-','')) for x in Dates_nh])
trackStartDate_nh = data_nh.variables['iTimeStart'][:]
trackStartMon_nh = mons_nh[trackStartDate_nh]
trackStartYear_nh = years_nh[trackStartDate_nh]
trackStartDate_nh = Dates_nh[trackStartDate_nh]

len_sh = len(trackLen_sh)
len_nh = len(trackLen_nh)
niters = np.max([len_sh,len_nh])
print(len_nh,len_sh,niters)

lifetimes_polar_nh = []
lifetimes_polar_sh = []
months_polar_nh = []
months_polar_sh = []
trackLen_polar_nh = []
trackLen_polar_sh = []

for x in range(0,niters):
#for x in range(0,60000):
  if x%10000 == 0:
    print ("On track {0}/{1}".format(x,niters))

  if x < len_nh:
  
      if trackLen_nh[x] < 8: # checking to make sure TPV track was longer than two days
          continue  
  
      lat_nh = data_nh.variables['latExtr'][x,:]
      if not ma.is_masked(lat_nh):
          per_life_in_arctic = float(np.where(lat_nh>=60)[0].shape[0])/float(lat_nh.shape[0]) # checking if TPV spent 60% of lifetime in Arctic
      else:
          per_life_in_arctic = float(np.where((lat_nh.data>=60)&(lat_nh.mask!=True))[0].shape[0])/float(np.where((lat_nh.mask!=True))[0].shape[0])
  
      if per_life_in_arctic*100 >= 60.:
          lifetimes_polar_nh.append(data_nh.variables['lenTrack'][x] / hinc)
	  months_polar_nh.append(trackStartMon_nh[x])
	  trackLen_polar_nh.append(data_nh.variables['lenTrack'][x])
  
  if x < len_sh:
      if trackLen_sh[x] < 8: # checking to make sure TPV track was longer than two days
          continue  
      
      lat_sh = data_sh.variables['latExtr'][x,:]
      if not ma.is_masked(lat_sh):   
          per_life_in_antarctic = float(np.where(lat_sh<=-60)[0].shape[0])/float(lat_sh.shape[0]) # checking if TPV spent 60% of lifetime in Antarctic
      else:
          per_life_in_antarctic = float(np.where((lat_sh.data<=-60)&(lat_sh.mask!=True))[0].shape[0])/float(np.where((lat_sh.mask!=True))[0].shape[0])
  
      if per_life_in_antarctic*100 >= 60.:
          lifetimes_polar_sh.append(data_sh.variables['lenTrack'][x] / hinc)
          months_polar_sh.append(trackStartMon_sh[x])
          trackLen_polar_sh.append(data_nh.variables['lenTrack'][x])

data_nh.close()
data_sh.close()

months_polar_nh = np.array(months_polar_nh)
months_polar_sh = np.array(months_polar_sh)
trackLen_polar_nh = np.array(trackLen_polar_nh)
trackLen_polar_sh = np.array(trackLen_polar_sh)

[winter_inds] = np.where( (months_polar_nh == 12) | (months_polar_nh <= 2) ) # only getting tpv tracks for winter months
trackLen_winter_nh = trackLen_polar_nh[winter_inds]
[summer_inds] = np.where( (months_polar_nh > 5) & (months_polar_nh < 9))
trackLen_summer_nh = trackLen_polar_nh[summer_inds]

del winter_inds, summer_inds
[summer_inds] = np.where( (months_polar_sh == 12) | (months_polar_sh <= 2) ) # only getting tpv tracks for winter months
trackLen_summer_sh = trackLen_polar_sh[summer_inds]
[winter_inds] = np.where( (months_polar_sh > 5) & (months_polar_sh < 9))
trackLen_winter_sh = trackLen_polar_sh[winter_inds]

lifetimes_nh = lifetimes_polar_nh
lifetimes_winter_nh = np.array(trackLen_winter_nh) / hinc
lifetimes_summer_nh = np.array(trackLen_summer_nh) / hinc

lifetimes_sh = lifetimes_polar_sh
lifetimes_winter_sh = np.array(trackLen_winter_sh) / hinc
lifetimes_summer_sh = np.array(trackLen_summer_sh) / hinc


#############################
# Get plot stuff ready
#############################
if plot_variable == 1:
    plotvar1 = lifetimes_nh
    plotvar2 = lifetimes_sh
    labeltext = ['Northern Hemisphere', 'Southern Hemisphere']
    save_name = 'tpv_lifetimes_annual_NH_vs_SH'   + '.png'
elif plot_variable == 2:
    plotvar1 = lifetimes_winter_sh
    plotvar2 = lifetimes_summer_sh    
    labeltext = ['Winter', 'Summer']
    save_name = 'tpv_lifetimes_winter_vs_summer_SH'   + '.png'
elif plot_variable == 3:
    plotvar1 = lifetimes_winter_nh
    plotvar2 = lifetimes_summer_nh       
    labeltext = ['Winter', 'Summer']
    save_name = 'tpv_lifetimes_winter_vs_summer_NH'   + '.png'
elif plot_variable == 4:
    plotvar1 = lifetimes_summer_nh
    plotvar2 = lifetimes_summer_sh       
    labeltext = ['Northern Hemisphere', 'Southern Hemisphere']
    save_name = 'tpv_lifetimes_summer_NH_vs_SH'   + '.png'
elif plot_variable == 5:
    plotvar1 = lifetimes_winter_nh
    plotvar2 = lifetimes_winter_sh       
    labeltext = ['Northern Hemisphere', 'Southern Hemisphere']
    save_name = 'tpv_lifetimes_winter_NH_vs_SH'   + '.png'

# Set the bin intervals (in days)
binvals = np.arange(0,81,1)

# These are used for the exponential curve fit on the plot
pplot1 = np.percentile(plotvar1, 99)
pplot2 = np.percentile(plotvar2, 99)

# These are just to print some statistics below
lq = np.percentile(plotvar1, 25)
uq = np.percentile(plotvar1, 75)
p90 = np.percentile(plotvar1, 90)
p95 = np.percentile(plotvar1, 95)
p99 = np.percentile(plotvar1, 99)
print "The 99th, 95th, 90th, and 75th lifetime percentiles are %.2f days, %.2f days, %.2f days, and %.2f days" % (p99, p95, p90, uq)

lq = np.percentile(plotvar2, 25)
uq = np.percentile(plotvar2, 75)
p90 = np.percentile(plotvar2, 90)
p95 = np.percentile(plotvar2, 95)
p99 = np.percentile(plotvar2, 99)
print "The 99th, 95th, 90th, and 75th lifetime percentiles are %.2f days, %.2f days, %.2f days, and %.2f days" % (p99, p95, p90, uq)


golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
############################################
fig = plt.figure(**figprops)   # New figure
############################################
ax1 = fig.add_subplot(1, 1, 1)
if plot_option == 1:
   alphaval = 0.75
   normval = 0
elif plot_option == 2:
   alphaval = 0.0
   normval = 1

n1, bins, patches = ax1.hist(plotvar1, bins=binvals, normed=normval, facecolor='blue', alpha=0, rwidth = 1.00) 

bincenters1 = 0.5*(bins[1:]+bins[:-1])
[p99inds] = np.where(bincenters1<=pplot1)
nn = n1[p99inds]
bincenters1b = bincenters1[p99inds]
nnn = np.log(nn)
inds = np.isfinite(nnn)    
m, b = np.polyfit(bincenters1b[inds],nnn[inds],1)               
plotvar_fit = m*bincenters1b + b    
exp_fit1 = np.exp(b)*np.exp(m*bincenters1b)
pearR = np.corrcoef(nnn[inds],exp_fit1[inds]); # Correlation coefficients
coefprint1 = pearR[0,1]

n2, bins, patches = ax1.hist(plotvar2, bins=binvals, normed=normval, facecolor='blue', alpha=0, rwidth = 1.00) 

bincenters2 = 0.5*(bins[1:]+bins[:-1])
del p99inds
[p99inds] = np.where(bincenters2<=pplot2)
nn = n2[p99inds]
bincenters2b = bincenters2[p99inds]
nnn = np.log(nn)
inds = np.isfinite(nnn)    
m, b = np.polyfit(bincenters2b[inds],nnn[inds],1)               
plotvar_fit = m*bincenters2b + b    
exp_fit2 = np.exp(b)*np.exp(m*bincenters2b)
pearR = np.corrcoef(nnn[inds],exp_fit2[inds]); # Correlation coefficients
coefprint2 = pearR[0,1]

n1[1:] = um.smooth_onedim(n1[1:],npasses_smooth)
n2[1:] = um.smooth_onedim(n2[1:],npasses_smooth)

plt.scatter(bincenters1,n1,color='b',s=50, marker='o', alpha=0.5, linewidths=None,label=labeltext[0])
plt.scatter(bincenters2,n2,color='r',s=50, marker='o', alpha=0.5, linewidths=None,label=labeltext[1])
colornow = 'b'
plt.plot(bincenters1b, exp_fit1, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint1))
colornow = 'r'
plt.plot(bincenters2b, exp_fit2, '-',linewidth=3,color=colornow,label="99%% exp. fit (r = %4.2f)"%(coefprint2))  
plt.yscale('log', nonposy='clip')
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True, fontsize=label_fontsize)
#plt.title(titletext1,fontsize=label_fontsize)
plt.xlim([binvals[0],binvals[-1]])
if plot_option == 1:
    plt.ylim([0,10000])

if plot_option == 1:
    plt.ylabel('Number',fontsize=label_fontsize)
    figname_suffix = 'histogram'
elif plot_option == 2:
    plt.ylabel('Probability',fontsize=label_fontsize)
    figname_suffix = 'pdf'

plt.xlabel('Lifetime (days)',fontsize=label_fontsize)
plt.savefig(imagedir + save_name, bbox_inches='tight')

plt.show()


