# plot_tpv_stats_pdfs_antarctica
#
# Plots histograms and probability distribution functions for TPV tracks in Southern Hemisphere.
# Will now also work for Northern Hemisphere.
# 
# Steven Cavallo
# February 2021
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
#from mpl_toolkits.basemap import Basemap, shiftgrid

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
# User options
############################
fpath_in = '/arctic3/datasets/tpv_SH/tracks/tracks_low_horizPlusVert.nc'
# fpath_in = '/arctic3/datasets/tpv/0.5x0.5/tracks_low.nc' # These are the northern hemisphere tracks 
imagedir = '' # Set this to a path on your local computer where you want to save your images, and uncomment out lines 276-277 below

# Variables in track file above: [can add to plot_variables below]
# circ, vortMean, ampMaxMin, rEquiv, thetaVol, ampMean, thetaExtr, latExtr, lonExtr

hemisphere = 'southern' # Either 'southern' or 'northern'
plot_type = 2 # 1 = histogram, 2 = pdf
plot_variable = 1 # 1 = minimum theta along track, 
		  # 2 = max amplitude along track, 
		  # 3 = median amplitude along track,
stat_sig_opt = 0 # 0 for K-S test, 
		 # 1 for difference in 2 samples based on confidence intervals set below, 
		 # 2 for bootstrap resampling based on confidence intervals set below
label_fontsize = 18
npasses_smooth = 1
legendtext = ['Annual','Winter','Summer']
confidence_interval_limits = [2.5,97.5]
#############################
# Load data in
#############################
data = netCDF4.Dataset(fpath_in,'r') # you will need to change the path to point to your data
trackLen = data.variables['lenTrack'][:]
Dates = data.variables['timeStamp'][:]
mons = np.array([int(str(x).replace('-','')[4:6]) for x in Dates])
years = np.array([int(str(x).replace('-','')[:4]) for x in Dates])
Dates = np.array([int(str(x).replace('-','')) for x in Dates])
trackStartDate = data.variables['iTimeStart'][:]
trackStartMon = mons[trackStartDate]
trackStartYear = years[trackStartDate]
trackStartDate = Dates[trackStartDate]

plotvar = []
plotvar_winter = []
plotvar_summer = []


if plot_type == 1:
   alphaval = 0.75
   normval = 0
elif plot_type == 2:
   alphaval = 0.0
   normval = 1

if plot_variable == 1:
   varname = 'thetaExtr'
   statmet = 'min'
   binvals = np.arange(250,320,1)
   binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
   bincenters = binvals[0:-1]
   xlabeltext = 'Kelvin'
   figname_prefix = 'tpv_mintheta_comparisons'
if plot_variable == 2:
   varname = 'ampMaxMin'
   statmet = 'max'
   binvals = np.arange(0,60,1)
   binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
   bincenters = binvals[0:-1]
   xlabeltext = 'Kelvin'
   figname_prefix = 'tpv_maxamplitude_comparisons'
if plot_variable == 3:
   varname = 'ampMaxMin'
   statmet = 'median'
   binvals = np.arange(0,60,1)
   binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
   bincenters = binvals[0:-1]
   xlabeltext = 'Kelvin'
   figname_prefix = 'tpv_medianamplitude_comparisons'

for x in range(0,len(trackLen)):
#for x in range(0,60000):
  if x%10000 == 0:
    print ("On track {0}/{1}".format(x,len(trackLen)))
  if trackLen[x] < 8: # checking to make sure TPV track was longer than two days
    continue

  lat = data.variables['latExtr'][x,:]
  if not ma.is_masked(lat):
     if hemisphere == 'northern':
        per_life_in_arctic = float(np.where(lat>=60)[0].shape[0])/float(lat.shape[0]) # checking if TPV spent 60% of lifetime in Arctic
     elif hemisphere == 'southern':
        per_life_in_arctic = float(np.where(lat<=-60)[0].shape[0])/float(lat.shape[0]) # checking if TPV spent 60% of lifetime in Antarctic
  else:
     if hemisphere == 'northern':
        per_life_in_arctic = float(np.where((lat.data>=60)&(lat.mask!=True))[0].shape[0])/float(np.where((lat.mask!=True))[0].shape[0])
     elif hemisphere == 'southern' :
        per_life_in_arctic = float(np.where((lat.data<=-60)&(lat.mask!=True))[0].shape[0])/float(np.where((lat.mask!=True))[0].shape[0])
  
  if per_life_in_arctic*100 < 60.:
    continue

  if statmet == 'min':
     plotvar.append(np.amin(data.variables[str(varname)][x,:]))
     if trackStartMon[x] == 12 or trackStartMon[x] <= 2: # only getting tpv tracks for winter months
       plotvar_summer.append(np.amin(data.variables[str(varname)][x,:]))
     elif trackStartMon[x] > 5 and trackStartMon[x] < 9:
       plotvar_winter.append(np.amin(data.variables[str(varname)][x,:]))# only getting tpv tracks for summer months
     else:
       continue
  if statmet == 'max':
     plotvar.append(np.amax(data.variables[str(varname)][x,:]))
     if trackStartMon[x] == 12 or trackStartMon[x] <= 2: # only getting tpv tracks for winter months
       plotvar_summer.append(np.amax(data.variables[str(varname)][x,:]))
     elif trackStartMon[x] > 5 and trackStartMon[x] < 9:
       plotvar_winter.append(np.amax(data.variables[str(varname)][x,:]))# only getting tpv tracks for summer months
     else:
       continue
  if statmet == 'median':
     plotvar.append(np.median(data.variables[str(varname)][x,:]))
     if trackStartMon[x] == 12 or trackStartMon[x] <= 2: # only getting tpv tracks for winter months
       plotvar_summer.append(np.median(data.variables[str(varname)][x,:]))
     elif trackStartMon[x] > 5 and trackStartMon[x] < 9:
       plotvar_winter.append(np.median(data.variables[str(varname)][x,:]))# only getting tpv tracks for summer months
     else:
       continue

data.close()

plotvar = np.array(plotvar)
plotvar_winter = np.array(plotvar_winter)
plotvar_summer = np.array(plotvar_summer)

mstats(plotvar)
mstats(plotvar_winter)
mstats(plotvar_summer)

############################
# Plots
############################
plotvar_cntl = plotvar
plotvar_exp1 = plotvar_winter
plotvar_exp2 = plotvar_summer

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

fig = plt.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

if plot_type == 1:
   n, bins, patches = ax1.hist(plotvar_cntl, binvals, normed=normval, facecolor='black', alpha=alphaval, rwidth = 1.00, align='mid',label=legendtext[0])
   n2, bins2, patches2 = ax1.hist(plotvar_exp1, binvals, normed=normval, facecolor='blue', alpha=alphaval, rwidth = 0.50, align='mid',label=legendtext[1])
   n3, bins3, patches3 = ax1.hist(plotvar_exp2, binvals, normed=normval, facecolor='red', alpha=alphaval, rwidth = 0.50, align='mid',label=legendtext[2]) 
    
elif plot_type == 2:
 
   n, bins  = np.histogram(plotvar_cntl, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
   n2, bins2 = np.histogram(plotvar_exp1, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)     
   n3, bins3 = np.histogram(plotvar_exp2, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
    
   n[1:] = um.smooth_onedim(n[1:],npasses_smooth)     
   n2[1:] = um.smooth_onedim(n2[1:],npasses_smooth)     
   n3[1:] = um.smooth_onedim(n3[1:],npasses_smooth) 
    
   ndiff = n2 - n
   nzeros = np.zeros_like(ndiff).astype('f')   

   if stat_sig_opt == 1:
      ul =  np.percentile(ndiff, confidence_interval_limits[1])
      ll =  np.percentile(ndiff, confidence_interval_limits[0])
      err1 = ndiff-ll
      err2 = ndiff-ul
   elif stat_sig_opt == 2:
      alpha = (confidence_interval_limits[0]*2.0)/100.
      boot_means,ci_limits_out = wm.bootstrap_twosamps(n,n2,10000,alpha)
      boot_means2,ci_limits_out2 = wm.bootstrap_twosamps(n,n3,10000,alpha)
      ci_limits = [100.0*alpha/2.0,100.0*(1.0-alpha/2.0)]
      print(ci_limits)

      ci_upper = np.percentile(boot_means,ci_limits[0])
      ci_lower = np.percentile(boot_means,ci_limits[1])
      print(ci_limits_out)
          
      err1a = n-np.abs(ci_limits_out[0])
      err1b = n+ci_limits_out[1]
        
      err2a = n2-np.abs(ci_limits_out[0])
      err2b = n2+ci_limits_out[1]
      
      err3a = n3-np.abs(ci_limits_out2[0])
      err3b = n3+ci_limits_out2[1]      

   if ( stat_sig_opt == 1 ):
      ax1.fill_between(bincenters, err1, err2, alpha=0.5, edgecolor='0.5', facecolor='0.5')
   elif ( stat_sig_opt == 2):
      ax1.fill_between(bincenters, err1a, err1b, alpha=0.5, edgecolor='0.5', facecolor='0.5')
      ax1.fill_between(bincenters, err2a, err2b, alpha=0.5, edgecolor='0.3', facecolor='0.3')	
      ax1.fill_between(bincenters, err3a, err3b, alpha=0.5, edgecolor='0.7', facecolor='0.7')	
   else: 
      p1, = ax1.plot(bincenters,n,'k',linewidth=3.0,label=legendtext[0])   
      p2, = ax1.plot(bincenters,n2,'0.3',linewidth=3.0,label=legendtext[1])   
      p3, = ax1.plot(bincenters,n3,'0.7',linewidth=3.0,label=legendtext[2])  
      
      ksstat1,pval1 = stats.ks_2samp(plotvar_cntl,plotvar_exp1)
      ksstat2,pval2 = stats.ks_2samp(plotvar_cntl,plotvar_exp2)
      print("K-S test between control and sample 1: ", ksstat1,pval1)
      print("K-S test between control and sample 2: ", ksstat2,pval2) 

ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
#plt.title(titletext1,fontsize=label_fontsize)

plt.xlim([binvals[0],binvals[-1]])
if plot_type == 1:
   plt.xlim([binvals[0],binvals[-1]]) 
   xticks_all = binvals[:]+0.5     
   xtickcint = 2
   xticks = xticks_all[::xtickcint]    
        
   xticklabels = binvals[::xtickcint]
      
if plot_type == 2:        
    
   xtickcint = 2
   if plot_variable == 7:    
      xtickcint = 100 
   if plot_variable == 12:
      xtickcint = 5
   plt.xlim([bincenters[0],bincenters[-1]])
   xticks = np.arange(binvals[0],binvals[-1],xtickcint)  
   xticklabels = np.arange(binvals[0],binvals[-1],xtickcint)  



ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticklabels, rotation = 90)

if plot_type == 1:
   plt.ylabel('Number',fontsize=label_fontsize)
   figname_suffix = 'histogram'
elif plot_type == 2:
   plt.ylabel('Probability',fontsize=label_fontsize)
   figname_suffix = 'pdf'
plt.xlabel(xlabeltext,fontsize=label_fontsize)
#save_name = figname_prefix + '_'  + figname_suffix + '.png'
#plt.savefig(imagedir + save_name, bbox_inches='tight')

plt.show()
