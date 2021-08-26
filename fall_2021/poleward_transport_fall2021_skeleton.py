#!/usr/bin/env python
# poleward_transport
#
# Steven Cavallo
# August 2021
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
from scipy import stats

import utilities_modules as um
from mstats import *

###########################
# User settings.  Students: set to your local paths!
###########################
datadir = ''
imagedir = ''
infile = 'slp_6590_annual_1948_2021.dat'

label_fontsize = 16
title_fontsize = 18
latitude = 60.

conversion_factor = 1000. # for converting meters to mm
analysis_years = [1979,2021]
plot_year = [2021]

ylims1 = [1006, 1022]
ylims2 =  [-6, 6]
ylims3 = [-2, 2]

month_labels = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
ndays = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]


###########################
# Constants
###########################
R_e = 6.37*10**6.0 # radius of Earth in meters
coslatr = np.cos(latitude*np.pi/180.) 
sinlatr = np.sin(latitude*np.pi/180.)
p0 = 100000.  # Pascals


###########################
# Read in the file and extract the subsets of the data that you need
###########################
ab = np.loadtxt(datadir+infile, skiprows=0)       
years = ab[:,0]
slp_full = ab[:,1:]

years = np.array(years)
slp_full = np.array(slp_full)
###########################
# Get the necessary indices
###########################
[ind1a] = np.where(years == analysis_years[0])
[ind1b] = np.where(years == analysis_years[1])
[ind1c] = np.where(years == plot_year[0])

ind1a = int(ind1a)
ind1b = int(ind1b)
ind1c = int(ind1c)

nyears,nmonths = np.shape(slp_full)

slp_full = um.filter_numeric_nans(slp_full,0,float('NaN'),'low')
slp_yearin = slp_full[ind1c,:]


###########################
# Calculate the monthly mean
###########################
slp_monthly_mean = np.nanmean(slp_full[ind1a:ind1b+1,:],0)

print(slp_monthly_mean)
print(nyears,nmonths)

######################
# mass transport calculations
######################
dslp_full = np.empty([nyears,nmonths])
dslp_yearin = np.empty(nmonths)
dslp_monthly_mean = np.empty(len(slp_monthly_mean))

# Take centered differences
dslp_full[:,0] = slp_full[:,1] - slp_full[:,0]
dslp_full[:,1:-1] = slp_full[:,2:] - slp_full[:,0:-2]
dslp_full[:,-1] = slp_full[:,-1] - slp_full[:,-2]

dslp_monthly_mean[0] = slp_monthly_mean[1] - slp_monthly_mean[0]
dslp_monthly_mean[1:-1] = slp_monthly_mean[2:] - slp_monthly_mean[0:-2]
dslp_monthly_mean[-1] = slp_monthly_mean[-1] - slp_monthly_mean[-2]
    
dslp_yearin[0] = slp_yearin[1] - slp_yearin[0]
dslp_yearin[1:-1] = slp_yearin[2:] - slp_yearin[0:-2]
dslp_yearin[-1] = slp_yearin[-1] - slp_yearin[-2]

# Calculate dp/dt
dpdt_monthly_mean = (dslp_monthly_mean*100.)/(np.array(ndays)*86400.) # Pa s-1    
dpdt_full = (dslp_full*100.)/(np.array(ndays)*86400.) # Pa s-1  
dpdt_yearin = (dslp_yearin*100.)/(np.array(ndays)*86400.) # Pa s-1  

h = R_e*(1-sinlatr)
vm_monthly_mean = h/(p0*coslatr)*(dpdt_monthly_mean)*conversion_factor
vm_full = h/(p0*coslatr)*(dpdt_full)*conversion_factor
vm_yearin = h/(p0*coslatr)*(dpdt_yearin)*conversion_factor

# For each year, what is the mean transport?
vm_yearly_means = h/(p0*coslatr)*(np.mean(dpdt_full[ind1a:ind1b+1,:],1))*conversion_factor
analysis_timeseries = np.arange(analysis_years[0],analysis_years[1]+1,1)

min_slp_yearly = np.min(slp_full[ind1a:ind1b+1,:],1)
max_slp_yearly = np.max(slp_full[ind1a:ind1b+1,:],1)
range_slp_yearly = np.ptp(slp_full[ind1a:ind1b+1,:],1)

# Pick a month for extracting an SLP timeseries
pickmonth_slp_yearly = slp_full[ind1a:ind1b+1,8]

# Pick the array to calculate a trend from
trend_plot_input = pickmonth_slp_yearly

# Calculate the trend.  Several ways to do this.  
# (1) Import the scipy stats module and use linregress OR
# (2) Use numpy polyfit with first order (linear) fit
#  slope, intercept = np.polyfit(analysis_timeseries[~np.isnan(trend_plot_input)],trend_plot_input[~np.isnan(trend_plot_input)],1)
#  trend_plot = slope*analysis_timeseries + intercept

slope, intercept, r_value, p_value, std_err = stats.linregress(analysis_timeseries[~np.isnan(trend_plot_input)],trend_plot_input[~np.isnan(trend_plot_input)])
trend_plot = slope*analysis_timeseries + intercept
print(p_value)


###########################
# Plot setup
###########################

tm = np.arange(1,13,1)
z0m = np.zeros(len(tm))

ty = np.arange(analysis_years[0],analysis_years[1]+1,1)
z0y = np.zeros(len(ty))

label_year_slp = 'Arctic SLP in ' + str(plot_year[0])

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 


###############################
# Figure 1
###############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(tm,slp_monthly_mean,'b',linewidth=3.0,label='Monthly mean Arctic SLP')
p2, = ax1.plot(tm,slp_yearin,'b--',linewidth=3.0,label=label_year_slp)
ax1.set_xlim([np.min(tm),np.max(tm)])
ax1.set_ylim([ylims1[0],ylims1[-1]])
ax1.xaxis.set_ticks(tm[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)

ax1.set_title('Monthly-mean sea level pressure over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('hPa',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

plt.savefig(imagedir + 'slp_65N_monthcomparisions' + '.png', bbox_inches='tight')
#plt.show()

###############################
# Figure 2: Students: Fill in the rest to make the plot
###############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)


plt.savefig(imagedir + 'vm_65N_monthcomparisions' + '.png', bbox_inches='tight')

#plt.show()

##############################
# Figure 3: Students: Fill in the rest to make the plot
###############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)


plt.savefig(imagedir + 'slp_65N_september_yearly' + '.png', bbox_inches='tight')

#plt.show()

