# poleward_transport
#
# Computes and plots the mass-averaged velocity
# Based on Wallace and Hobbs Problem 1.2.1
#
# Steven Cavallo
# August 2023
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
# User settings
###########################
datadir = '/' # Students: Need to put the path to where your data file is locally on your computer
imagedir = '/' # Students: Need to put the path to where you want the plots to save locally on your computer
# Example of mine below:
datadir = '/Users/scavallo/Documents/data/'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'

infile = 'slp_6590_annual_1948_2024.dat' # Name of data file

label_fontsize = 16 # Fontsize for plot labels
title_fontsize = 18 # Fontsize for plot titles
latitude = 65. # 

analysis_years = [1979,2024] # [start_year,end_year] of the calculations
plot_year = [2021] # Year for overlaying the comparisons
month_index_pick = 5 # Note that Python indices start with 0.  So 0 means January.
month_label_pick = 'June transport'

ylims1 = [1002,1026] # y-axis limits for monthly-mean SLP plot
# Students: Consider setting the axis limits for the other figures too
###########################
# Constants (i.e., do not change!)
###########################
month_labels = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
ndays = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
conversion_factor = 1000. # for converting meters to mm
R_e = 6.37*10**6.0 # radius of Earth in meters
coslatr = np.cos(latitude*np.pi/180.) 
sinlatr = np.sin(latitude*np.pi/180.)
p0 = 100000.  # Pascals

###########################
# Read in the file and extract the subsets of the data that you need
###########################
ab = np.loadtxt(datadir+infile, skiprows=0)       
years = ab[:,0]
slp_full_in = ab[:,1:]

years = np.array(years)

# Get the necessary indices
[ind1a] = np.where(years == analysis_years[0])
[ind1b] = np.where(years == analysis_years[1])
[ind1c] = np.where(years == plot_year[0])
[ind1d] = np.where(years == 2015)

# Convert the indices to integers so they can be referenced
ind1a = int(ind1a)
ind1b = int(ind1b)
ind1c = int(ind1c)
ind1d = int(ind1d)

# Pull in just the slp data that are necessary into numpy arrays
slp_full = np.array(slp_full_in[ind1a:ind1b+1,:])
slp_full = um.filter_numeric_nans(slp_full,0,float('NaN'),'low')
slp_yearin = slp_full_in[ind1c,:]
nyears,nmonths = np.shape(slp_full)


slp_19792010 = np.array(slp_full_in[ind1a:ind1d,:])
slp_20152024 = np.array(slp_full_in[ind1d:-1,:])

# Calculate the monthly mean and standard deviation and limits for the plots
slp_monthly_mean = np.nanmean(slp_full,0)
slp_monthly_std = np.nanstd(slp_full,0)

slp_monthly_mean = np.array(slp_monthly_mean)
slp_monthly_std = np.array(slp_monthly_std)

lowerlim_slp1 = slp_monthly_mean - slp_monthly_std
upperlim_slp1 = slp_monthly_mean + slp_monthly_std
lowerlim_slp2 = np.nanmin(slp_full,0)
upperlim_slp2 = np.nanmax(slp_full,0)

# Compute monthly-mean SLP for the subsets of SLP data
slp_monthly_mean_19792010 = np.nanmean(slp_19792010,0)
slp_monthly_mean_20152024 = np.nanmean(slp_20152024,0)

######################
# mass transport calculations
######################
dslp_full = np.empty([nyears,nmonths])
dslp_yearin = np.empty(nmonths)
dslp_monthly_mean = np.empty(len(slp_monthly_mean))

# Calculate the changes in SLP using centered differences
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

# Calculate mass-averaged velocity across the given latitude using the equation derived on the homework
h = R_e*(1-sinlatr)
vm_monthly_mean = h/(p0*coslatr)*(dpdt_monthly_mean)*conversion_factor
vm_full = h/(p0*coslatr)*(dpdt_full)*conversion_factor
vm_yearin = h/(p0*coslatr)*(dpdt_yearin)*conversion_factor

# Compute the monthly standard deviations so that variability can be displayed on the plots
vm_monthly_std = np.nanstd(vm_full,0)
lowerlim_vm1 = vm_monthly_mean - vm_monthly_std
upperlim_vm1 = vm_monthly_mean + vm_monthly_std
lowerlim_vm2 = np.nanmin(vm_full,0)
upperlim_vm2 = np.nanmax(vm_full,0)

# Need to make a custom vector for x to pass into stats functions below
analysis_timeseries = np.arange(analysis_years[0],analysis_years[1]+1,1)

min_slp_yearly = np.min(slp_full,1)
max_slp_yearly = np.max(slp_full,1)
range_slp_yearly = np.ptp(slp_full,1)

# Pick a month for extracting timeseries
pickmonth_slp_yearly = slp_full[:,month_index_pick]
pickmonth_vm_yearly = vm_full[:,month_index_pick]


# Pick the array to calculate a trend from
#trend_plot_input = pickmonth_slp_yearly
trend_plot_input = pickmonth_vm_yearly
if trend_plot_input[-1] == -999:
    trend_plot_input[-1] = np.float('NaN')


# Calculate the trend.  Several ways to do this.  
# (1) Import the scipy stats module and use linregress OR
# (2) Use numpy polyfit with first order (linear) fit
#  slope, intercept = np.polyfit(analysis_timeseries[~np.isnan(trend_plot_input)],trend_plot_input[~np.isnan(trend_plot_input)],1)
#  trend_plot = slope*analysis_timeseries + intercept

slope, intercept, r_value, p_value, std_err = stats.linregress(analysis_timeseries[~np.isnan(trend_plot_input)],trend_plot_input[~np.isnan(trend_plot_input)])
trend_plot = slope*analysis_timeseries + intercept
print("The p-value for the trend in SLP is %7.4f" %p_value)

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
p1 = plt.fill_between(tm,lowerlim_slp1,upperlim_slp1,alpha=0.4,facecolor='0.3',edgecolor='0.3',label=r'$\pm$ 1 standard deviation')
p2 = plt.fill_between(tm,lowerlim_slp2,upperlim_slp2,alpha=0.4,facecolor='0.6',edgecolor='0.6',label='Range')
p3, = ax1.plot(tm,slp_monthly_mean,'b',linewidth=4.0,label='Monthly mean Arctic SLP')
p4, = ax1.plot(tm,slp_monthly_mean_19792010,'r',linewidth=3.0,label='1979-2010')
p5, = ax1.plot(tm,slp_monthly_mean_20152024,'r--',linewidth=3.0,label='2015-2024')

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

###############################
# Figure 2
###############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

# Students: Fill in the rest for the monthly mean mass average velocity across 65 deg N latitude.  
# Hint: You can cut/copy/paste much of the above code but change the variable names and labels to the appropriate ones


##############################
# Figure 3
###############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

# Students: Fill in the rest for the February monthly-mean v_m as a function of year (1979 to 2023)

plt.show()


