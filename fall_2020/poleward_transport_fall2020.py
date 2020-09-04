#!/usr/bin/env python
#
# poleward_transport_fall2020
# Steven Cavallo
# August 2020
#
# Version info
# >>conda -V
# conda 4.3.21
# >>python --version
# Python 3.5.2 :: Anaconda custom (x86_64)


###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

# Make sure you link to these in your working directory!  They are available in the 'utils' directory on github:
# https://github.com/scavallo/classcode/tree/master/utils
import utilities_modules as um
from mstats import *

###########################
# User settings
###########################
datadir = '/Users/scavallo/Documents/data/'
#imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2018/homework/images/'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
infile = 'slp_6590_annual_1948_2020.dat'

label_fontsize = 16
title_fontsize = 18
latitude = 60.

analysis_years = [1979,2020]
plot_year = [2020]
plot_year2 = [2019]

ylims1 = [1004, 1018]
ylims2 =  [-6, 6]

month_labels = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
ndays = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]
###########################
# Constants
###########################
R_e = 6.37*10**6.0 # radius of Earth in meters
coslatr = np.cos(latitude*np.pi/180.) 
sinlatr = np.sin(latitude*np.pi/180.)
p0 = 100000.  # Pascals
conversion_factor = 1000. # for converting meters to mm

###########################
# Read in the file and extract the subsets of the data that you need
###########################
ab = np.loadtxt(datadir+infile, skiprows=0)       
years = ab[:,0]
slp_full = ab[:,1:]

###########################
# Get the necessary indices
###########################
[ind1a] = np.where(years == analysis_years[0])
[ind1b] = np.where(years == analysis_years[1])
[ind1c] = np.where(years == plot_year[0])
[ind1d] = np.where(years == plot_year2[0])

ind1a = np.int(ind1a)
ind1b = np.int(ind1b)
ind1c = np.int(ind1c)
ind1d = np.int(ind1d)

nyears,nmonths = np.shape(slp_full)

slp_full = um.filter_numeric_nans(slp_full,0,float('NaN'),'low')
slp_yearin = slp_full[ind1c,:]
slp_yearin2 = slp_full[ind1d,:]

slp_monthly_mean = np.nanmean(slp_full[ind1a:ind1b+1],0)
print(nyears,nmonths)

######################
# mass transport calculations
######################
dslp_full = np.zeros([nyears,nmonths])
dslp_yearin = np.zeros(nmonths)
dslp_yearin2 = np.zeros(nmonths)
dslp_monthly_mean = np.zeros(len(slp_monthly_mean))

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

dslp_yearin2[0] = slp_yearin2[1] - slp_yearin2[0]
dslp_yearin2[1:-1] = slp_yearin2[2:] - slp_yearin2[0:-2]
dslp_yearin2[-1] = slp_yearin2[-1] - slp_yearin2[-2]


# Calculate dp/dt
dpdt_monthly_mean = (dslp_monthly_mean*100.)/(np.array(ndays)*86400.) # Pa s-1    
dpdt_full = (dslp_full*100.)/(np.array(ndays)*86400.) # Pa s-1  
dpdt_yearin = (dslp_yearin*100.)/(np.array(ndays)*86400.) # Pa s-1  
dpdt_yearin2 = (dslp_yearin2*100.)/(np.array(ndays)*86400.) # Pa s-1  

# mass transport
h = R_e*(1-sinlatr)
vm_monthly_mean = h/(p0*coslatr)*(dpdt_monthly_mean)*conversion_factor
vm_full = h/(p0*coslatr)*(dpdt_full)*conversion_factor
vm_yearin = h/(p0*coslatr)*(dpdt_yearin)*conversion_factor
vm_yearin2 = h/(p0*coslatr)*(dpdt_yearin2)*conversion_factor
######################
# PLot setup
######################
tm = np.arange(1,13,1)
z0m = np.zeros(len(tm))

ty = np.arange(analysis_years[0],analysis_years[1]+1,1)
z0y = np.zeros(len(ty))

label_year_slp = 'Arctic SLP in ' + str(plot_year[0])
label_year_slp2 = 'Arctic SLP in ' + str(plot_year2[0])

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
#p3, = ax1.plot(tm,slp_yearin2,'r--',linewidth=3.0,label=label_year_slp2)
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

p1, = ax1.plot(tm,vm_monthly_mean,'b',linewidth=3.0,label='Monthly mean')
p2, = ax1.plot(tm,vm_yearin,'b--',linewidth=3.0,label=str(plot_year[0]))
#p3, = ax1.plot(tm,vm_yearin2,'r--',linewidth=3.0,label=str(plot_year2[0]))
ax1.set_xlim([np.min(tm),np.max(tm)])
ax1.xaxis.set_ticks(tm[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower right', shadow=True)
p0, = ax1.plot(tm,np.zeros(len(tm)),'k',linewidth=2.0)

ax1.set_title('Mass-weighted velocity transport over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('mm s-1',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)

plt.savefig(imagedir + 'vm_65N_monthcomparisions' + '.png', bbox_inches='tight')
plt.show()

