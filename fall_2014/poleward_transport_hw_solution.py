# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

# <codecell>

###########################
# User settings
###########################
datadir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2014/homework/'
infile1 = 'slp_sept_090.dat'
infile2 = 'slp_nov_090.dat'
label_fontsize = 14
# Get data at:
# http://www.esrl.noaa.gov/psd/cgi-bin/data/timeseries/timeseries1.pl
# Variable: Sea level pressure
# Analysis level: surface
# Latitude: Enter range; Longitude: 0-360
# Seasonal Average
# First month: August, Second month: November
# Area weight grids: No
# Raw data values
####################################

# <codecell>

###########################
# Read in data
###########################
datain1 = np.loadtxt(datadir+infile1, skiprows=0)
year1 = datain1[:,0]
value1 = datain1[:,1]

datain2 = np.loadtxt(datadir+infile2, skiprows=0)
year2 = datain2[:,0]
value2 = datain2[:,1]
###########################
# filter bad values out
###########################
# How do we know what's bad?  
#    Know what range of values you should expect and look in the file.
#    For sea level pressure, you should probably not get values below 800 hPa.
thresh = 800;
inds = np.argwhere(value1<thresh) 	
value1[inds] = float('NaN')    
inds = np.argwhere(value2<thresh) 	
value2[inds] = float('NaN') 

# <codecell>

meanslp_sept = value1[~np.isnan(value1)].mean()
meanslp_nov = value2[~np.isnan(value2)].mean()

dslp_nh = value2 - value1
mean_dslp = mean(dslp_nh[:-1])
# Take the mean.  The :-1 notation means to use all values except the very last entry.  
# We do this because we know it is NaN.
print 'Mean change in sea level pressure (SLP) from September to November is ' , mean_dslp, ' hPa'
print 'The mean change in SLP from September to November PER MONTH is ', mean_dslp/3, ' hPa'

# <codecell>

###########################
# Plot 1
###########################
t = year1

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,value1,'r',linewidth=3.0,label='September')
p2, = ax1.plot(t,value2,'b',linewidth=3.0,label='November')
ax1.set_xlim([np.min(year1),np.max(year1)])
ax1.set_ylim([1010,1016.5])
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.set_ylabel('Sea level pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Year',fontsize=label_fontsize)


plt.show()

# <codecell>

###########################
# Plot 2
###########################
y0 = np.zeros_like(dslp_nh).astype('f')
ym = np.zeros_like(dslp_nh).astype('f')
ym[:] = mean(dslp_nh[:-1])


golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,dslp_nh,'r',linewidth=3.0,label='Change in SLP')
p2, = ax1.plot(t,y0,'k',linewidth=3.0)
p1, = ax1.plot(t,ym,'r--',linewidth=3.0,label='Long-term mean change in SLP')
ax1.set_xlim([np.min(year1),np.max(year1)])
ax1.set_ylim([-6,6])
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.set_ylabel('Change in sea level pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Year',fontsize=label_fontsize)


plt.show()

# <codecell>


