# Steven Cavallo
# METR 5004 Homework 1 plotting code
# September 1, 2015


###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
from mstats import *

###########################
# User settings
###########################
datadir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2015/homework/'
imagedir = '/Users/scavallo/scripts/python_scripts/images/'
infile = 'slp_monthly_6090.dat'
label_fontsize = 14
# Get data at:
# http://www.esrl.noaa.gov/psd/cgi-bin/data/timeseries/timeseries1.pl
# Variable: Sea level pressure
# Analysis level: surface
# Latitude: 60-90; Longitude: 0-360
# Monthly
# Area weight grids: No
# Raw data values
# Save to a text file
####################################


###########################
# Read in data
###########################
datain1 = np.loadtxt(datadir+infile, skiprows=0)
years = datain1[:,0] # Read in all rows of the first column
# Hereafter, I am always skipping the first column because it is the year
vals = datain1[:,1:]
apr = datain1[:,4] # year is index 0, january 1, february 2, march 3, april 4, ...
july = datain1[:,7]

###########################
# filter bad values out
###########################
# How do we know what's bad?  
#    Know what range of values you should expect and look in the file.
#    For sea level pressure, you should probably not get values below 800 hPa.
thresh = 800;
inds = np.argwhere(vals<thresh) 	
vals[inds] = float('NaN')   

slp_monthly_mean = np.nanmean(vals,0)

[maxind] = np.where(slp_monthly_mean == max(slp_monthly_mean))
[minind] = np.where(slp_monthly_mean == min(slp_monthly_mean))
print "The maximum SLP is %2.1f in month %d" % (slp_monthly_mean[maxind], maxind)
print "The minimum SLP is %2.1f in month %d" % (slp_monthly_mean[minind], minind)

slp_maxminusmin = slp_monthly_mean[maxind] - slp_monthly_mean[minind]

inds = np.argwhere(apr<thresh) 	
apr[inds] = float('NaN')    


inds = np.argwhere(july<thresh) 	
july[inds] = float('NaN')


slp_diff_maxmin = apr - july
mean_slpdiff = np.mean(slp_diff_maxmin[:-1])


# declare an array of zeros
dslp = np.zeros_like(slp_monthly_mean)

dslp[0:-1] = slp_monthly_mean[1:] - slp_monthly_mean[0:-1]
dslp[-1] = slp_monthly_mean[0] - slp_monthly_mean[-1]

ym1 = np.zeros_like(slp_monthly_mean).astype('f')
ym2 = np.zeros_like(slp_diff_maxmin).astype('f')
ym1[:] = slp_maxminusmin
ym2[:] = mean_slpdiff


###########################
# Plot 1
###########################
#t = year1
t = np.arange(1,13,1)

#numticks = 9
#yticks1 = np.linspace(1010,1018,num=numticks)
#yticks2 = np.linspace(-4,4,num=numticks)

cint = 2
yticks2 = np.arange(-10,10+(cint/2),cint)
numticks = len(yticks2)
yticks1 = np.linspace(1010,1018,num=numticks)

y0 = np.zeros_like(dslp)

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,slp_monthly_mean,'b',linewidth=3.0,label='Sea level pressure')

ax1.set_ylim([1010,1016.5])
ax1.grid(True, linestyle='-')

ax2 = ax1.twinx()
p2, = ax2.plot(t,dslp,'r',linewidth=3.0,label='Change in sea level pressure')
p1, = ax2.plot(t,ym1,'g--',linewidth=3.0,label='Long-term mean difference in SLP')
p0, = ax2.plot(t,y0,'k',linewidth=3.0)
ax2.grid(True, linestyle='--')
ax2.set_ylim([-4,4])

ax1.yaxis.set_ticks(yticks1)
ax1.set_yticklabels(yticks1)
ax2.yaxis.set_ticks(yticks2)
ax2.set_yticklabels(yticks2)

ax1.set_xlim([np.min(t),np.max(t)])
ax2.set_xlim([np.min(t),np.max(t)])

legend = ax1.legend(loc='lower left', shadow=True)
legend = ax2.legend(loc='upper left', shadow=True)
ax1.set_ylabel('Sea level pressure (hPa)',fontsize=label_fontsize)
ax2.set_ylabel('Change in sea level pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

save_name = imagedir + "monthly_mean_slp_arctic_1948_2015.png"
plt.savefig(save_name, bbox_inches='tight')

#plt.show()


###########################
# Plot 2
###########################
y0 = np.zeros_like(slp_diff_maxmin).astype('f')


t = years

cint = 3
yticks2 = np.arange(-12,12+(cint/2),cint)
numticks = len(yticks2)
yticks3 = np.linspace(1004,1024,num=numticks)

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,slp_diff_maxmin,'g',linewidth=3.0,label='Difference in SLP')
p2, = ax1.plot(t,y0,'k',linewidth=3.0)
p1, = ax1.plot(t,ym2,'g--',linewidth=3.0,label='Long-term mean difference in SLP')


ax2 = ax1.twinx()
p2, = ax2.plot(t,apr,'b',linewidth=3.0,label='April SLP')
p3, = ax2.plot(t,july,'r',linewidth=3.0,label='July SLP')
p0, = ax2.plot(t,y0,'k',linewidth=3.0)
ax2.grid(True, linestyle='--')
ax2.set_ylim([yticks3[0],yticks3[-1]])

ax1.set_xlim([np.min(years),np.max(years)])
ax1.set_ylim([yticks2[0],yticks2[-1]])
ax1.yaxis.set_ticks(yticks2)
ax1.set_yticklabels(yticks2)
ax1.grid(True, linestyle='-')
ax1.set_ylabel('Change in sea level pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Year',fontsize=label_fontsize)

ax2.yaxis.set_ticks(yticks3)
ax2.set_yticklabels(yticks3)


legend = ax1.legend(loc='upper left', shadow=True)
legend = ax2.legend(loc='upper right', shadow=True)

save_name = imagedir + "april_vs_july_slp_arctic_1948_2015.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()
