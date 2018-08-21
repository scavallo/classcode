###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

import utilities_modules as um
from mstats import *


###########################
# User settings
###########################
datadir = ''
imagedir = ''
infile = 'slp_6090_annual_1948_2017.dat'

plot_option = 2 # 1 for SLP, 2 for mass-averaged velocity
calculate_running_mean = 'True'
label_fontsize = 16
title_fontsize = 18
latitude = 60.

conversion_factor = 1000. # for converting meters to mm
climo_years = [1981,2010]
analysis_years1 = [2016,2016]
analysis_years2 = [2017,2017]

ylims1 = [1005, 1020]
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

ab = np.loadtxt(datadir+infile, skiprows=1)       
years = ab[:,0]
slp_jan = ab[:,0]
slp_july = ab[:,6]
slp_full = ab[:,1:]

slp_full = um.filter_numeric_nans(slp_full,0,float('NaN'),'low')

[ind1a] = np.where(years == analysis_years1[0])
[ind1b] = np.where(years == analysis_years1[1])
[ind2a] = np.where(years == analysis_years2[0])
[ind2b] = np.where(years == analysis_years2[1])

[climind1a] = np.where(years == climo_years[0])
[climind1b] = np.where(years == climo_years[1])

slp_monthly_mean = np.nanmean(slp_full,0)

slp_subset1 = np.nanmean(slp_full[ind1a:ind1b+1,:],0)
slp_subset2 = np.nanmean(slp_full[ind2a:ind2b+1,:],0)
slp_climatology = np.nanmean(slp_full[climind1a:climind1b+1,:],0)

slp_anomaly1 = slp_subset1 - slp_climatology
slp_anomaly2 = slp_subset2 - slp_climatology


# mass transport calculations
dslp_monthly_mean = np.zeros(len(slp_monthly_mean))
dslp_climatology = np.zeros(len(slp_climatology))
dslp_subset1 = np.zeros(len(slp_subset1))
dslp_subset2 = np.zeros(len(slp_subset2))

dslp_monthly_mean[-1] = slp_monthly_mean[0] - slp_monthly_mean[-1]
dslp_monthly_mean[0:-1] = slp_monthly_mean[1:] - slp_monthly_mean[0:-1]

dslp_climatology[-1] = slp_climatology[0] - slp_climatology[-1]
dslp_climatology[0:-1] = slp_climatology[1:] - slp_climatology[0:-1]

dslp_subset1[-1] = slp_subset1[0] - slp_subset1[-1]
dslp_subset1[0:-1] = slp_subset1[1:] - slp_subset1[0:-1]

dslp_subset2[-1] = slp_subset2[0] - slp_subset2[-1]
dslp_subset2[0:-1] = slp_subset2[1:] - slp_subset2[0:-1]


dpdt_monthly_mean = (dslp_monthly_mean*100.)/(np.array(ndays)*86400.) # Pa s-1    
dpdt_climatology = (dslp_climatology*100.)/(np.array(ndays)*86400.) # Pa s-1  
dpdt_subset1 = (dslp_subset1*100.)/(np.array(ndays)*86400.) # Pa s-1  
dpdt_subset2 = (dslp_subset2*100.)/(np.array(ndays)*86400.) # Pa s-1  

h = R_e*(1-sinlatr)
vm_monthly_mean = h/(p0*coslatr)*(dpdt_monthly_mean)*conversion_factor
vm_climatology = h/(p0*coslatr)*(dpdt_climatology)*conversion_factor
vm_subset1 = h/(p0*coslatr)*(dpdt_subset1)*conversion_factor
vm_subset2 = h/(p0*coslatr)*(dpdt_subset2)*conversion_factor


vm_subset1_anomaly = vm_subset1 - vm_climatology
vm_subset2_anomaly = vm_subset2 - vm_climatology



t = np.arange(1,13,1)
z0 = np.zeros(len(t))

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 


###########################
# Plots
###########################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,slp_monthly_mean,'b',linewidth=3.0,label='Monthly mean Arctic SLP (1948-2017)')
p2, = ax1.plot(t,slp_climatology,'b--',linewidth=3.0,label='Monthly mean Arctic SLP climatology')
p3, = ax1.plot(t,slp_subset1,'r',linewidth=2.0,label='2016')
p4, = ax1.plot(t,slp_subset2,'k',linewidth=2.0,label='2017')
ax1.set_xlim([np.min(t),np.max(t)])
ax1.set_ylim([ylims1[0],ylims1[-1]])
ax1.xaxis.set_ticks(t[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)

ax1.set_title('Monthly-mean sea level pressure over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('hPa',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

plt.savefig(imagedir + 'slp_60N_comparisions' + '.png', bbox_inches='tight')



fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,dslp_monthly_mean,'b',linewidth=3.0,label=r'$\Delta$SLP')
p2, = ax1.plot(t,dslp_subset1,'r',linewidth=3.0,label=r'$\Delta$SLP (2016)')
p3, = ax1.plot(t,dslp_subset2,'g',linewidth=3.0,label=r'$\Delta$SLP (2017)')
ax1.set_xlim([np.min(t),np.max(t)])
ax1.xaxis.set_ticks(t[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
p0, = ax1.plot(t,np.zeros(len(t)),'k',linewidth=1.0)

ax1.set_title('Monthly-mean sea level pressure tendency over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('hPa',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

plt.savefig(imagedir + 'deltaslp_60N_comparisions' + '.png', bbox_inches='tight')



fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,vm_monthly_mean,'b',linewidth=3.0,label='Monthly mean (1948-2017)')
p2, = ax1.plot(t,vm_climatology,'b--',linewidth=3.0,label='Climatological (1981-2010)')
p3, = ax1.plot(t,vm_subset1,'r',linewidth=2.0,label='2016')
p4, = ax1.plot(t,vm_subset2,'g',linewidth=2.0,label='2017')
ax1.set_xlim([np.min(t),np.max(t)])
ax1.xaxis.set_ticks(t[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower right', shadow=True)
p0, = ax1.plot(t,np.zeros(len(t)),'k',linewidth=2.0)

ax1.set_title('Mass-weighted velocity transport over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('mm s-1',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)
plt.savefig(imagedir + 'vm_60N_comparisions' + '.png', bbox_inches='tight')


fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,vm_subset1_anomaly,'r',linewidth=2.0,label='2016 anomaly')
p2, = ax1.plot(t,vm_subset2_anomaly,'g',linewidth=2.0,label='2017 anomaly')
ax1.set_xlim([np.min(t),np.max(t)])
ax1.xaxis.set_ticks(t[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower right', shadow=True)
p0, = ax1.plot(t,np.zeros(len(t)),'k',linewidth=2.0)

ax1.set_title('Mass-weighted velocity transport anomaly over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('mm s-1',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)
plt.savefig(imagedir + 'vm_anomaly_60N_comparisions' + '.png', bbox_inches='tight')
plt.show()

