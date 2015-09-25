
# coding: utf-8

# In[1]:

# plot_sounding_data.py
# 
# Example plotting script to read in sounding data and plot.
#
# Obtain sounding data from the University of Wyoming website:
# http://weather.uwyo.edu/upperair/sounding.html
# The function readsounding (located within weather_modules.py) will read data 
#    in the format given when type of plot is 'Text: List' from the dropdown menu
# 
# Steven Cavallo
# University of Oklahoma
# September 2015
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

# readsounding is in weather_modules.py
#import weather_modules as wm
import imp
wm = imp.load_source('weather_modules','/Users/scavallo/scripts/python_scripts/weather_modules.py')
#sm = imp.load_source('sounding_modules','/Users/scavallo/scripts/python_scripts/sounding_modules.py')
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
from mstats import *


# In[2]:

###########################
# User settings
###########################
datadir = '/Users/scavallo/data/'
imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2015/homework/images/'
#infile1 = 'sounding_data_dtx_2014010700.dat'
#infile1 = 'sounding_data_mpx_2014010400.dat'
infile1 = 'sounding_data_hw3_fall2015.dat'
label_fontsize = 16
plot_logscale = 'false'
Lv = 2.5*10**6 # J kg-1
cp = 1004 # J K-1 kg-1
Rd = 287. # J K-1 kg-1
Rv = 461. # J K-1 kg-1
epsi = Rd/Rv


# In[3]:

###########################
# Read in data
###########################
data,fieldnames = wm.readsounding(datadir+infile1)
# Print the plotting options
print fieldnames
 
 # Read in and separate out the fields into their own arrays  
pres=data["pres"]
temp=data["temp"]
dwpt=data["dwpt"]
hght=data["hght"]
theta=data["thta"]
thetae=data["thte"]
mixr=data["mixr"]
thetav=data["thtv"]


# In[4]:

# To compute relative humidity, will need saturation mixing ratio.  
# Saturation mixing ratio depends on saturation vapor pressure and total pressure.
es = wm.claus_clap(temp+273.15)
esi = wm.claus_clap_ice(temp+273.15)
ws = wm.satur_mix_ratio(es, pres*100)
wsi = wm.satur_mix_ratio(esi,pres*100)

# Calculate vertical derivatives.  If you didn't know how to do this, there are numerous resources, including examples within weather_modules.py.  
# A quick search in Google will also give you some ways to do it. Be resourceful!
# Below, I am using forward differencing in z
dtemp_dz = ((temp[1::]-temp[0:-1]) / (hght[1::]-hght[0:-1])*10**3)*-1 # K km-1
dtheta_dz = (theta[1::]-theta[0:-1]) / (hght[1::]-hght[0:-1])*10**3 # K km-1
dthetae_dz = (thetae[1::]-thetae[0:-1]) / (hght[1::]-hght[0:-1])*10**3 # K km-1

# Dry and moist adiabatic lapse rates 
sat_lapse = np.zeros_like(temp).astype('f')
sat_lapse = (9.81/cp)*( (1 + (Lv*ws)/(Rd*(temp+273.15)))) / (1 + (ws*(Lv**2)/(cp*Rv*(temp+273.15)**2)) )
dry_lapse = np.zeros_like(sat_lapse).astype('f')
dry_lapse[:] = 9.81/cp
sat_lapse = sat_lapse*10**3.0
dry_lapse = dry_lapse*10**3.0

# Relative humidity
e = ((mixr/1000.)/((mixr/1000.) + epsi))*(pres*100)
rh = (e/es)*100
rhi =(e/esi)*100

rh_eff = np.zeros_like(rh).astype('f')
for ii in range(0, len(rh)-1):
   if (temp[ii]<-0):
       rh_eff[ii] = rhi[ii]
   else:
       rh_eff[ii] = rh[ii]   

# Note that I did NOT ask you to find the dendritic growth zone.  
# I am doing this for illustration.
rh_mask = np.zeros_like(rh).astype('f')
rh_mask2 = np.zeros_like(rh).astype('f')
inds = np.argwhere(rh_mask==0)   
rh_mask[inds] = float('NaN')
inds = np.argwhere(rh_mask2==0)   
rh_mask2[inds] = float('NaN')
inds2 = np.where( (temp<-12) & (temp>-18))
for ii in range(0, len(rh)-1):
   if (temp[ii]<-12 and temp[ii]>-18):
       rh_mask[ii] = rh[ii]
       rh_mask2[ii] = rhi[ii]


# In[5]:

# These are the levels I want to show on my plot explicitly       
#yticks = [1000.,850.,700.,600.,500.,300.,100.]

pcint = 100.
yticks = np.arange(1000.0,100.-(pcint/2),-pcint)
r100 = np.zeros_like(rh).astype('f')
r100[:] = 100.
###########################
# Plot 1
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(rh,pres,'g',linewidth=3.0,label='rh')
p2, = ax1.plot(rh_eff,pres,'b',linewidth=3.0,label=r'rh$_i$')
p3, = ax1.plot(rh_mask,pres,'r',linewidth=3.0,label='Dendritic growth zone')
p4, = ax1.plot(rh_mask2,pres,'r',linewidth=3.0)
p0, = ax1.plot(r100,pres,'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')

    ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)

ax1.set_xlim([20.,105.])


ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Relative humidity (Percent)',fontsize=label_fontsize)
save_name = imagedir + "hw3_rh.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[6]:

cint = 5
xticks2 = np.arange(-30,30+(cint/2),cint)
numticks = len(xticks2)
xticks1 = np.linspace(260,365,num=numticks)
xticks1 = np.round(xticks1)

x0 = np.zeros_like(dtheta_dz).astype('f')
###########################
# Plot 2
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
#for tick in ax1.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(label_fontsize)
#    tick.label1.set_fontweight('bold')
p1, = ax1.plot(theta,pres,'r',linewidth=3.0,label=r'$\theta$')
p2, = ax1.plot(thetae,pres,'b',linewidth=3.0,label=r'$\theta_e$')

ax2 = ax1.twiny()
p3, = ax2.plot(dtheta_dz,pres[1:],'m',linewidth=3.0,label=r'$\partial{\theta}/\partial{z}$')
p4, = ax2.plot(dthetae_dz,pres[1:],'g',linewidth=3.0,label=r'$\partial{\theta_e}/\partial{z}$')
p0, = ax2.plot(x0,pres[0:-1],'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')

ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)

ax2.set_ylim([1000.,100.])
ax2.yaxis.set_ticks(yticks)
ax2.set_yticklabels(yticks)

#ax1.set_xlim([260.,360.])
ax1.set_xlim([xticks1[0],xticks1[-1]])
ax1.xaxis.set_ticks(xticks1)
ax1.set_xticklabels(xticks1)
ax2.set_xlim([xticks2[0],xticks2[-1]])
ax2.xaxis.set_ticks(xticks2)
ax2.set_xticklabels(xticks2)

ax1.grid(True, linestyle='-')
ax2.grid(True, linestyle='--')



legend = ax1.legend(loc='upper left', shadow=True)
legend = ax2.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Potential temperature (K)',fontsize=label_fontsize)
ax2.set_xlabel('Lapse rate (K km-1)',fontsize=label_fontsize)
save_name = imagedir + "hw3_theta_thetae.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[7]:

###########################
# Plot 3
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(dtemp_dz,pres[1:],'g',linewidth=2.5,label=r'$\Gamma$')
p2, = ax1.plot(sat_lapse,pres,'b',linewidth=2.5,label=r'$\Gamma_s$')
p3, = ax1.plot(dry_lapse,pres,'r',linewidth=2.5,label=r'$\Gamma_d$')
p0, = ax1.plot(x0,pres[0:-1],'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
ax1.set_xlim([-20,20])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Lapse rate (K Km-1)',fontsize=label_fontsize) 
save_name = imagedir + "hw3_lapes_rates.png"
#plt.savefig(save_name, bbox_inches='tight')
plt.show()


# In[ ]:



