# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

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
# September 2014
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

# readsounding is in weather_modules.py
import weather_modules as wm
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
from mstats import *
###########################
# User settings
###########################
datadir = '/home/scavallo/classes/metr_5004/fall_2013/classcode/'
imagedir = '/home/scavallo/classes/metr_5004/fall_2014/homework/'
infile1 = 'sounding_data_oun_2014102212.dat'
label_fontsize = 16
plot_logscale = 'false'
Lv = 2.5*10**6 # J kg-1
cp = 1004 # J K-1 kg-1

###########################
# End user settings
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
winddir=data["drct"]
windspeed_kt=data["sknt"]
Spd = windspeed_kt*0.514444


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
dws_dz = (ws[1::]-ws[0:-1]) / (hght[1::]-hght[0:-1])
des_dz = (es[1::]-es[0:-1]) / (hght[1::]-hght[0:-1])

temp_k = temp + 273.15 # Convert to Kelvin
theta_k = theta + 273.15
temp_avg = ( temp_k[1::] + temp_k[0:-1] ) / 2;
delta_z = (hght[1::]-hght[0:-1]) # meters
delta_theta_k = (theta_k[1::]-theta_k[0:-1]) 

u = -Spd * np.sin(winddir * (np.pi/180.))
v = -Spd * np.cos(winddir * (np.pi/180.))
delta_u = (u[1::]-u[0:-1]) 
delta_v = (v[1::]-v[0:-1]) 
u_avg = ( u[1::] + u[0:-1] ) / 2;
v_avg = ( v[1::] + v[0:-1] ) / 2;

Ri = ( ((9.81/temp_avg)*delta_theta_k)*delta_z ) / (delta_u**2 + delta_v**2);

du_dz = (u[1::]-u[0:-1]) / (hght[1::]-hght[0:-1])*10**3 # m s-1 km-1
dv_dz = (v[1::]-v[0:-1]) / (hght[1::]-hght[0:-1])*10**3 # m s-1 km-1

# Just for fun I am calculating the dry and moist adiabatic lapse rates for comparison.  You didn't need to do this.
sat_lapse = (9.81/cp) + (Lv/cp)*dws_dz
dry_lapse = np.zeros_like(sat_lapse).astype('f')
dry_lapse[:] = 9.81/cp
sat_lapse = sat_lapse*10**3
dry_lapse = dry_lapse*10**3

# Relative humidity
rh = ((mixr/1000)/ws)*100
rhi =((mixr/1000)/wsi)*100

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

# These are the levels I want to show on my plot explicitly       
yticks = [1000.,850.,700.,650.,500.,300.,100.]

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
p1, = ax1.plot(theta,pres,'r',linewidth=3.0,label='Potential temperature')
#p2, = ax1.plot(thetae,pres,'b',linewidth=3.0,label='Equivalent potential temperature')
if plot_logscale == 'true':
    ax1.set_yscale('log')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_xlim([260.,360.])
ax1.set_ylim([1000.,100.])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper left', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Potential temperature (K)',fontsize=label_fontsize)
save_name = imagedir + "theta_sounding.png"
#plt.savefig(save_name, bbox_inches='tight')

x0 = np.zeros_like(dtheta_dz).astype('f')
###########################
# Plot 2
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(dtheta_dz,pres[1:],'r',linewidth=2.5,label=r'$\partial{\theta}/\partial{z}$')
#p2, = ax1.plot(dthetae_dz,pres[1:],'b',linewidth=2.5,label=r'$\partial{\theta_e}/\partial{z}$')
p0, = ax1.plot(x0,pres[0:-1],'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
#ax1.set_xlim([-20,20])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Lapse rate (K km-1)',fontsize=label_fontsize)
save_name = imagedir + "dtheta_dthetae_sounding.png"
#plt.savefig(save_name, bbox_inches='tight')


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
p1, = ax1.plot(du_dz,pres[1:],'b',linewidth=2.5,label=r'$\partial{u}/\partial{z}$')
p2, = ax1.plot(dv_dz,pres[1:],'r',linewidth=2.5,label=r'$\partial{v}/\partial{z}$')
p0, = ax1.plot(x0,pres[0:-1],'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
#ax1.set_xlim([-20,20])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Lapse rate (K km-1)',fontsize=label_fontsize)
save_name = imagedir + "dudv_sounding.png"
#plt.savefig(save_name, bbox_inches='tight')


Ri_c = np.zeros_like(dtheta_dz).astype('f')
Ri_c[:] = 0.25
###########################
# Plot 4
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(Ri,pres[1:],'r',linewidth=2.5,label='Richardson Number')
p0, = ax1.plot(Ri_c,pres[0:-1],'k--',linewidth=3.5,label='Critical Richardson Number')
#p0, = ax1.plot(x0,pres[0:-1],'k',linewidth=3.5)
if plot_logscale == 'true':
    ax1.set_yscale('log')
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([1000.,100.])
ax1.set_xlim([-10,10])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper left', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Richardson Number',fontsize=label_fontsize)
save_name = imagedir + "richardson_sounding.png"
#plt.savefig(save_name, bbox_inches='tight')


plt.show()

# <codecell>


