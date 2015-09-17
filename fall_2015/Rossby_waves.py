
# coding: utf-8

# In[ ]:

# Rossby_waves.py
#
# Script to plot Rossby wave phase speed and group velocity, given a background flow, wavenumber, and latitude
#
# Steven Cavallo
# University of Oklahoma
# September 2015


# In[1]:

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
from mstats import *


# In[2]:

###########################
# User parameters
###########################
ubar = 22.0
num_waves_merid = 0.0
num_waves_zonal = 5.0
delta_latitude = 0.5;
lat_info = 35.0

label_fontsize = 16
imagedir = '/Users/scavallo/scripts/python_scripts/images/'


# In[3]:

x = np.linspace(0,2*np.pi,100)
omeg_e = (2*np.pi) / (24*3600)
latitudes = np.arange(-90,90+(delta_latitude/2),delta_latitude)

dlat = np.zeros_like(latitudes).astype('f')
dlat[1:-1] = (latitudes[2:] - latitudes[:-2])/2
dlat[0] = (latitudes[1]-latitudes[0]) 
dlat[-1] = (latitudes[-1]-latitudes[-2])

f = 2.0*omeg_e*np.sin(latitudes*(np.pi/180))
dy = 111.0*dlat*1000

beta = np.zeros_like(latitudes).astype('f')
beta[1:-1] = (f[2:] - f[:-2])/(2*dy[1:-1])
beta[0] = (f[1] - f[0])/(dy[1]) 
beta[-1] = (f[-1] - f[-2])/(dy[-1]) 


# In[4]:

pole_to_pole = 180.0*111000.0
circums_earth = 360.*111000.*np.cos(latitudes*(np.pi/180.))
ind = np.where(latitudes==lat_info)
ind_eq = np.where(latitudes==0.0)
print "Beta at latitude %2.1f degrees is %e" % (latitudes[ind], beta[ind])
print "The distance around the Earth along latitude %2.1f degrees is %2.2f km" %(latitudes[ind], circums_earth[ind]/1000)
print "The distance around the Earth at the equator is %2.2f km" %(circums_earth[ind_eq]/1000)


# In[5]:

k = num_waves_zonal*((2.*np.pi)/circums_earth)
#print k[ind]
l = num_waves_merid*((2.*np.pi)/pole_to_pole)
disp = ubar*k - ((k*beta)/(k**2.0 + l**2.0))
c = disp/k
cg = ubar - (beta*(l**2.0 - k**2.0))/((k**2.0 + l**2.0)**2.0) 

cg_over_c = np.abs(cg/c)

ubar_crit = beta/(k**2+l**2)
print "The critical background wind at latitude %2.1f degrees and wave number %d is %2.2f m s-1" %(latitudes[ind], num_waves_zonal, ubar_crit[ind])


# In[6]:

ubar_jja = 12.5
ubar_djf = 22

k_jja = np.sqrt(beta[ind]/ubar_jja)
k_djf = np.sqrt(beta[ind]/ubar_djf)
num_waves_zonal_jja = k_jja/((2.*np.pi)/circums_earth[ind])
num_waves_zonal_djf = k_djf/((2.*np.pi)/circums_earth[ind])
#print num_waves_zonal_jja, num_waves_zonal_djf
print "There are %2.1f zonal waves in JJA and %2.1f zonal waves in DJF" %(num_waves_zonal_jja, num_waves_zonal_djf)


# In[7]:

#mstats(c)
#mstats(cg)


# In[8]:

xticks = [-90,-60,-30,0,30,60,90]
yticks = np.arange(-60,60.5,10)

latline = np.zeros(np.size(beta))
latline[:] = lat_info

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,c,'r',linewidth=3.0,label=r'c')
p2, = ax1.plot(latitudes,cg,'b',linewidth=3.0,label=r'$c_g$')



ax1.grid(linestyle='-')

ax1.set_xlim([-90.,90.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

ax1.set_ylim(yticks[0],yticks[-1])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)

ax1.set_ylabel(r'Phase speed (m s$^{-1}$)',fontsize=label_fontsize)
ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + "rossby_waves.png"
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


# In[9]:

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,cg_over_c,'r',linewidth=3.0,label=r'$\frac{c_g}{c}$')
ax1.grid(linestyle='-')

ax1.set_xlim([-90.,90.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

ax1.set_ylim([1.,5.])

ax1.set_ylabel(r'$\frac{c_g}{c}$',fontsize=label_fontsize)
ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)



save_name = imagedir + "rossby_waves_cg_over_c.png"
#plt.savefig(save_name, bbox_inches='tight')
plt.show()

