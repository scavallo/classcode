
# coding: utf-8

# In[1]:

# moist_ascent_latent_heating.py
# 
#
# Computes moist adiabatic parcel ascent and 
# latent heating rates
# 
# Steven Cavallo
# University of Oklahoma
# September 2017
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

import weather_modules as wm
#import imp
#wm = imp.load_source('weather_modules','/Users/scavallo/scripts/python_scripts/weather_modules.py')
from mstats import *


###########################
# User options
###########################
pres1 = 70000. # Pa
pres2 = 30000. # Pa
omega = -1.0 # Pa s-1
reference_latitude = 45.

read_temperature_fromfile = 'True'
label_fontsize = 16
datadir = ''
imagedir = ''
infile1 = 'nnrp_t700_vs_lat_monthly.dat'


# In[3]:

###########################
# constants
###########################
cp = 1004.  # J K-1 kg-1
cv = 717.    # J K-1 kg-1
Rd = 287.    # J K-1 kg-1
Rv = 461.    # J K-1 kg-1
g  = 9.81    # m s-2
p0 = 100000. # Pa
epsi = Rd/Rv
Lv = 2.5*10**(6.0) # J kg-1
p0 = 100000. # Pa


# In[4]:

gammad = g/cp # Dry adiabatic lapse rate is a constant

if read_temperature_fromfile == 'True':
    ab = np.loadtxt(datadir+infile1, skiprows=1)       
    lats = ab[:,0]
    T1a = ab[:,1]
    T1b = ab[:,7]
    [refind] = np.where(lats==reference_latitude)
else:
    T1 = np.arange(-40,41,0.25)+273.15 # Virtual temperatures
    T1a = T1
    [refind] = np.where(T1a == 273.15)
es1 = wm.claus_clap(T1a) # saturation vapor pressure
ws1a = (epsi*es1)/(pres1-es1) # saturation mixing ratio
es1 = wm.claus_clap(T1b) # saturation vapor pressure
ws1b = (epsi*es1)/(pres1-es1) # saturation mixing ratio

# Compute moist adiabatic lapse rate, first in z-coords
fa = (( (1 + (Lv*ws1a)/(Rd*T1a))) / (1.0 + (ws1a*(Lv**2)/(cp*Rv*T1a**2.0)) ) )
fb = (( (1 + (Lv*ws1b)/(Rd*T1b))) / (1.0 + (ws1b*(Lv**2)/(cp*Rv*T1b**2.0)) ) )
sat_lapse_a = gammad*fa
sat_lapse_b = gammad*fb

Rho_a=pres1/(Rv*T1a) # Density, moist air
Rho_b=pres1/(Rv*T1b) 
wa = -omega/(Rho_a*g)
print('w in profile A is %8.5f m s-1' %(wa[refind]))
wb = -omega/(Rho_b*g)
print('w in profile B is %8.5f m s-1' %(wb[refind]))

deltap = pres2-pres1
deltaz_a = (-1.0*deltap)/(Rho_a*g)
deltaz_b = (-1.0*deltap)/(Rho_b*g)
print('Cloud depth is %7.2f m in profile A and %7.2f in profile B' % (deltaz_a[refind],deltaz_b[refind]))

#sat_lapse_isobaric = -sat_lapse/(Rho*9.81) # Convert to p-coords

# Temperature at top of ascent layer (T2) is equal to what it was initially (T1) plus slope (-gamma_m) * deltaz
#T2 = T1 - (sat_lapse_isobaric)*(pres2-pres1) 
T2a = T1a + (-1.0*sat_lapse_a*deltaz_a)
T2b = T1b + (-1.0*sat_lapse_b*deltaz_b)
print('Lower- and upper-level temperatures in profile A are %7.2f K and %7.2f' %(T1a[refind],T2a[refind]))
print('Lower- and upper-level temperatures in profile B are %7.2f K and %7.2f' %(T1b[refind],T2b[refind]))

# Convert to Celsius (needed for plot if not reading from file)
#T1C = T1 - 273.15
#T2C = T2 - 273.15

# Convert to K/km 
sat_lapse_a = sat_lapse_a*1000. # per km
sat_lapse_b = sat_lapse_b*1000. # per km

# Fractional departure from dry adiabatic
fd_a = 1.0 - fa
fd_b = 1.0 - fb


# In[5]:

exner = (((pres1+pres2)/2.0)/p0)**(Rd/cp)

###################
# Profile A
###################
es2 = wm.claus_clap(T2a)
ws2 = (epsi*es2)/(pres2-es2)

qs1 = ws1a/(1+ws1a)
qs2 = ws2/(1+ws2)

deltaqs = (qs2 - qs1)
deltap = pres2 - pres1
deltat = deltap/omega # Time it takes parcel to travel deltap (through cloud)

dthetadt_a = -1.0*(Lv/(exner*cp))*(deltaqs/deltat)
dthetadt_perday_a = dthetadt_a * 86400.

###################
# Profile B
###################
es2 = wm.claus_clap(T2b)
ws2 = (epsi*es2)/(pres2-es2)

qs1 = ws1b/(1+ws1b)
qs2 = ws2/(1+ws2)

deltaqs = (qs2 - qs1)
deltap = pres2 - pres1
deltat = deltap/omega # Time it takes parcel to travel deltap (through cloud)

dthetadt_b = -1.0*(Lv/(exner*cp))*(deltaqs/deltat)
dthetadt_perday_b = dthetadt_b * 86400.


# In[6]:

y0 = np.zeros_like(T1a).astype('f') #
x = np.linspace(0,len(T1a),np.size(T1a))
if read_temperature_fromfile == 'True':
    xplot = lats
    xlabel = r'Latitude ($^{\circ}$N)'
    xlims = [np.min(lats), np.max(lats)]
else:
    xplot = T1C
    xlims = [-40.,40.]
    xlabel = r'Temperature ($^{\circ}$C)'

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2)


# In[7]:

nticks = 9
yticks_ax1 = np.linspace(2,10,nticks)
yticks_ax2 = np.linspace(0,0.8,nticks)  

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
ax2 = ax1.twinx()
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1a, = ax1.plot(xplot,sat_lapse_a,'r',linewidth=3.0,label=r'$\Gamma_m$ (January)')
p1a, = ax1.plot(xplot,sat_lapse_b,'r--',linewidth=3.0,label=r'$\Gamma_m$ (July)')
p2a, = ax2.plot(xplot,fd_a,'g',linewidth=3.0,label='1-f (January)')
p2b, = ax2.plot(xplot,fd_b,'g--',linewidth=3.0,label='1-f (July)')

legend = ax1.legend(loc='lower left', shadow=True)
legend = ax2.legend(loc='upper left', shadow=True)

ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])

ax1.set_ylim([yticks_ax1[0],yticks_ax1[1]])
ax2.set_ylim([yticks_ax2[0],yticks_ax2[1]])
ax1.yaxis.set_ticks(yticks_ax1)
ax1.set_yticklabels(yticks_ax1)
ax2.yaxis.set_ticks(yticks_ax2)
ax2.set_yticklabels(yticks_ax2)

ax1.set_ylabel('Moist adiabatic lapse rate',fontsize=label_fontsize)
ax2.set_ylabel('Fractional departure from dry adiabatic',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + "moist_adiabatic_hw3_2017.png"
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


# In[8]:

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1a, = ax1.plot(xplot,fd_a,'r',linewidth=3.0,label='1-f (January)')
p1b, = ax1.plot(xplot,fd_b,'r--',linewidth=3.0,label='1-f (July)' )
legend = ax1.legend(loc='lower left', shadow=True)
ax1.axvline(xplot[refind],linewidth=4, color='k')
ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])
ax1.set_ylabel('Fractional departure from dry adiabatic',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + "moist_adiabatic_factor_hw3_2017.png"
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


# In[9]:

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(xplot,dthetadt_perday_a,'b',linewidth=3.0,label=r'$\theta$ latent heating rate (January)')
p2, = ax1.plot(xplot,dthetadt_perday_b,'b--',linewidth=3.0,label=r'$\theta$ latent heating rate (July)')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])
ax1.set_ylabel('Latent heating rate (K day-1)',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + "latent_heating_rates_hw3_2017.png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()

