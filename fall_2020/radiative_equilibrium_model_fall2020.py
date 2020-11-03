#!/usr/bin/env python
#
# radiative_equilibrium_model
# 
# Simple model of any planet's radiative equilibrium state
#
# Steven Cavallo
# University of Oklahoma
# November 2020
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

###########################
# User settings and options
###########################
label_fontsize = 16
plot_greenhouse_effect = 'false'
albedo_planet = 0.3
emiss_obs = 0.77 # for plotting purposes only
astronimcal_units = 1. # Earth = 1; Ross 128-b = 1/20 
solar_luminosity_factor =1. # Earth = 1, Ross 128-b = 0.00362
Rd_planet = 287. # Earth = 287, Titan = 296.8, Ross 128-b = ?
g_planet = 1.12*9.81 # Earth = 9.81, Titan = 1.35, Ross 128-b = 1.12*g_earth 
p0_planet = 1.35*1000. # Earth = 1000, Titan = 1467, Ross 128-b = 1.35*g_earth 

imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
###########################
# Constants and setup
###########################
emiss = np.arange(0.0,1.01,0.01)
sigma = 5.67*10**(-8)
luminosity_sun = 3.867*10**26.0
earth_sun_dist = 1.5*10**11.0
planet_sun_dist = astronimcal_units*earth_sun_dist 

solar_luminosity = solar_luminosity_factor*luminosity_sun
###########################
# Calculations
###########################
Fs_planet = solar_luminosity/(4.0*np.pi*(planet_sun_dist)**2.0)
print(Fs_planet)

emiss_find = 0.77
eind = np.where(emiss==emiss_find)

# Equivalent blackbody temperature
Te = (((1.0-albedo_planet)*Fs_planet)/(4.0*sigma))**(1.0/4.0)
print("The equiv. blackbody temperature is %.2f K" % (Te))
Ts =   ((2.0/(2.0-emiss)))**(1.0/4.0)*Te
Ta = (  (emiss/(2.0-emiss)))**(1.0/4.0)*Te
OLR = sigma*Te**4

# Lapse rate
lapse = Ts - Ta
GHE = sigma*Ts**(4.0) - OLR

scale_height = (Rd_planet*((Ta+Ts)/2))/g_planet
phalf = p0_planet / 2.
z = -scale_height*np.log((p0_planet/2.0)/p0_planet)
print("The height corresponding to %.2f hPa is %.2f m" % (phalf, z[eind]))

lapse_rate = (lapse/z)*1000
print("The lapse rate is %.2f K/km" % (lapse_rate[eind]))

print("The Outgoing longwave radiation (OLR) is %.2f W m-2" % (OLR))
print("The surface and atmosphere temperatures corresponding to an atmospheric emissivity of %.2f are %.2f K and %.2f K" % (emiss[eind], Ts[eind], Ta[eind]))
print("Ta minus Ts and greenhouse effect corresponding to an atmospheric emissivity of %.2f are %.2f K and %.2f W m-2 " % (emiss[eind], lapse[eind], GHE[eind]))
print(lapse_rate[eind])


# Climate sensitivity calculations; not used in homework
clim_sens = 0.7
deltaF_sens = 10.
a_sens = 0.3
Fin_sens = 1368.

deltaT_sens = clim_sens*deltaF_sens
print(deltaT_sens)

epsi_sensi = 2.*(1.0-((1.0-a_sens)/(sigma*(288.+deltaT_sens)**(4.0)))*(Fin_sens/4.0))
print(epsi_sensi)


# plot setup
[inds] = np.where(Ts >= 288.)
print(emiss[inds[0]])
###########################
# Plots
###########################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(emiss,Ts,'r',linewidth=3.0,label=r'$T_s$')
p2, = ax1.plot(emiss,Ta,'b',linewidth=3.0,label=r'$T_a$')
plt.axvline(emiss_obs,color='k',linewidth=3.0,linestyle='dashed',label='Earth atmospheric emissivity')
ax1.set_ylabel('Temperature (Kelvin)',fontsize=label_fontsize)

plt.grid(True, linestyle='-')
ax2 = ax1.twinx()
if plot_greenhouse_effect == 'false':
    p4 = ax2.plot(emiss,lapse_rate,'r--',linewidth=3.0,label=r'$\Gamma$')
    ax2.set_ylabel(r'$\Gamma$ (K km$^{-1}$)',fontsize=label_fontsize)
    legend = ax1.legend(loc='center right', shadow=True)
    legend = ax2.legend(loc='lower left', shadow=True)    
else:
    p4 = ax2.plot(emiss,GHE,'g--',linewidth=3.0,label='GHE')
    ax2.set_ylabel(r'Greenhouse effect (W m$^{-2}$)',fontsize=label_fontsize)
    legend = ax1.legend(loc='lower right', shadow=True)
    legend = ax2.legend(loc='lower left', shadow=True)    

ax1.set_ylim([150.,330.])
ax1.set_xlim([0.,1.])
ax1.set_xlabel(r'Emissivity, $\varepsilon$',fontsize=label_fontsize)
save_name = imagedir + 'radiative_equilibrium_model_Ross128b.png'
plt.savefig(save_name, bbox_inches='tight')  
plt.show()




