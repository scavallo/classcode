# microphysics_homework.py
# 
#
# Steven Cavallo
# University of Oklahoma
# November 2018
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

import weather_modules as wm
from mstats import *
###########################
# User settings
###########################
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
label_fontsize = 16
plot_logscale = 'true'

###########################
# End user settings
###########################
# set a range of radii in microns
r = np.arange(0,12.1,0.01)
# Temperature and density
T = 283.15 # K
rho_L = 10**(3) # kg m-3
rho_I = 0.9167*10**(3.0)
sigma = 0.076164 # J m-2
ival = 2

#msol = 10**(-22)
Msol = (22.9898+35.453)*10**(-3) # NaCl
msol = 10**(-19)
#Msol = (132.14)*10**(-3) # ammonium sulfate 


rmeters = r*10**(-6)

#sigma = (rho_L*wm.Rv*T*a)/2
a = (2*sigma)/(rho_L*wm.Rv*T)

S_curvature = np.exp(a/rmeters)
rh_curvature = S_curvature*100

b = (4.3*ival*10**(-6)*msol)/Msol
S_solute = 1 - (b/(rmeters**3))
rh_solute = S_solute*100

S_total = S_solute*S_curvature
rh_total = S_total*100

# coefficients for GiSi:
aa = -2.4579*10**(-13.0)
bb = 1.8*10**(-10.0)
cc = -4.38*10**(-8.0)
dd = 3.5435*10**(-6.0)
TT = np.arange(-40,0.5,0.5)
[ind1] = np.argwhere(TT==-5)
[ind2] = np.argwhere(TT==-15)
TC = TT
TK = TC + 273.15

GiSi = (aa*TK**3) + (bb*TK**2) + (cc*TK) + dd
print(GiSi[ind1])
print(GiSi[ind2])

radius1 = 1*10**(-3.0)
h = 10*10**(-6.0)
radius2 = radius1 + (4*GiSi*3600)/(np.pi*rho_I*h)
radius1_mm = radius1*10**(3.0)
radius2_mm = radius2*10**(3.0)
delta_radius = radius2_mm - radius1_mm
radius_difference = radius2_mm[ind2] - radius2_mm[ind1]
mass_temp1 = rho_I*h*np.pi*((radius2[ind1])**(2.0))
mass_temp2 = rho_I*h*np.pi*((radius2[ind2])**(2.0))
mass_difference = mass_temp2 - mass_temp1

print("The radius of droplet growing at temperature %.2f deg C is %.2f mm" % (TC[ind1], radius2_mm[ind1]))
print("The radius of droplet growing at temperature %.2f deg C is %.2f mm" % (TC[ind2], radius2_mm[ind2]))
print("The increase in radius of the droplet growing at temperature %.2f deg C is %.2f mm" % (TC[ind1], delta_radius[ind1]))
print("The increase in radius of the droplet growing at temperature %.2f deg C is %.2f mm" % (TC[ind2], delta_radius[ind2]))
print("The radius of the droplet growing at temperature %.2f deg C is %.2f mm larger than the droplet growing at temperature %.2f deg C" % (TC[ind2], radius_difference, TC[ind1]))
print("The mass of droplet growing at temperature %.2f deg C is %.4g kg" % (TC[ind1], mass_temp1))
print("The mass of droplet growing at temperature %.2f deg C is %.4g kg" % (TC[ind2], mass_temp2))
print("The difference in mass between the droplet growing at temperature %.2f deg C and the droplet growing at %.2f deg C is %.4g kg" % (TC[ind1], TC[ind2], mass_difference))



#xticks = [0.05,0.1,1,5,10]
xticks = [0.1,0.5,1,5,10]
###########################
# Plot 1
###########################
x0 = np.zeros_like(rh_curvature).astype('f')
x0[:] = 100.
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(r,rh_curvature,'b',linewidth=2.5,label='Curvature')
p2, = ax1.plot(r,rh_solute,'r',linewidth=2.5,label='Solution')
p3, = ax1.plot(r,rh_total,'g',linewidth=2.5,label='Both')
p4, = ax1.plot(r,x0,'k',linewidth=4.0)
if plot_logscale == 'true':
    ax1.set_xscale('log')
    
ax1.set_ylim([99,101])
ax1.set_xlim([0.05,10.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Relative humidity (%) ',fontsize=label_fontsize)
ax1.set_xlabel('Radius (micrometers)',fontsize=label_fontsize)
save_name = imagedir + "hw3_curvature.png"
plt.savefig(save_name, bbox_inches='tight')


###########################
# Plot 2
###########################
x0 = np.zeros_like(rh_curvature).astype('f')
x0[:] = 100.
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(TC,delta_radius,'b',linewidth=2.5,label='Curvature')
ax1.set_xlim([-40,0.])
plt.gca().invert_xaxis()
ax1.grid(True, linestyle='-')
#legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Increase in radius size (mm) ',fontsize=label_fontsize)
ax1.set_xlabel(r'Temperature ($^{\circ}$C)',fontsize=label_fontsize)
save_name = imagedir + "hw3_deposition.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()


