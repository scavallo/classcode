# microphysics_homework.py
# 
#
# Steven Cavallo
# University of Oklahoma
# October 2020
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
r = np.arange(0,12.1,0.001)
# Temperature and density
T = 283.15 # K
rho_L = 10**(3) # kg m-3
rho_I = 0.9167*10**(3.0)
sigma = 0.076164 # J m-2
ival = 2.0

msol = 10**(-19.) 
Msol_NaCl = (22.9898+35.453)*10**(-3.) # NaCl
# Enter this
Msol = # ammonium sulfate 

rmeters = r*10**(-6)

a = (2*sigma)/(rho_L*wm.Rv*T)

# Curvature effect
S_curvature = # Enter this
rh_curvature = S_curvature*100

# Solute effect
b = (4.3*ival*10**(-6)*msol)/Msol
b_NaCl = (4.3*ival*10**(-6)*msol)/Msol_NaCl
b_NaCl = (4.3*ival*msol*10**(-6.))/(Msol_NaCl)

S_solute = # Enter this
rh_solute = S_solute*100.

# Combination of curvature and solute effect
S_total = # Enter this
rh_total = S_total*100.

S_solute_NaCl = 1.0 - (b_NaCl/(rmeters**3.0))
rh_solute_NaCl = S_solute_NaCl*100.
S_total_NaCl = S_solute_NaCl*S_curvature
rh_total_NaCl = S_total_NaCl*100.

xticks = [0.01,0.1,0.5,1,5,10]
###########################
# Figure
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
p1, = ax1.plot(r,rh_curvature,'b',linewidth=4,label='Curvature')
p2a, = ax1.plot(r,rh_solute_NaCl,'r--',linewidth=4,label='Solute Effect (Sodium Chloride)')
p2b, = ax1.plot(r,rh_solute,'r',linewidth=4,label='Solute Effect (Ammonium Sulfate)')
# Add plot overlays here

p5, = ax1.plot(r,x0,'k',linewidth=4.0)
if plot_logscale == 'true':
    ax1.set_xscale('log')
    
ax1.set_ylim([99,101])
ax1.set_xlim([0.01,10.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Relative humidity (%) ',fontsize=label_fontsize)
ax1.set_xlabel('Radius (micrometers)',fontsize=label_fontsize)
save_name = imagedir + "hw5_kohler.png"
plt.savefig(save_name, bbox_inches='tight')


plt.show()


