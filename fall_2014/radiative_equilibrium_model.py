# radiative_equilbrium_model.py
# 
# Simple one-layer radiative equilbrium model
#
# Steven Cavallo
# University of Oklahoma
# November 2014
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime


###########################
# User settings
###########################
imagedir = ''
label_fontsize = 16
plot_greenhouse_effect = 'false'

a = 0.3
Fin = 1368
sigma = 5.67*10**(-8)
emiss = np.arange(0.0,1.01,0.01)
emiss_obs = 0.77
emiss_find = 0.77
###########################
# END user settings
###########################
eind = np.where(emiss==emiss_find)


Te = (((1-a)*Fin)/(4*sigma))**(1.0/4.0)
Ts =   ((2/(2-emiss)))**(1.0/4.0)*Te
Ta = (  (emiss/(2-emiss)))**(1.0/4.0)*Te
OLR = sigma*Te**4


lapse = Ts - Ta
GHE = sigma*Ts**(4.0) - OLR

scale_height = (287.0*Ta)/9.81
p0 = 1000.
phalf = 500.
z = -scale_height*np.log(phalf/p0)

lapse_rate = (lapse/z)*1000

print "Outgoing longwave radiation is %.2f W m-2" % (OLR)
print "The height corresponding to %.2f hPa is %.2f m" % (phalf, z[eind])
print "The lapse rate corresponding to an atmospheric emissivity of %.2f is %.2f K km-1" % (emiss[eind], lapse_rate[eind])
print "The surface and atmosphere temperatures corresponding to an atmospheric emissivity of %.2f are %.2f K and %.2f K" % (emiss[eind], Ts[eind], Ta[eind])
print "Ta minus Ts and greenhouse effect corresponding to an atmospheric emissivity of %.2f are %.2f K and %.2f W m-2 " % (emiss[eind], lapse[eind], GHE[eind])
print lapse_rate[eind]


golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(emiss,Ts,'r',linewidth=3.0,label=r'$T_s$')
p2, = ax1.plot(emiss,Ta,'b',linewidth=3.0,label=r'$T_a$')
plt.axvline(emiss_obs,color='k',linewidth=3.0,linestyle='dashed',label='Actual emissivity')
#p3, = ax1.plot(emiss,OLR,'g',linewidth=3.0,label='OLR')
ax1.set_ylabel('Temperature (Kelvin)',fontsize=label_fontsize)

plt.grid(True, linestyle='-')
ax2 = ax1.twinx()
if plot_greenhouse_effect == 'false':
    p4 = ax2.plot(emiss,lapse,'r--',linewidth=3.0,label=r'$T_s - T_a$')
    save_name = imagedir + "radiative_equilibrium_model.png"
    ax2.set_ylabel(r'$T_s - T_a$ (Kelvin)',fontsize=label_fontsize)
    legend = ax1.legend(loc='center right', shadow=True)
    legend = ax2.legend(loc='lower left', shadow=True)    
else:
    p4 = ax2.plot(emiss,GHE,'g--',linewidth=3.0,label='GHE')
    save_name = imagedir + "radiative_equilibrium_model_GHE.png"
    ax2.set_ylabel(r'Greenhouse effect (W m$^{-2}$)',fontsize=label_fontsize)
    legend = ax1.legend(loc='lower right', shadow=True)
    legend = ax2.legend(loc='lower left', shadow=True)    
    
ax1.set_ylim([0.,310.])
#ax1.set_xlim([0.2,1.0])

ax1.set_xlabel(r'Emissivity, $\varepsilon$',fontsize=label_fontsize)
#ax2.set_ylabel(r'$T_s - T_a$ (Kelvin)',fontsize=label_fontsize)
#plt.savefig(save_name, bbox_inches='tight')


plt.show()
