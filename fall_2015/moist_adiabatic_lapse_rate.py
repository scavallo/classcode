
# coding: utf-8

# In[ ]:

# moist_adiabatic_lapse_rate.py
# 
# 
# Steven Cavallo
# University of Oklahoma
# September 2015


# In[217]:

import matplotlib.pyplot as plt
import numpy as np
import pylab
import os, datetime

# readsounding is in weather_modules.py
#import weather_modules as wm
import imp
wm = imp.load_source('weather_modules','/Users/scavallo/scripts/python_scripts/weather_modules.py')
sm = imp.load_source('sounding_modules','/Users/scavallo/scripts/python_scripts/sounding_modules.py')
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
from mstats import *


# In[218]:

imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2015/homework/images/'
label_fontsize = 16
plot_logscale = 'false'
g = 9.81 # m s-2
Lv = 2.5*10**6 # J kg-1
cp = 1004 # J K-1 kg-1
Rd = 287. # J K-1 kg-1
Rv = 461. # J K-1 kg-1
epsi = Rd/Rv


# In[219]:

T = np.arange(-50,51,0.25)+273.15 # Convert to Kelvin
es = wm.claus_clap(T)
pres = 1000;  # hPa
#pressures = np.arange(100,1050,100)
#pressures = [50,70,100,250,500,700,850,1000]
pressures = [500]

#ws = wm.satur_mix_ratio(es, pres*100)


# In[220]:

#sat_lapse = (9.81/cp)*( (1 + (Lv*ws)/(Rd*T))) / (1 + (ws*(Lv**2)/(cp*Rv*T**2)) )
#sat_lapse = sat_lapse * 1000. # per km


# In[221]:

###########################
# Plot
###########################
y0 = np.zeros_like(T).astype('f') # create a vector of (floating point) zeros that is the length of the T vector 
y0[:] = g/cp*1000
y00 = np.zeros_like(T).astype('f')

x = T - 273.15 # convert values on x axis to degrees C


golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')

count = 0
nz = len(pressures)
ym = np.zeros_like(pressures).astype('f')
for ii in pressures:
    ws = wm.satur_mix_ratio(es, ii*100)  
    sat_lapse = (9.81/cp)*( (1 + (Lv*ws)/(Rd*T))) / (1 + (ws*(Lv**2)/(cp*Rv*T**2)) ) 
    sat_lapse = sat_lapse*1000
    
    ym[count] = np.median(sat_lapse)
    
    numeric_now = float(count)/float(nz)
    colnow = [numeric_now,numeric_now,numeric_now]
    
    pref = r'$\Gamma_s$ at '
    pref2 = r'$\Gamma_d$'

    textstrnow = pref + np.str(ii) + ' hPa'
    p1, = ax1.plot(sat_lapse,x, color=colnow,linewidth=3.0,label = textstrnow )
    p0, = ax1.plot(y0,x,'b',linewidth=3.0, label = pref2)
    count += 1

y00[:] = np.mean(ym)
print 'The mean saturated adiabatic lapse rate is %5.2f K km-1'  %np.mean(ym)
pref3 = r'Mean $\Gamma_s$ '
p00, = ax1.plot(y00,x,'r',linewidth=3.0, label=pref3)    
ax1.set_xlim([1,11])
ax1.set_ylim([np.min(x),np.max(x)])
plt.gca().invert_yaxis()
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower right', shadow=True)
ax1.set_xlabel(r'Saturated adiabatic lapse rate (K km$^{-1}$)',fontsize=label_fontsize)
ax1.set_ylabel(r'Temperature ($^{\circ}$C)',fontsize=label_fontsize)
save_name = imagedir + "hw3_saturated_adiabatic_lapes_rates.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()


# In[ ]:



