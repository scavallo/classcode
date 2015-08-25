# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

# <codecell>

A = (2.53*10**8)*10**3 # convert to Pa
B = 5.42*10**3 # Kelvin
label_fontsize = 14

# <codecell>

T = arange(-40,41,0.25)+273.15 # Convert to Kelvin
es = A*exp(-B/T) # Calculation in SI units
es = es/10**2 # convert to hPa for plotting

# <codecell>

###########################
# Plot
###########################
y0 = np.zeros_like(T).astype('f') # create a vector of (floating point) zeros that is the length of the T vector 
x = T - 273.15 # convert values on x axis to degrees C

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(x,es,'r',linewidth=3.0)
p2, = ax1.plot(x,y0,'k',linewidth=3.0)
ax1.set_xlim([np.min(x),np.max(x)])
ax1.set_ylim([np.min(es),np.max(es)])
ax1.grid(True, linestyle='-')
ax1.set_ylabel('Saturation vapor pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Temperature (degrees C)',fontsize=label_fontsize)


plt.show()

# <codecell>


