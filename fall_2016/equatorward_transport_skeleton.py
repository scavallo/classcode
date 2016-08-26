###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

from mstats import *

###########################
# Constants
###########################
R_e = 6.37*10**6.0 # radius of Earth in meters
coslatr = np.cos(latitude*np.pi/180.) 
p0 = 100000.  # Pascals

###########################
# How to read in a file and extract the years of analysis
###########################
ab = np.loadtxt(datadir+infile1, skiprows=1)       
years_april = ab[:,0]
slp_april = ab[:,1]

inds = np.where( (years_april>=analysis_years[0]) & (years_april<=analysis_years[1]))
tmp = slp_april[inds]
del slp_april
slp_april = tmp
del tmp

years = years_april[inds]

# to check statistics, make sure mstats.m is in this directory and type:
mstats(slp_april)

###########################
# If you want to calculate a running mean, do something like this:
###########################
#new_array = np.zeros(len(slp_april)) # create a new array with all zeros
#new_array[0:nyears] = float('NaN') # change to zeros in the years where you can't compute a 5-y mean to NaNs
#for ii in range(nyears,len(new_array)):
#        new_array[ii] = np.nanmean(slp_april[ii:ii+nyears+1])

###########################
# To calculate the long-term climatology, try something like:
###########################
#inds_climo = np.where( (years>=climo_years[0]) & (years<=climo_years[1]))
#myarray_climo = np.nanmean(myarray[inds_climo])

#myanomaly = myarray - myarray_climo # x = xbar + xprime --> xprime = x - xbar

###########################
# Plot example
###########################
# Define your x-axis
t = years

# These lines only need to be done before the first graph
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

# Below is unique to each graph
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,myarray,'r',linewidth=3.0,label='My label')
p2, = ax1.plot(t,myotherarray,'b',linewidth=3.0,label='My other label')
ax1.set_xlim([np.min(years),np.max(years)])
ax1.set_ylim([ylims1[0],ylims1[-1]])

ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower left', shadow=True)

ax1.set_title('My title',fontsize=title_fontsize)
ax1.set_ylabel('hPa',fontsize=label_fontsize)
ax1.set_xlabel('Year',fontsize=label_fontsize)
figname1 = 'myfigurename.png'
plt.savefig(imagedir + figname1, bbox_inches='tight')