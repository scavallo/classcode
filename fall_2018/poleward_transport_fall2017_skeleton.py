
# coding: utf-8
#
# Version info
# >>conda -V
# conda 4.3.21
# >>python --version
# Python 3.5.2 :: Anaconda custom (x86_64)

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

# Make sure you link to these in your working directory!  They are available in the 'utils' directory on github:
# https://github.com/scavallo/classcode/tree/master/utils
import utilities_modules as um
from mstats import *


# In[2]:

###########################
# User settings
# YOU need to set the path to the location of your data input file (datadir) and image save directory (imagedir) on your local machine
###########################
datadir = '/'
imagedir = '/'
infile = 'slp_6090_annual_1948_2017.dat'

label_fontsize = 16
title_fontsize = 18
latitude = 60.

conversion_factor = 1000. # for converting meters to mm
climo_years = [1981,2010]
analysis_years1 = [1996,1996]
analysis_years2 = [1997,1997]

ylims1 = [1005, 1020]
ylims2 =  [-6, 6]

month_labels = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
ndays = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]


# In[3]:

###########################
# Constants
###########################
R_e = 6.37*10**6.0 # radius of Earth in meters
coslatr = np.cos(latitude*np.pi/180.) 
sinlatr = np.sin(latitude*np.pi/180.)
p0 = 100000.  # Pascals


# In[4]:

# Reading in the file
ab = np.loadtxt(datadir+infile, skiprows=1)       
years = ab[:,0]
slp_jan = ab[:,1]
slp_july = ab[:,7]
slp_full = ab[:,1:]

# In the data file, you will notice that -999.999 implies NaN.  The below filters those out.
slp_full = um.filter_numeric_nans(slp_full,0,float('NaN'),'low')


# In[5]:

# Extract the subsets of the data that you need
[ind1a] = np.where(years == analysis_years1[0])
[ind1b] = np.where(years == analysis_years1[1])
[ind2a] = np.where(years == analysis_years2[0])
[ind2b] = np.where(years == analysis_years2[1])

[climind1a] = np.where(years == climo_years[0])
[climind1b] = np.where(years == climo_years[1])

slp_subset1 = np.nanmean(slp_full[ind1a:ind1b+1,:],0)
slp_climatology = np.nanmean(slp_full[climind1a:climind1b+1,:],0)


# In[6]:

######################
# Plot
######################
t = np.arange(1,13,1) 
z0 = np.zeros(len(t))

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(t,slp_climatology,'b--',linewidth=3.0,label='Monthly mean Arctic SLP climatology')
ax1.set_xlim([np.min(t),np.max(t)])
ax1.set_ylim([ylims1[0],ylims1[-1]])
ax1.xaxis.set_ticks(t[:])
ax1.set_xticklabels(month_labels[:],rotation=90)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)

ax1.set_title('Monthly-mean sea level pressure climatology over the Arctic',fontsize=title_fontsize)
ax1.set_ylabel('hPa',fontsize=label_fontsize)
ax1.set_xlabel('Month',fontsize=label_fontsize)

# Comment the line below back in if you want to save the plot
#plt.savefig(imagedir + 'unique_name_of_your_figure' + '.png', bbox_inches='tight')

# Need the below if you want to display the plot to your screen
plt.show()

