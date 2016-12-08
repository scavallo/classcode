
# coding: utf-8

# In[1]:

# Steven Cavallo
# University of Oklahoma
# December 2016
###########################
# imports
###########################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
import pylab
import os, datetime
import imp

#wm = imp.load_source('weather_modules','/Users/scavallo/Documents/scripts/python_scripts/weather_modules.py')
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
import weather_modules as wm
import utilities_modules as um
from mstats import *


# In[2]:

###########################
# User settings
###########################
label_fontsize = 16
title_fontsize = 24
#imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2015/homework/images/'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'


# In[3]:

ncints = 100
a = np.linspace(-2,2,ncints)
b = np.linspace(-2,2,ncints)
deltat = 1.0

x = a*deltat
y = b*deltat
[X,Y] = np.meshgrid(x,y)


# In[4]:

# Backward Euler
z = 1.0/((1-X)**(2.0) + (Y)**(2.0) )


# In[5]:

# Set global figure properties
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden),dpi=128)
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top = 0.93, wspace=0.2, hspace=0.2)

cint = 0.1
clevs = np.arange(0.0,2.0+(cint/2),cint)
cmap_opt = plt.cm.RdBu_r

z = um.filter_numeric_nans(z,np.max(clevs),np.max(clevs),'high')


# In[6]:

fig = plt.figure(**figprops)   # New figure   
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])

CS1 = plt.contourf(X,Y,z,cmap=cmap_opt,levels=clevs)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.10)


CS2 = ax1.contour(X,Y,z,levels=clevs, colors='k',extend='both',zorder=1,linewidths=2)
plt.clabel(CS2, inline=1, fontsize=10)

CS3 = plt.contour(X,Y,z,(1.0,),colors='r',extend='both',zorder=1,linewidths=4)
plt.xlabel(r'a$\Delta$t',fontsize=label_fontsize)
plt.ylabel(r'b$\Delta$t',fontsize=label_fontsize)

plt.show()

