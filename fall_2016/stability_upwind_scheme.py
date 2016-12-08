
# coding: utf-8

# In[1]:

# Steven Cavallo
# University of Oklahoma
# November 2016
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
wm = imp.load_source('weather_modules','/Users/scavallo/Documents/scripts/python_scripts/weather_modules.py')
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
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

ncints = 21
mu = np.linspace(-2,2,ncints)
kdx = np.linspace(0.0,2.0*np.pi,ncints)
courant = mu   


# In[4]:

amp_upwind = np.empty([np.size(mu), np.size(kdx)])
amp_laxwendroff = np.empty([np.size(mu), np.size(kdx)])
for ii in range (0,np.size(mu)):
    munow = mu[ii]
    for jj in range(0,np.size(kdx)):
        amp_upwind[ii,jj] = np.sqrt(1.0-2.0*munow*(1.0-munow)*(1.0-np.cos(kdx[jj])))
        amp_laxwendroff[ii,jj] = np.sqrt( (1-munow**2.0*(1-np.cos(kdx[jj])))**2.0 + (munow**2.0)*((np.sin(kdx[jj]))**(2.0)))
 
        
cint = 0.1
clevs = np.arange(0,2+(cint/2),cint)
cmap_opt = plt.cm.RdBu_r


ninds = np.where(amp_upwind<clevs[0])
amp_upwind[ninds] = clevs[0]
pinds = np.where(amp_upwind>clevs[-1])
amp_upwind[pinds] = clevs[-1]

ninds = np.where(amp_laxwendroff<clevs[0])
amp_laxwendroff[ninds] = clevs[0]
pinds = np.where(amp_laxwendroff>clevs[-1])
amp_laxwendroff[pinds] = clevs[-1]


# In[5]:

t = np.linspace(0,1,50)
phi_nplus1a = np.zeros_like(t).astype('f')
phi_na = np.zeros_like(t).astype('f')
phi_nplus1b = np.zeros_like(t).astype('f')
phi_nb = np.zeros_like(t).astype('f')

phi_na[:] = 1.0
phi_nb[:] = 1.0
mu1 = 0.5
mu2 = 1.5
for ii in range(1,np.size(phi_na)-1):    
    phi_nplus1a[ii] = (1.0-mu1)*phi_na[ii-1]
    phi_na[ii] = phi_nplus1a[ii]
    
    phi_nplus1b[ii] = (1.0-mu2)*phi_nb[ii-1]
    phi_nb[ii] = phi_nplus1b[ii]    
    


# In[6]:


# Set global figure properties
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden),dpi=128)
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top = 0.93, wspace=0.2, hspace=0.2)


# In[7]:

fig = plt.figure(**figprops)   # New figure   
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
CS1 = plt.contourf(kdx,courant,amp_upwind,cmap=cmap_opt,levels=clevs)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.10) 
CS2 = ax1.contour(kdx,courant,amp_upwind,levels=clevs, colors='k',extend='both',zorder=1,linewidths=2)
plt.clabel(CS2, inline=1, fontsize=10)
plt.xlabel(r'k$\Delta$x',fontsize=label_fontsize)
plt.ylabel('Courant Number',fontsize=label_fontsize)
ax1.set_xlim([0,np.pi])
ax1.set_title('Upstream', fontsize=title_fontsize)
save_name = imagedir + "upstream_stability.png"
#plt.savefig(save_name, bbox_inches='tight')
#plt.show()


# In[8]:

fig = plt.figure(**figprops)   # New figure   
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
#CS = ax1.contour(k,mu,amp,levels=clevs, extend='both',zorder=1,linewidths=2)
#plt.clabel(CS, inline=1, fontsize=10)
CS1 = plt.contourf(kdx,courant,amp_laxwendroff,cmap=cmap_opt,levels=clevs)
#CS1.cmap.set_under('w')
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.10) 
CS2 = ax1.contour(kdx,courant,amp_laxwendroff,levels=clevs, colors='k',extend='both',zorder=1,linewidths=2)
plt.clabel(CS2, inline=1, fontsize=10)
plt.xlabel(r'k$\Delta$x',fontsize=label_fontsize)
plt.ylabel('Courant Number',fontsize=label_fontsize)
ax1.set_xlim([0,np.pi])
ax1.set_title('Lax-Wendroff', fontsize=title_fontsize)
save_name = imagedir + "lax_wendroff_stability.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[9]:

fig = plt.figure(**figprops)   # New figure   
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
CS1 = plt.contour(kdx,courant,amp_upwind,(1.0,),colors='r',extend='both',zorder=1,linewidths=4)
CS2 = plt.contour(kdx,courant,amp_laxwendroff,(1.0,),colors='b',extend='both',zorder=1,linewidths=4)
plt.clabel(CS1, inline=1, fontsize=10)
plt.clabel(CS2, inline=1, fontsize=10)
plt.xlabel(r'k$\Delta$x',fontsize=label_fontsize)
plt.ylabel('Courant Number',fontsize=label_fontsize)
ax1.set_xlim([0,np.pi])
ax1.set_ylim([-5,5])

save_name = imagedir + "upwind_vs_lax_wendroff_stability.png"


# In[10]:

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(t[1:],phi_nplus1a[1:],'r',linewidth=3.0,label=r'$\mu=0.5$')
p2, = ax1.plot(t[1:],phi_nplus1b[1:],'b',linewidth=3.0,label=r'$\mu=1.5$')

plt.show()


# In[ ]:



