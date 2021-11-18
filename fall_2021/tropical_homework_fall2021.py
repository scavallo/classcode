# tropical_homework.py
# 
#
# Steven Cavallo
# University of Oklahoma
# November 2021
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime
import matplotlib
#matplotlib.use('TkAgg')
matplotlib.interactive(True)


import weather_modules as wm
from mstats import *
###########################
# User settings
###########################
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
label_fontsize = 16
H = 10000.
g = 9.81
Rearth = 6.37*10**6
delta_latitude = 0.25;
latitude_beta = 0
omeg_e = (2.*np.pi) / (24.*3600.)
latitudes = np.arange(-90,90+(delta_latitude/2),delta_latitude)

###########################
# End user options
###########################
dlat = np.zeros_like(latitudes).astype('f')
dlat[:] = delta_latitude


# Easiest to define wavelengths on a logarithmic scale here in order to plot the wide spectrum of wavenumbers
'''
Numpy logspace function:
numpy.logspace(start, stop, num=50, endpoint=True, base=10.0, dtype=None)[source]

    Return numbers spaced evenly on a log scale.

    In linear space, the sequence starts at base ** start (base to the power of start) and ends with base ** stop (see endpoint below).
    Parameters:	

    start : float
        base ** start is the starting value of the sequence.
    stop : float
        base ** stop is the final value of the sequence, unless endpoint is False. In that case, num + 1 values are spaced over the interval in log-space, of which all but the last (a sequence of length num) are returned.
    num : integer, optional
        Number of samples to generate. Default is 50.
    endpoint : boolean, optional
        If true, stop is the last sample. Otherwise, it is not included. Default is True.
    base : float, optional
        The base of the log space. The step size between the elements in ln(samples) / ln(base) (or log_base(samples)) is uniform. Default is 10.0.
    dtype : dtype
        The type of the output array. If dtype is not given, infer the data type from the other input arguments.
    
    Returns:	
    samples : ndarray
        num samples, equally spaced on a log scale.
'''
lambdas = np.logspace(-1, 9, num=100, endpoint=True, base=10.0) # Range is from 10^-1 to 10^9 meters with 100 points

k = (2.0*np.pi)/lambdas

f = 2.0*omeg_e*np.sin(latitudes*(np.pi/180))
dy = 111.0*dlat*1000.0
midlatitude_error_mag = 2.0*omeg_e*np.cos(latitudes*(np.pi/180.0))

beta = (2.0*omeg_e*np.cos(latitudes*(np.pi/180)))/Rearth

reflat = 45.
yscale = Rearth*((reflat-latitudes)*(np.pi/180.))
f_full = 2.0*omeg_e*np.sin(latitudes*(np.pi/180.))
f_term1 = np.zeros_like(latitudes).astype('f')
f_term1[:] = 2.0*omeg_e*np.sin(reflat*(np.pi/180.))
f_term2 = ((2.0*omeg_e*np.cos(reflat*(np.pi/180.)))/Rearth)*yscale
f_betaplane = f_term1 + f_term2

betaplane_scale = (yscale/Rearth)*(np.cos(latitudes*(np.pi/180.))/np.sin(latitudes*(np.pi/180.)))

eqind = np.argwhere(latitudes==latitude_beta)  
beta_val = beta[eqind]
print("Beta at %2.f degrees N is %6.3e s-1 m-1" % (latitude_beta, beta_val))

c = np.sqrt(g*H)
fac = 1.0

mu = k*(c/(fac*beta_val))**(1.0/2.0)
nu = np.zeros((4,np.size(k)))
nu_ponc_pos = np.zeros((4,np.size(k)))
nu_ponc_neg = np.zeros((4,np.size(k)))
nu_mrg_pos = np.zeros((4,np.size(k)))
nu_mrg_neg = np.zeros((4,np.size(k)))
for ii in range(0, 4):
    n = ii
    
    omega = (-beta_val*k)/( k**2.0 + (2.0*n+1)*(beta_val/c))
    omega_ponc_pos = np.sqrt(((2.0*n+1)*beta_val*c) + ( (k**2)*(c**2)))
    omega_ponc_neg = -1.0*np.sqrt(((2.0*n+1)*beta_val*c) + ( (k**2)*(c**2)))
    omega_mrg_pos = ((k*c)/2.0)+np.sqrt( ( (k**2.0) * (c**2.0) )/4.0 + beta_val*c)
    omega_mrg_neg = ((k*c)/2.0)-np.sqrt( ( (k**2.0) * (c**2.0) )/4.0 + beta_val*c)
    
    nu[ii,:] = omega/((fac*c*beta_val)**(1.0/2.0))    
    nu_ponc_pos[ii,:] = omega_ponc_pos/((fac*c*beta_val)**(1.0/2.0))     
    nu_ponc_neg[ii,:] = omega_ponc_neg/((fac*c*beta_val)**(1.0/2.0))  
    nu_mrg_pos[ii,:] = omega_mrg_pos/((fac*c*beta_val)**(1.0/2.0))     
    nu_mrg_neg[ii,:] = omega_mrg_neg/((fac*c*beta_val)**(1.0/2.0))  

omega_kelvin = c*k
nu_kelvin =  omega_kelvin/((fac*c*beta_val)**(1.0/2.0))  

xticks = [-90,-60,-30,0,30,60,90]

yticks_beta = np.arange(-3,3.1,1.0)
yticks_beta = yticks_beta*(10**(-11.0))
yticks_beta[3] = 0.
    
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

##################################################################
# Figure 
##################################################################
xticks2 = np.arange(0,4.5,1)
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p00, = ax1.plot(mu[0],np.zeros(np.size(mu[0])),'k',linewidth=3.0)
# Mixed Rossby gravity waves
p0p, = ax1.plot(mu[0],nu_mrg_pos[0,:],'m',linewidth=3.0,label=r'$n=0$ (MRG)')
p0p, = ax1.plot(mu[0],nu_mrg_neg[0,:],'m',linewidth=3.0)
# Equatorial Rossby waves
p1, = ax1.plot(mu[0],nu[1,:],'c',linewidth=2.5,label=r'$n=1,2,3$ (ER)')
p2, = ax1.plot(mu[0],nu[2,:],'c',linewidth=2.5)
p3, = ax1.plot(mu[0],nu[3,:],'c',linewidth=2.5)
# Poincare waves
p1a, = ax1.plot(mu[0],nu_ponc_pos[1,:],'b',linewidth=2.5,label=r'$n=1$ (Poincare)')
p2a, = ax1.plot(mu[0],nu_ponc_pos[2,:],'r',linewidth=2.5,label=r'$n=2$ (Poincare)')
p3a, = ax1.plot(mu[0],nu_ponc_pos[3,:],'g',linewidth=2.5,label=r'$n=3$ (Poincare)')

p1b, = ax1.plot(mu[0],nu_ponc_neg[1,:],'b',linewidth=2.5,label=r'$n=1$ (Poincare)')
p2b, = ax1.plot(mu[0],nu_ponc_neg[2,:],'r',linewidth=2.5,label=r'$n=2$ (Poincare)')
p3b, = ax1.plot(mu[0],nu_ponc_neg[3,:],'g',linewidth=2.5,label=r'$n=3$ (Poincare)')

# Add in your plot for Kelvin waves here


legend = ax1.legend(loc='lower right', shadow=True,fontsize=10)
ax1.set_xlim([0.,4.])
ax1.set_ylim([-4,4])
ax1.xaxis.set_ticks(xticks2)
ax1.set_xticklabels(xticks2)
ax1.set_ylabel(r'$\nu$',fontsize=label_fontsize)
ax1.set_xlabel(r'$\mu$',fontsize=label_fontsize)
save_name = imagedir + "tropical_dispersion.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()







