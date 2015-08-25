# tropical_homework.py
# 
#
# Steven Cavallo
# University of Oklahoma
# October 2014
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

# readsounding is in weather_modules.py
import weather_modules as wm
# The module below will tell me statistics of an array by typing:
#    mstats(array_name) 
from mstats import *
###########################
# User settings
###########################
datadir = ''
imagedir = ''
label_fontsize = 16
H = 10000.
g = 9.81
delta_latitude = 0.5;
latitude_beta = 0
omeg_e = (2*np.pi) / (24*3600)
latitudes = np.arange(-90,90+(delta_latitude/2),delta_latitude)

###########################
# End user options
###########################

dlat = np.zeros_like(latitudes).astype('f')
dlat[1:-1] = (latitudes[2:] - latitudes[:-2])/2
dlat[0] = (latitudes[1]-latitudes[0]) 
dlat[-1] = (latitudes[-1]-latitudes[-2])

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

k = (2*np.pi)/lambdas

f = 2.0*omeg_e*np.sin(latitudes*(np.pi/180))
dy = 111.0*dlat*1000

beta = np.zeros_like(latitudes).astype('f')
beta[1:-1] = (f[2:] - f[:-2])/(2*dy[1:-1])
beta[0] = (f[1] - f[0])/(dy[1]) 
beta[-1] = (f[-1] - f[-2])/(dy[-1]) 

eqind = np.argwhere(latitudes==latitude_beta)  
beta_val = beta[eqind]
print "Beta at the %2.f degrees N is %6.3e s-1 m-1" % (latitude_beta, beta_val)

c = np.sqrt(g*H)

mu = k*(c/(2.0*beta_val))**(1.0/2.0)
nu = np.zeros((4,np.size(k)))
nu_ponc_pos = np.zeros((4,np.size(k)))
nu_ponc_neg = np.zeros((4,np.size(k)))
for ii in range(0, 4):
    n = ii
    
    omega = (-beta_val*k)/( k**2.0 + (2.0*n+1)*(beta_val/c))
    omega_ponc_pos = np.sqrt(((2.0*n+1)*beta_val*c) + ( (k**2)*(c**2)))
    omega_ponc_neg = -1.0*np.sqrt(((2.0*n+1)*beta_val*c) + ( (k**2)*(c**2)))
    
    nu[ii,:] = omega/((2.0*c*beta_val)**(1.0/2.0))    
    nu_ponc_pos[ii,:] = omega_ponc_pos/((2.0*c*beta_val)**(1.0/2.0))     
    nu_ponc_neg[ii,:] = omega_ponc_neg/((2.0*c*beta_val)**(1.0/2.0))  


omega_kelvin = c*k
nu_kelvin =  omega_kelvin/((2.0*c*beta_val)**(1.0/2.0))  

xticks = np.arange(-90,91,10)

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(latitudes,f,'r',linewidth=3.0,label='f')

ax2 = ax1.twinx()
p2, = ax2.plot(latitudes,beta,'b',linewidth=3.0,label=r'$\beta$')

ax1.grid(linestyle='-')

ax1.set_xlim([-90.,90.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_ylabel('Coriolis parameter (f)',fontsize=label_fontsize)
ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)
ax2.set_ylabel(r'$\beta = \partial{f}/\partial{y}$',fontsize=label_fontsize)

legend = ax1.legend(loc='upper left', shadow=True)
legend = ax2.legend(loc='lower right', shadow=True)
save_name = imagedir + "tropical_f_vs_beta.png"
plt.savefig(save_name, bbox_inches='tight')


xticks2 = np.arange(0,4.5,1)
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p00, = ax1.plot(mu[0],np.zeros(np.size(mu[0])),'k',linewidth=3.0)
#p0, = ax1.plot(mu[0],nu[0,:],'k--',linewidth=3.0,label=r'$n=0$')
p1, = ax1.plot(mu[0],nu[1,:],'b',linewidth=2.5,label=r'$n=1$')
p2, = ax1.plot(mu[0],nu[2,:],'r',linewidth=2.5,label=r'$n=2$')
p3, = ax1.plot(mu[0],nu[3,:],'g',linewidth=2.5,label=r'$n=3$')

p1a, = ax1.plot(mu[0],nu_ponc_pos[1,:],'b',linewidth=2.5)
p2a, = ax1.plot(mu[0],nu_ponc_pos[2,:],'r',linewidth=2.5)
p3a, = ax1.plot(mu[0],nu_ponc_pos[3,:],'g',linewidth=2.5)

p1b, = ax1.plot(mu[0],nu_ponc_neg[1,:],'b',linewidth=2.5)
p2b, = ax1.plot(mu[0],nu_ponc_neg[2,:],'r',linewidth=2.5)
p3b, = ax1.plot(mu[0],nu_ponc_neg[3,:],'g',linewidth=2.5)

p4 = ax1.plot(mu[0],nu_kelvin[0],'k--',linewidth=2.5, label='Kelvin wave')

legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_xlim([0.,4.])
ax1.set_ylim([-3,3])
ax1.xaxis.set_ticks(xticks2)
ax1.set_xticklabels(xticks2)
ax1.set_ylabel(r'$\nu$',fontsize=label_fontsize)
ax1.set_xlabel(r'$\mu$',fontsize=label_fontsize)
save_name = imagedir + "tropical_dispersion.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()





