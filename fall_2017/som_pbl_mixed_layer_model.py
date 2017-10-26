# som_pbl_mixed_layer_model
# 
#
# Steven Cavallo
# University of Oklahoma
# October 2017
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

from IPython.display import display,clear_output
from IPython.display import Image
import time as Time

###########################
# Set directories, filenames
###########################
datadir = '/Users/scavallo/data/'
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
label_fontsize = 16
label_fontsize_subplot = 8
figname1 = "abl2017_fig1"
figname2 = "abl2017_fig2"
figname3 = "abl2017_fig3"
figname4 = "abl2017_fig4"


###########################
# Parameters to set
###########################
CD = 0.06   # Drag coefficient corresponding to surface roughness about 2.  See Table 9.2 in W&H for various values.
zi0 = 1000. # initial boundary layer height
diurnal_range = 30. # diurnal temperature range
theta0_min = 280.  # daily minimum potential temperature
flux_sponge = 0.5 # Set to a value between 0 and 1.  0 for the largest surface heat fluxes, 1 for none.

T = 86400. # Period of simulation in seconds.  Have not tested changing this, so don't!
deltat = 60. # timestep increment in seconds.  Have not tested changing this, so don't!
deltaz = 100. # vertical grid spacing in meters

Afac = 0.2 # Mixing parameter, A
U10 = 5.0 # 10-m wind speed magnitude in m s-1
capping_strength = 10.0 # Strength of capping inversion, in Kelvin
dthetadz = 5.0 # theta lapse rate, in Kelvin km-1
max_we = 0.5 # Maximum allowable entrainment velocities.   Typically don't even get above 0.1 m s-1 (Angevine et al. 1998). 
###########################
# Constants.  Do not change!
###########################
cp = 1004.
rho = 1.0
Rd = 287.
g = 9.81
k = 0.4
p0 = 101325.

###########################
# Begin calculations
###########################
# An empirical relation to get from drag coefficient to heat coefficient
CH = 0.147*CD + 0.0004
print('The heat coefficient is CH = %8.4f' %(CH))

# Time array.  Not a good idea to change unless you know what you're doing.
t = np.arange(0,86401.0,deltat)
# Vertical grid
zlevs = np.arange(0.,10000.,deltaz)
thetavert = np.zeros_like(zlevs).astype('f')

if (flux_sponge < 0):
    flux_sponge = 0.0
if (flux_sponge > 1):
    flux_sponge = 1.0

# Assume the surface heats up a lot more than the first atmospheric layer above it.
# Functions below assume t=0 is sunrise
c = 2.0*np.pi/T
theta0_max = theta0_min + diurnal_range
theta0_avg = (theta0_max+theta0_min)/2.0
theta_0m = (diurnal_range/2.0)*np.sin(c*t )+theta0_avg
theta_2m = ((diurnal_range*0.5)/2.0)*np.sin(c*t)+theta0_avg

# Compute kinematic heat flux
dtheta = theta_2m - theta_0m
FHS = -1.0*CH*U10*dtheta

###########################
# Some initializations 
###########################
nplots = 13
plottimes = np.round(np.linspace(0,len(FHS)+1,nplots))

we = np.zeros_like(FHS).astype('f')
thetabar = np.zeros_like(FHS).astype('f')
zi = np.zeros_like(FHS).astype('f')


thetabar[0] = theta_2m[0]
zi[0] = zi0
thetafree = thetabar[0] + capping_strength
print(thetabar[0],thetafree)

###########################
# Animation
###########################
count = 0
plotcount = 0
plottime_now = plottimes[plotcount]

f, ax =plt.subplots(figsize=(8,6))
display(f)

while ( count<=len(FHS)-1 ):
    thrs = (t[count])/3600.
    #print('Hours: %7.2f' %(thrs))
   
    # Note we are not parameterizing FHzi.  The parameter to account for FHzi is A.
    # Since A=-FHzi/FHS, and in this model FHzi is always directed downward (negative),
    # then A must be positive when FHS>0 and A must be negative when FHS<0
    if ( FHS[count] >= 0 ):
        A = 1.0*Afac
    else:
        A = -1.0*Afac
    
    capping_strength = thetafree - thetabar[count]  
    # Superadiabatic PBL top not allowed
    if capping_strength < 0:
        capping_strength = 0.0     
    we[count] = (A*FHS[count])/capping_strength
    # Model breaks down when capping inversion is too weak. Let it finish but give a warning.
    if (we[count]>=max_we):
        print('Warning!  High entrainment velocities of %8.5f m s-1' %(we[count]))
        we[count] = max_we
        dthetabar_prev = dthetabar
        
    if count < 0:
        print('%d, %8.5f, %8.5f, %8.5f, %8.5f' %(count, A, FHS[count],we[count],capping_strength))
    dthetabar = ((1.0+A)*FHS[count])/zi[count]
    # If entrainment velocities too high, keep tendencies constant from last time step
    if (we[count]>=max_we):
        dthetabar = dthetabar_prev    
        
    # Numerical integration to next time step
    if ( count < len(FHS)-1 ):

        thetabar[count+1] = thetabar[count] + (deltat*dthetabar)
        zi[count+1] = zi[count] + (deltat*we[count])

    if (we[count]>=max_we):
        [zinds] = np.where(thetavert>thetabar[count])
        try:
            pblind = zinds[0]
        except:
            pblind = len(zlevs)-1
        zi[count+1] = zlevs[pblind-1]
    else:
        [zinds] = np.where(zlevs>zi[count])
        try:
            pblind = zinds[0]
        except:
            pblind = len(zlevs)-1 

    # Only save the animation every so often, say every hour or two
    if ( (count == plottime_now) or (count == len(FHS)-1)):
    
        # These are just for illustration in the figure and are not used in for calculations.
        if (we[count]>=max_we):
            [zinds] = np.where(thetavert>thetabar[count])
            try:
                pblind = zinds[0]
            except:
                pblind = len(zlevs)-1
        else:
            [zinds] = np.where(zlevs>zi[count])
            try:
                pblind = zinds[0]
            except:
            	pblind = len(zlevs)-1
        thetavert[0:pblind] = thetabar[count]
        thetavert[pblind] = thetabar[count]+capping_strength
        for kk in range(pblind+1,len(zlevs)-1):
            thetavert[kk] = thetavert[kk-1]+(dthetadz*10**(-3.0))*deltaz       

 
        # If you want to speed up the rate of the plots, increase the denominator in the line below
        Time.sleep(1./2.)
        ax.clear()
        labeltext = (r'$\theta$ at t = %7.2f hours' %(thrs))
        ax.plot(thetavert,zlevs,'r-',linewidth=3.0,label=labeltext)
        clear_output(wait=True)
        ax.legend(loc='upper right',shadow=True,fontsize=12)
        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel('z',rotation=0)
        ax.grid(True, linestyle='-')
        ax.set_xlim((290.,330.))
        ax.set_ylim((0.,5000.))
        display(f)
        if (plotcount < 10):
            figname_now = 'animation_0' + str(plotcount) + '.png'
        elif ( (plotcount >= 10) or (plotcount <100) ):
            figname_now = 'animation_' + str(plotcount) + '.png'
        print(figname_now)
        save_name = imagedir + figname_now
        plt.savefig(save_name, bbox_inches='tight')   
        
        plotcount += 1
        if (count < len(FHS)-1):
            plottime_now = plottimes[plotcount]
 

        
    count+=1
    

###########################
# Set up plots
###########################
thrs = t/(60.0*60.0)
xticks = np.arange(0,25,6)
zeroarr = np.zeros_like(thrs).astype('f')     

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

###########################
# Figure 1
###########################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thrs,theta_0m,'r',linewidth=3.0,label=r'$\theta_{sst}$')
p2, = ax1.plot(thrs,theta_2m,'b',linewidth=3.0,label=r'$\theta_{2m}$')
legend = ax1.legend(loc='upper right', shadow=True)

ax1.grid(True, linestyle='-')
ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_ylabel('Surface skin temperature (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname1
plt.savefig(save_name, bbox_inches='tight')

###########################
# Figure 2
###########################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thrs,we,'r',linewidth=3.0,label=r'$w_e$')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.grid(True, linestyle='-')

ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_ylabel('Entrainment velocity',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname2
plt.savefig(save_name, bbox_inches='tight')

###########################
# Figure 3
###########################
zimax = np.ceil(np.max(zi))
ziceil = zimax
for ii in range(0,1000):
    if (np.mod(ziceil,1000) == 0):
        ymaxval = ziceil
        break
    ziceil = ziceil+1

    
yticks_ax1 = np.arange(zi0,ziceil,500)
nticks = len(yticks_ax1)
yticks_ax2 = np.linspace(290,312,nticks)

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
ax2 = ax1.twinx()
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thrs,zi,'b',linewidth=3.0,label=r'$z_i$')
p2, = ax2.plot(thrs,thetabar,'r',linewidth=3.0,label=r'$\overline{\theta}$')
legend = ax2.legend(loc='lower right', shadow=True)
legend = ax1.legend(loc='upper left', shadow=True)

ax1.grid(True, linestyle='-')


ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

#ax1.set_yticks(yticks_ax1)
#ax1.set_yticklabels(yticks_ax1,fontsize=16)
#ax2.set_yticks(yticks_ax2)
#ax2.set_yticklabels(yticks_ax2,fontsize=16)

ax1.set_ylabel('PBL height (m)',fontsize=label_fontsize,fontweight='bold')
ax2.set_ylabel(r'$\overline{\theta}$ (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname3
plt.savefig(save_name, bbox_inches='tight')
plt.show()






