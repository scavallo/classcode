# mixed_layer_model
# 
#
# Steven Cavallo
# University of Oklahoma
# November 2021
###########################
# imports
###########################
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import sys
import pylab
import os, datetime

from IPython.display import display,clear_output
from IPython.display import Image
import time as Time

###########################
# User settings
###########################
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
label_fontsize = 16
label_fontsize_subplot = 8
figname1 = "abl2021_vertical_profile_firstlast_comparison_part1"
figname2 = "abl2021_theta2m_thetasst_part1"
figname3 = "abl2021_heatflux_part1"
figname4 = "abl2021_entrainment_velocity_part1"
figname5 = "abl2021_pblheight_thetabar_part1"

#######################
# Frequently changed parameters 
#######################
roughness_length = 0.005
zi0 = 500. # initial boundary layer height
diurnal_range = 15. # diurnal range in Kelvin
theta_min = 280.  # daily minimum potential temperature
flux_option = 1 # 1 to compute flux based on 2-m and skin temperature.  Set sst_inflation below.
                # 2 to compute skin temp from flux.  Set max_sensible_heatflux below.
                # 3 for a constant skin temperature.  Set below in sst_tempk_const 
sst_inflation = 2.8 # inflation factor for surface skin temperature; Only used if flux_option == 1; initial setting = 2.8
max_sensible_heatflux = 400. # W m-2; Only used if flux_option == 2
sst_tempk_const = 273. # Kelvin; Only used for flux_option == 3

Afac = 0.2 # Mixing parameter
U10 = 5.0 # 10-meter wind speed magnitude
capping_strength = 10.0 # capping strength, in Kelvin
dthetadz = 5.0 # theta lapse rate of free troposphere, in Kelvin km-1

#######################
# Less frequently changed parameters 
#######################
T = 86400. # Period of simulation in seconds.  
deltat = 60. # timestep increment in seconds
deltaz = 10. # vertical grid spacing in meters

#######################
# Plot settings
#######################
xlims = [theta_min-5,330.]
ylims = [0.,6000.] # in meters

#######################
# Constants.  Do not change!
#######################
cp = 1004.
rho = 1.0
Rd = 287.
g = 9.81
k = 0.4
p0 = 101325.

coef1 = 0.147
coef2 = 0.0004

# What is the maximum physically-allowable entrainment velocity?  
#     we = Delta z_i / Delta t.
#     In 12 hours (0.5 days), Delta z_i should not be more than the depth of the troposphere, or 10 km:
#     ==> we ~ (10 km) / (0.5 days) = 0.23
max_we = 0.23 # Maximum allowable entrainment velocities.  Typically don't even get above 0.1 m s-1 (Angevine et al. 1998).
              


#######################
# Calculations
#######################
CD = k**2.0/(np.log(10.0/np.array(roughness_length)))**2.0
print('The drag coefficient CD is %8.4f' %(CD))
CH = coef1*CD + coef2
print('The heat coefficient CH is %8.4f' %(CH))

# Define times and vertical levels
t = np.arange(0,86401.0,deltat)
zlevs = np.arange(0.,10000.,deltaz)
thetavert = np.zeros_like(zlevs).astype('f')

c = 2.0*np.pi/T
theta_max = theta_min + diurnal_range
theta_avg = (theta_max+theta_min)/2.0
theta_2m = (diurnal_range/2.0)*np.sin(c*t)+theta_avg
if flux_option == 1:
    theta_2m = ((diurnal_range)/2.0)*np.sin((c*t)-np.pi/4.0)+theta_avg
    theta_0m = ((diurnal_range*sst_inflation)/2.0)*np.sin((c*t))+theta_avg
    dtheta = theta_2m - theta_0m
    FHS = -1.0*CH*U10*dtheta
    QHS = rho*cp*FHS
elif flux_option == 2:
    QHS = max_sensible_heatflux*np.sin(c*t)
    FHS = QHS/(rho*cp)
    theta_0m = theta_2m + (FHS/(CH*U10))
elif flux_option == 3:
    theta_2m = ((diurnal_range)/2.0)*np.sin((c*t)-np.pi/4.0)+theta_avg
    theta_0m = np.zeros_like(theta_2m)
    theta_0m[:] = sst_tempk_const
    dtheta = theta_2m - theta_0m
    FHS = -1.0*CH*U10*dtheta
    QHS = rho*cp*FHS    
 


#######################
# Set up animation
#######################

nplots = 13
plottimes = np.round(np.linspace(0,len(FHS)+1,nplots))

we = np.zeros_like(FHS).astype('f')
thetabar = np.zeros_like(FHS).astype('f')
zi = np.zeros_like(FHS).astype('f')
FHzi = np.zeros_like(FHS).astype('f')
cap = np.zeros_like(FHS).astype('f')

thetabar[0] = theta_2m[0]
zi[0] = zi0
thetafree = thetabar[0] + capping_strength


#######################
# Iterate though the day
#######################

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
    # then A must be set to positive when FHS>0 and A must be set to negative when FHS<0
    if ( FHS[count] >= 0 ):
        A = 1.0*Afac
    else:
        A = -1.0*Afac       
    
    capping_strength = thetafree - thetabar[count]   
    # Superadiabatic PBL not supported with this model
    if capping_strength < 0:
        capping_strength = 0.0 
    we[count] = (A*FHS[count])/capping_strength
    if (we[count]>=max_we):
        print('Warning!  High entrainment velocities of %8.5f m s-1' %(we[count]))
        we[count] = max_we
        dthetabar_prev = dthetabar
         
    FHzi[count] = -1.0*A*FHS[count]
    cap[count] = capping_strength
    
    if count < 0:
        print('%d, %8.5f, %8.5f, %8.5f, %8.5f' %(count, A, FHS[count],we[count],capping_strength))
    
    dthetabar = ((1.0+A)*FHS[count])/zi[count]
    if (we[count]>=max_we):
        dthetabar = dthetabar_prev   
        we[count] = max_we
    
    if ( count < len(FHS)-1 ):       
        thetabar[count+1] = thetabar[count] + (deltat*dthetabar)
        zi[count+1] = zi[count] + (deltat*we[count])

    [zinds] = np.where(zlevs>zi[count])
    try:
        pblind = zinds[0]
    except:
        pblind = len(zlevs)-1 

    thetavert[0:pblind] = thetabar[count]
    thetavert[pblind] = thetabar[count]+capping_strength
    if plotcount == 0:
        thetavert[pblind] = thetabar[count]+capping_strength
    else:
        thetavert[pblind] = thetavert_init[pblind]
        
    thetafree = thetavert[pblind]
    for kk in range(pblind+1,len(zlevs)-1):
        thetavert[kk] = thetavert[kk-1]+(dthetadz*10**(-3.0))*deltaz
        
        
    if ( plotcount == 0 ):
        thetavert_init = np.zeros_like(zlevs).astype('f')
        thetavert_init[:] = thetavert[:]
    if thrs == 3:
        thetavert_save_3h = np.zeros_like(zlevs).astype('f')
        thetavert_save_3h[:] = thetavert[:]    
    if thrs == 6:
        thetavert_save_6h = np.zeros_like(zlevs).astype('f')
        thetavert_save_6h[:] = thetavert[:]
    if thrs == 9:
        thetavert_save_9h = np.zeros_like(zlevs).astype('f')
        thetavert_save_9h[:] = thetavert[:]
    if thrs == 12:
        thetavert_save_12h = np.zeros_like(zlevs).astype('f')
        thetavert_save_12h[:] = thetavert[:]
    
    if ( (count == plottime_now) or (count == len(FHS)-1)): 
        Time.sleep(1./2.)
        ax.clear()
        labeltext = (r'$\theta$ at t = %7.2f hours' %(thrs))
        ax.plot(thetavert,zlevs,'r-',linewidth=3.0,label=labeltext)
        clear_output(wait=True)
        ax.legend(loc='upper right',shadow=True,fontsize=12)
        ax.set_xlabel(r'$\overline{\theta}$ (K)',fontsize=label_fontsize,fontweight='bold')
        ax.set_ylabel('Height (meters)',fontsize=label_fontsize,fontweight='bold')
        ax.grid(True, linestyle='-')
        ax.set_xlim((xlims[0],xlims[1]))
        ax.set_ylim((ylims[0],ylims[1]))
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
    


#######################
# Set up plots
#######################
tsunset_ind = int(np.round(len(t)/2.0))
print('Boundary layer height at sunset is %7.2f m' %(zi[tsunset_ind]))
print('Boundary layer theta at sunset is %7.2f K' %(thetabar[tsunset_ind]))

thrs = t/(60.0*60.0)
xticks = np.arange(0,25,6)
zeroarr = np.zeros_like(thrs).astype('f')     

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 


#######################
# Figure 1
#######################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thetavert_init,zlevs,'r',linewidth=3.0,label=r'$\theta_{t=0}$')
p2, = ax1.plot(thetavert_save_3h,zlevs,'m--',linewidth=3.0)
p3, = ax1.plot(thetavert_save_6h,zlevs,'m--',linewidth=3.0)
p4, = ax1.plot(thetavert_save_9h,zlevs,'m--',linewidth=3.0)
p5, = ax1.plot(thetavert_save_12h,zlevs,'r--',linewidth=3.0,label=r'$\theta_{t=12 \, hrs}$')
ax1.legend(loc='upper right',shadow=True,fontsize=12)
ax1.set_xlabel(r'$\overline{\theta}$ (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_ylabel('Height (meters)',fontsize=label_fontsize,fontweight='bold')
ax1.grid(True, linestyle='-')
ax1.set_xlim((xlims[0],xlims[1]))
#ax1.set_ylim((0.,6000.))
ax1.set_ylim((ylims[0],ylims[1]))

save_name = imagedir + figname1
plt.savefig(save_name, bbox_inches='tight')


#######################
# Figure 2
#######################
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
ax1.set_ylabel(r'$\theta$ (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname2
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


#######################
# Figure 3
#######################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thrs,rho*cp*FHS,'r',linewidth=3.0,label=r'$Q_{HS}$')
p2, = ax1.plot(thrs,rho*cp*FHzi,'b',linewidth=3.0,label=r'$Q_{Hzi}$')
legend = ax1.legend(loc='upper right', shadow=True)

ax1.grid(True, linestyle='-')
ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_ylabel('Sensible heat flux (W m-2)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname3
plt.savefig(save_name, bbox_inches='tight')


#######################
# Figure 4
#######################
yticks_ax1 = np.arange(0,0.25,0.05)
nticks = len(yticks_ax1)
yticks_ax2 = np.linspace(0,20,nticks)

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
ax2 = ax1.twinx()
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(thrs,we,'r',linewidth=3.0,label=r'$w_e$')
p2, = ax2.plot(thrs,cap,'b',linewidth=3.0,label=r'$\Delta\theta$')
legend = ax1.legend(loc='upper right', shadow=True)
legend = ax2.legend(loc='lower right', shadow=True)

ax1.grid(True, linestyle='-')


ax1.set_yticks(yticks_ax1)
ax1.set_yticklabels(yticks_ax1,fontsize=16)
ax2.set_yticks(yticks_ax2)
ax2.set_yticklabels(yticks_ax2,fontsize=16)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%4.2f'))

ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_ylabel('Entrainment velocity',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax2.set_ylabel(r'$\Delta\theta$ (K)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname4
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


#######################
# Figure 5
#######################
zimax = np.ceil(np.max(zi))
ziceil = zimax
for ii in range(0,1000):
    if (np.mod(ziceil,1000) == 0):
        ymaxval = ziceil
        break
    ziceil = ziceil+1

    
yticks_ax1 = np.arange(zi0,ziceil,500)
nticks = len(yticks_ax1)
yticks_ax2 = np.linspace(theta_min-5,312,nticks)

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
#p0, = ax1.plot(CH,y0,'k',linewidth=3.0)

ax1.grid(True, linestyle='-')


ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

ax1.set_ylim((ylims[0],ylims[1]))
#ax1.set_yticks(yticks_ax1)
#ax1.set_yticklabels(yticks_ax1,fontsize=16)
#ax2.set_yticks(yticks_ax2)
#ax2.set_yticklabels(yticks_ax2,fontsize=16)

ax1.set_ylabel('PBL height (m)',fontsize=label_fontsize,fontweight='bold')
ax2.set_ylabel(r'$\overline{\theta}$ (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + figname5
plt.savefig(save_name, bbox_inches='tight')
plt.show()






