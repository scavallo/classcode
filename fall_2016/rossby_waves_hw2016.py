# Rossby_waves_hw2016.py
#
# Script to plot Rossby wave phase speed and group velocity, given a background flow, wavenumber, and latitude
#
# Steven Cavallo
# University of Oklahoma
# September 2016

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab

###########################
# User parameters
###########################
umax = 20.0
jet_option = 3 # 1 for no jet/ubar=0, 2 for constant jet/ubar=constant=umax, 3 for zonally varying jet with umax amp
jet_thinning_factor = 1.0
beta_plane = 'False'
constant_wavenumber = 'False'
num_waves_merid = 0.0
num_waves_zonal = 5.0
delta_latitude = 0.5;
lat_info = 45.0

label_fontsize = 16
imagedir = '/Users/scavallo/scripts/python_scripts/images/'

###########################
# Be careful before changing the rest...
###########################
x = np.linspace(0,2*np.pi,100)
omeg_e = (2*np.pi) / (24*3600)
latitudes = np.arange(-90,90+(delta_latitude/2),delta_latitude)

dlat = np.zeros_like(latitudes).astype('f')
dlat[1:-1] = (latitudes[2:] - latitudes[:-2])/2
dlat[0] = (latitudes[1]-latitudes[0]) 
dlat[-1] = (latitudes[-1]-latitudes[-2])


dy = 111.0*dlat*1000
f = 2.0*omeg_e*np.sin(latitudes*(np.pi/180))


beta = np.zeros_like(latitudes).astype('f')
beta[1:-1] = (f[2:] - f[:-2])/(2*dy[1:-1])
beta[0] = (f[1] - f[0])/(dy[1]) 
beta[-1] = (f[-1] - f[-2])/(dy[-1]) 


if beta_plane == 'True':
    [lind] = np.where(latitudes==lat_info)
    tmp = beta
    tmp[:] = beta[lind]
    del beta
    beta = tmp
    del tmp

ubar = np.zeros_like(latitudes).astype('f')
latr = latitudes*np.pi/180.
latd_shift = 0.0
latr_shift = latd_shift*np.pi/180.
if jet_option == 1:
    ubar[:] = 0.0
elif jet_option == 2:
    ubar[:] = umax
elif jet_option == 3:
    ubar = -umax*np.cos(jet_thinning_factor*4.0*latr + latr_shift) + (umax/4.0)
    
    
ind = np.where(latitudes==lat_info)
pole_to_pole = 180.0*111000.0
circums_earth = 360.*111000.*np.cos(latitudes*(np.pi/180.))

k = num_waves_zonal*((2.*np.pi)/circums_earth)
l = num_waves_merid*((2.*np.pi)/pole_to_pole)

if constant_wavenumber == 'True':
    kval = k[ind]
    k[:] = kval

disp = ubar*k - ((k*beta)/(k**2.0 + l**2.0))    
c = disp/k
cg = ubar - (beta*(l**2.0 - k**2.0))/((k**2.0 + l**2.0)**2.0) 

cg_over_c = np.abs(cg/c)
ubar_minus_c = ubar - c

###########################
# Get ready to plot
###########################
xticks = [-90,-60,-30,0,30,60,90]
yticks = np.arange(-60,60.5,10)

latline = np.zeros(np.size(beta))
latline[:] = lat_info

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2)     

###########################
# Plots: Figure 1
###########################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,ubar,'r',linewidth=3.0,label=r'$\overline{u}$')
p2, = ax1.plot(latitudes,ubar_minus_c,'g',linewidth=3.0,label=r'$\overline{u} - c_x$')

ax1.grid(linestyle='-')

ax1.set_xlim([-90.,90.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

ax1.set_ylim(yticks[0],yticks[-1])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)

ax1.set_ylabel(r'Jet speed (m s$^{-1}$)',fontsize=label_fontsize)
ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + "ubar_latitude_RossbyWaves.png"
plt.savefig(save_name, bbox_inches='tight')

###########################
# Plots: Figure 2
###########################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,c,'r',linewidth=3.0,label=r'$c_x$')
p2, = ax1.plot(latitudes,cg,'b',linewidth=3.0,label=r'$c_{gx}$')



ax1.grid(linestyle='-')

ax1.set_xlim([-90.,90.])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

ax1.set_ylim(yticks[0],yticks[-1])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)

ax1.set_ylabel(r'Phase speed (m s$^{-1}$)',fontsize=label_fontsize)
ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + "phase_speed_group_velocity_latitude_Rossbywaves.png"
plt.savefig(save_name, bbox_inches='tight')

###########################
# Plots: Figure 3
###########################
if 1 == 1:
    fig = pylab.figure(**figprops)   # New figure  
    ax1 = fig.add_subplot(1, 1, 1)

    p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
    ax1.axvline(lat_info,linewidth=4, color='k')
    p1, = ax1.plot(latitudes,k,'r',linewidth=3.0,label='Zonal wavenumber')
    ax1.grid(linestyle='-')

    ax1.set_xlim([-90.,90.])
    ax1.xaxis.set_ticks(xticks)
    ax1.set_xticklabels(xticks)
    ax1.set_ylim([0.,5e-5])
    ax1.set_ylabel(r'$k$',fontsize=label_fontsize)
    ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)
    save_name = imagedir + "zonal_wavenumber_latitude_RossbyWaves.png"
    plt.savefig(save_name, bbox_inches='tight')
    plt.show()
    