# Rossby_waves_hw2017.py
#
# Script to plot Rossby wave phase speed and group velocity, given a background flow, wavenumber, and latitude
#
# Steven Cavallo
# University of Oklahoma
# September 2017


# In[2]:

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
# Default settings are:
#    umax = 20.0
#    jet_thinning_factor = 1.0
#    beta_plane = 'False'
#    constant_wavenumber = 'False' 
#    num_waves_merid = 0.0
#    num_waves_zonal = 5.0
#    lat_info = 45.0
#    delta_latitude = 0.5
#
# Set yours below this point!
umax = 20.0
jet_option = 4 # 1 for no jet/ubar=0, 
               # 2 for constant jet/ubar=constant=umax, 
               # 3 for zonally varying jet with umax amp
               # 4 to determine jet/ubar by RW phase speed 
phasespeed_jet_option4 = 0.0 # Frequency of RWs if jet_option == 4
jet_thinning_factor = 1.0
beta_plane = 'False'
constant_wavenumber = 'True' # 'True' : Wavenumber at every latitude is constant.  This implies the number of waves will vary with latitude
                             # 'False': Wavenumber varies with latitude, but number of waves remains constant
num_waves_merid = 0.0
num_waves_zonal = 5.0
lat_info = 45.0

delta_latitude = 0.5
label_fontsize = 16
image_suffix = 'jetoption3_ubar20'
imagedir = ''
#######################
# END User parameters
#######################

R_earth = 6371200;

x = np.linspace(0,2*np.pi,100)
omeg_e = (2*np.pi) / (24*3600)
latitudes = np.arange(-90,90+(delta_latitude/2),delta_latitude)

dlat = np.zeros_like(latitudes).astype('f')
dlat[1:-1] = (latitudes[2:] - latitudes[:-2])/2
dlat[0] = (latitudes[1]-latitudes[0]) 
dlat[-1] = (latitudes[-1]-latitudes[-2])

dy = 111.0*dlat*1000
f = 2.0*omeg_e*np.sin(latitudes*(np.pi/180))
beta = (2.0*omeg_e*np.cos(latitudes*(np.pi/180))) / R_earth

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
if ( (jet_option == 1) or (jet_option == 4)) :
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

num_waves_zonal_vec = np.zeros_like(latitudes).astype('f')
num_waves_merid_vec = np.zeros_like(latitudes).astype('f')

if constant_wavenumber == 'True':
    k[:] = k[ind]
    l[:] = l[ind]
wavelen = ((2.0*np.pi)/k)/1000.
print('Wavelength is %7.2f km' %(wavelen[ind]))
num_waves_zonal_vec = k/((2.*np.pi)/circums_earth)
num_waves_zonal_merid = l/((2.*np.pi)/circums_earth)


freq = ubar*k - ((k*beta)/(k**2.0 + l**2.0))   
cx = freq/k
if jet_option == 4:
    cx[:] = phasespeed_jet_option4
    freq = cx*k
    ubar = (freq/k) + (beta/(k**2.0+l**2.0))
     

cgx = ubar - (beta*(l**2.0 - k**2.0))/((k**2.0 + l**2.0)**2.0) 

cgx_over_cx = np.abs(cgx/cx)
ubar_minus_cx = ubar - cx


#######################
# Get stuff ready for plotting
#######################

xticks = [-90,-60,-30,0,30,60,90]
yticks = np.arange(-60,60.5,10)

latline = np.zeros(np.size(beta))
latline[:] = lat_info

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 


#######################
# New figure
#######################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,ubar,'r',linewidth=4.0,label=r'$\overline{u}$')
p2, = ax1.plot(latitudes,ubar_minus_cx,'g',linewidth=2.0,label=r'$\overline{u} - c_x$')

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
save_name = imagedir + 'ubar_latitude_RossbyWaves_' + image_suffix + '.png'
plt.savefig(save_name, bbox_inches='tight')

#######################
# New figure
#######################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
ax1.axvline(lat_info,linewidth=4, color='k')

p1, = ax1.plot(latitudes,cx,'r',linewidth=4.0,label=r'$c_x$')
p2, = ax1.plot(latitudes,cgx,'b',linewidth=2.0,label=r'$c_{gx}$')



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
save_name = imagedir + 'phase_speed_group_velocity_latitude_Rossbywaves_' + image_suffix + '.png'
plt.savefig(save_name, bbox_inches='tight')
#plt.show()


#######################
# New figure
#######################
if 1 == 1:
    nticks = 5
    yticks_ax1 = np.linspace(0,8000,nticks)
    yticks_ax2 = np.linspace(0,8,nticks)    
    
    fig = pylab.figure(**figprops)   # New figure  
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twinx()

    p0, = ax1.plot(latitudes,np.zeros(np.size(beta)),'k',linewidth=4.0)
    ax1.axvline(lat_info,linewidth=4, color='k')
    p1, = ax1.plot(latitudes,wavelen,'r',linewidth=3.0,label='Wavelength')
    p2, = ax2.plot(latitudes,num_waves_zonal_vec,'g',linewidth=3.0,label='Number of waves')
    legend = ax1.legend(loc='center left', shadow=True)
    legend = ax2.legend(loc='upper left', shadow=True)
    ax1.grid(linestyle='-')
    ax2.grid(True, linestyle='-')

    ax1.set_xlim([-90.,90.])
    ax1.xaxis.set_ticks(xticks)
    ax1.set_xticklabels(xticks)

    ax1.set_ylim([yticks_ax1[0],yticks_ax1[1]])
    ax2.set_ylim([yticks_ax2[0],yticks_ax2[1]])
    ax1.yaxis.set_ticks(yticks_ax1)
    ax1.set_yticklabels(yticks_ax1)
    ax2.yaxis.set_ticks(yticks_ax2)
    ax2.set_yticklabels(yticks_ax2)
    
    ax1.set_ylabel(r'$\lambda$ (km)',fontsize=label_fontsize)
    ax2.set_ylabel('Number of zonal waves',fontsize=label_fontsize)
    ax1.set_xlabel('Latitude (degrees)',fontsize=label_fontsize)
    save_name = imagedir + 'zonal_wavelength_latitude_RossbyWaves_' + image_suffix + '.png'
    plt.savefig(save_name, bbox_inches='tight')
    plt.show()

