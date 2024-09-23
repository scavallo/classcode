# moist_ascent_latent_heating.py
# 
# Computes moist adiabatic parcel ascent and 
# latent heating rates
# 
# Steven Cavallo
# University of Oklahoma
# September 2024
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

import weather_modules as wm


###########################
# User options
###########################
pres1 = 70000. # Pa
pres2 = 30000. # Pa
omega = -1.0 # Pa s-1
reference_latitude = 45.
months2compare = [1,7] # [1,7] will compare January and July, [2,8] will compare February and August, etc.
monthlabels = ['January','July']

label_fontsize = 16
datadir = ''
imagedir = ''
imagename_suff = 'a_vs_b'

infile1 = 'nnrp_t700_vs_lat_monthly.dat'
infile2 = 'nnrp_t700_vs_lat_monthly_2024.dat'

###########################
# constants
###########################
cp = 1004.  # J K-1 kg-1
cv = 717.    # J K-1 kg-1
Rd = 287.    # J K-1 kg-1
Rv = 461.    # J K-1 kg-1
g  = 9.81    # m s-2
p0 = 100000. # Pa
epsi = Rd/Rv
Lv = 2.5*10**(6.0) # J kg-1
p0 = 100000. # Pa

###########################
# Be very careful if editing below
###########################
gammad = g/cp # Dry adiabatic lapse rate is a constant
deltap = pres2-pres1
exner = (((pres1+pres2)/2.0)/p0)**(Rd/cp)

indat = np.loadtxt(datadir+infile1, skiprows=1)       
lats = indat[:,0]
T1a = indat[:,months2compare[0]]
T1b = indat[:,months2compare[1]]
[refind] = np.where(lats==reference_latitude)
    
del indat
indat = np.loadtxt(datadir+infile2, skiprows=1)       
lats2 = indat[:,0]
T2a = indat[:,months2compare[0]]
T2b = indat[:,months2compare[1]]
[refind2] = np.where(lats==reference_latitude)  

for ii in range(1,5):
    if ii == 1:
        Tnow = T1a
    if ii == 2:
        Tnow = T1b
    if ii == 3:
        Tnow = T2a
    if ii == 4:
        Tnow = T2b
    esnow = wm.claus_clap(Tnow) # saturation vapor pressure
    wsnow = (epsi*esnow)/(pres1-esnow) # saturation mixing ratio
    rhonow = pres1/(Rv*Tnow) # Density
    
    
    # Compute moist adiabatic lapse rate
    fnow = (( (1 + (Lv*wsnow)/(Rd*Tnow))) / (1.0 + (wsnow*(Lv**2)/(cp*Rv*Tnow**2.0)) ) )
    sat_lapse_now = gammad*fnow
    
    # Vertical velocity, hydrostatic
    wnow = -omega/(rhonow*g)
    print('w in profile %d is %8.5f m s-1' %(ii,wnow[refind]))   
    
    deltaz_now = (-1.0*deltap)/(rhonow*g)
    print('Cloud depth is %7.2f m in profile %d' % (deltaz_now[refind],ii))
    
    # Temperature at top of ascent layer (T2) is equal to what it was initially (T1) plus slope (-gamma_m) * deltaz
    T2now = Tnow + (-1.0*sat_lapse_now*deltaz_now)
    rho2now=pres1/(Rv*T2now) 
    es2now = wm.claus_clap(T2now)
    ws2now = (epsi*es2now)/(pres2-es2now)   
    qs1now = wsnow/(1+wsnow)
    qs2now = ws2now/(1+ws2now)
    
    print('Lower- and upper-level temperatures in profile %d are %7.2f K and %7.2f' %(ii,Tnow[refind],T2now[refind]))

    # Convert to K/km 
    sat_lapse_now = sat_lapse_now*1000. # per km

    # Fractional departure from dry adiabatic
    fd_now = 1.0 - fnow

    deltaqs = (qs2now - qs1now)
    deltat = deltap/omega # Time it takes parcel to travel deltap (through cloud)

    # You will need to fill in the 2 lines below
    dthetadt_now = 
    dthetadt_perday_now = 
    
    if ii == 1:
        temp_profile1 = Tnow
        sat_lapse_profile1 = sat_lapse_now
        fraction_profile1 = fnow
        fractional_departure_profile1 = fd_now
        dthetadt_perday_profile1 = dthetadt_perday_now
    if ii == 2:
        temp_profile2 = Tnow
        sat_lapse_profile2 = sat_lapse_now
        fraction_profile2 = fnow
        fractional_departure_profile2 = fd_now
        dthetadt_perday_profile2 = dthetadt_perday_now
    if ii == 3:
        temp_profile3 = Tnow
        sat_lapse_profile3 = sat_lapse_now
        fraction_profile3 = fnow
        fractional_departure_profile3 = fd_now
        dthetadt_perday_profile3 = dthetadt_perday_now
    if ii == 4:
        temp_profile4 = Tnow
        sat_lapse_profile4 = sat_lapse_now
        fraction_profile4 = fnow
        fractional_departure_profile4 = fd_now
        dthetadt_perday_profile4 = dthetadt_perday_now
    
###########################
# Setup figures
###########################
y0 = np.zeros_like(temp_profile1).astype('f') #
x = np.linspace(0,len(temp_profile1),np.size(temp_profile1))

xplot = lats
xlabel = r'Latitude ($^{\circ}$N)'
xlims = [np.min(lats), np.max(lats)]

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2)

nticks = 9
yticks_ax1 = np.linspace(2,10,nticks)
yticks_ax2 = np.linspace(0,0.8,nticks)  
###########################
# Figure 1
###########################
fig = pylab.figure(**figprops)  
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1a, = ax1.plot(xplot,sat_lapse_profile1,'b',linewidth=3.0,label=r'$\Gamma_m$ (' + monthlabels[0] + ')')
p1b, = ax1.plot(xplot,sat_lapse_profile3,'b--',linewidth=3.0,label=r'$\Gamma_m$ (' + monthlabels[0] + ' 2024)')
p1c, = ax1.plot(xplot,sat_lapse_profile2,'r',linewidth=3.0,label=r'$\Gamma_m$ (' + monthlabels[1] + ')')
p1c, = ax1.plot(xplot,sat_lapse_profile4,'r--',linewidth=3.0,label=r'$\Gamma_m$ (' + monthlabels[1] + '2024)')


legend = ax1.legend(loc='upper left', shadow=True)
ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])

ax1.set_ylim([yticks_ax1[0],yticks_ax1[1]])
ax1.yaxis.set_ticks(yticks_ax1)
ax1.set_yticklabels(yticks_ax1)
ax1.set_ylabel('Moist adiabatic lapse rate',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + 'moist_adiabatic_' + imagename_suff + '.png'
plt.savefig(save_name, bbox_inches='tight')

###########################
# Figure 2
###########################
fig = pylab.figure(**figprops)   
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(xplot,fraction_profile1,'b',linewidth=3.0,label='f (' + monthlabels[0] + ')')
p2, = ax1.plot(xplot,fraction_profile3,'b--',linewidth=3.0,label='f (' + monthlabels[0] + ' 2024)')
p3, = ax1.plot(xplot,fraction_profile2,'r',linewidth=3.0,label='f (' + monthlabels[1] + ')' )
p4, = ax1.plot(xplot,fraction_profile4,'r--',linewidth=3.0,label='f (' + monthlabels[1] + ' 2024)' )
legend = ax1.legend(loc='upper left', shadow=True)
ax1.axvline(xplot[refind],linewidth=4, color='k')
ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])
ax1.set_ylabel(r'Fraction of $\Gamma_d$',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + 'moist_adiabatic_factor_' + imagename_suff + '.png'
plt.savefig(save_name, bbox_inches='tight')

###########################
# Figure 3
###########################
fig = pylab.figure(**figprops)   
ax1 = fig.add_subplot(1, 1, 1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1, = ax1.plot(xplot,dthetadt_perday_profile1,'b',linewidth=3.0,label=r'$\theta$ LH rate (' + monthlabels[0] + ')')
p2, = ax1.plot(xplot,dthetadt_perday_profile3,'b--',linewidth=3.0,label=r'$\theta$ LH rate (' + monthlabels[0] + ' 2024)')
p3, = ax1.plot(xplot,dthetadt_perday_profile2,'r',linewidth=3.0,label=r'$\theta$ LH rate (' + monthlabels[1] + ')')
p4, = ax1.plot(xplot,dthetadt_perday_profile4,'r--',linewidth=3.0,label=r'$\theta$ LH rate (' + monthlabels[1] + ' 2024)')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])
ax1.set_ylabel('Latent heating rate (K day-1)',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + 'latent_heating_rates_' + imagename_suff +'.png'
plt.savefig(save_name, bbox_inches='tight')
plt.show()

