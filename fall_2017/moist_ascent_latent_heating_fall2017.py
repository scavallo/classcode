# moist_ascent_latent_heating.py
# 
#
# Computes moist adiabatic parcel ascent and 
# latent heating rates.
# 
# This WILL NOT work out of the box for anybody.  It is designed this way for the homework.
# Use this as code skeleton to guide you through the calculations that you need to make.  
# Pay attention to the comments for further direction.
# 
# Steven Cavallo
# University of Oklahoma
# September 2017
###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime

# I am importing weather modules and calling it `wm'
# So if I want to use a function in weather modules later, I would refer to it by `wm'.  
# For example:
#   new_value = wm.function_name_in_weather_modules(arguments_needed_for_function)
import weather_modules as wm

###########################
# User options
###########################
pres1 = 70000. # Base of cloud in Pa
pres2 = -99999. # Top of cloud in Pa.  You need to enter this in based on the problem
omega = -1.00 # Pa s-1
reference_latitude = 45. # Not critical, just used for printing values
read_temperature_fromfile = 'True' # If you are reading in temperature data from a file, set to 'True'

label_fontsize = 16
datadir = ''  # Directory where your data are located
imagedir = '' # Directory where you want to save your images
infile1 = 'nnrp_t700_vs_lat_monthly.dat' # Name of the input data file
###########################
# constants
###########################
cp = 1004.   # J K-1 kg-1
cv = 717.    # J K-1 kg-1
Rd = 287.    # J K-1 kg-1
Rv = 461.    # J K-1 kg-1
g  = 9.81    # m s-2
p0 = 100000. # Pa
epsi = Rd/Rv
Lv = 2.5*10**(6.0) # J kg-1
p0 = 100000. # Pa
###########################
# Calculations 
###########################
gammad = g/cp # Dry adiabatic lapse rate is a constant

if read_temperature_fromfile == 'True':
    ab = np.loadtxt(datadir+infile1, skiprows=1)       
    lats = ab[:,0]
    T1a = ab[:,1]
    # You need to read in the July temperature profile here
    
    
    [refind] = np.where(lats==reference_latitude)
else:
    T1 = np.arange(-40,41,0.25)+273.15 # Virtual temperatures
    T1a = T1
    [refind] = np.where(T1a == 273.15)
    
# weather_modules.py has a function to compute saturation vapor pressure.  
# If weather_modules.py does not work for you, then use that function as a 
# guide for placing the calculation below.  
es1 = wm.claus_clap(T1a) # saturation vapor pressure.  
ws1a = (epsi*es1)/(pres1-es1) # saturation mixing ratio
es1 = wm.claus_clap(T1b) # saturation vapor pressure
ws1b = (epsi*es1)/(pres1-es1) # saturation mixing ratio

# Compute moist adiabatic lapse rate, first in z-coords
# You need to compute fa and fb
sat_lapse_a = gammad*fa
sat_lapse_b = gammad*fb

Rho_a=pres1/(Rv*T1a) # Density, moist air
Rho_b=pres1/(Rv*T1b) 

deltap = pres2-pres1
deltaz_a = (-1.0*deltap)/(Rho_a*g)
deltaz_b = (-1.0*deltap)/(Rho_b*g)
print('Cloud depth is %7.2f m in profile A and %7.2f in profile B' % (deltaz_a[refind],deltaz_b[refind]))


# Temperature at top of ascent layer (T2) is equal to what it was initially (T1) plus slope (-gamma_m) * deltaz
T2a = T1a + (-1.0*sat_lapse_a*deltaz_a)
# You still need to do this again for T1b
print('Lower- and upper-level temperatures in profile A are %7.2f K and %7.2f' %(T1a[refind],T2a[refind]))
print('Lower- and upper-level temperatures in profile B are %7.2f K and %7.2f' %(T1b[refind],T2b[refind]))


# Convert to K/km 
sat_lapse_a = sat_lapse_a*1000. # per km
sat_lapse_b = sat_lapse_b*1000. # per km

# Fractional departure from dry adiabatic
# You need to compute this

# Exner function
exner = (((pres1+pres2)/2.0)/p0)**(Rd/cp)
###################
# Profile A
###################
es2 = wm.claus_clap(T2a)
ws2 = (epsi*es2)/(pres2-es2)

qs1 = ws1a/(1+ws1a)
qs2 = ws2/(1+ws2)

# You need to calculate deltaqs and deltat here

dthetadt_a = -1.0*(Lv/(exner*cp))*(deltaqs/deltat)
dthetadt_perday_a = dthetadt_a * 86400.

###################
# Profile B
###################
# You need to repeat calculations for Profile B

###################
# Get things ready for the plots
###################
y0 = np.zeros_like(T1a).astype('f') #
x = np.linspace(0,len(T1a),np.size(T1a))
if read_temperature_fromfile == 'True':
    xplot = lats
    xlabel = r'Latitude ($^{\circ}$N)'
    xlims = [np.min(lats), np.max(lats)]
else:
    xplot = T1C
    xlims = [-40.,40.]
    xlabel = r'Temperature ($^{\circ}$C)'

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2)

###################
# Plot
###################

nticks = 9
yticks_ax1 = np.linspace(2,10,nticks)
yticks_ax2 = np.linspace(0,0.8,nticks)  

fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
ax2 = ax1.twinx()
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(label_fontsize)
    tick.label1.set_fontweight('bold')
p1a, = ax1.plot(xplot,sat_lapse_a,'r',linewidth=3.0,label=r'$\Gamma_m$ (January)')
p1a, = ax1.plot(xplot,sat_lapse_b,'r--',linewidth=3.0,label=r'$\Gamma_m$ (July)')
p2a, = ax2.plot(xplot,fd_a,'g',linewidth=3.0,label='1-f (January)')
p2b, = ax2.plot(xplot,fd_b,'g--',linewidth=3.0,label='1-f (July)')

legend = ax1.legend(loc='lower left', shadow=True)
legend = ax2.legend(loc='upper left', shadow=True)

ax1.grid(True, linestyle='-')
ax1.set_xlim([xlims[0],xlims[1]])

ax1.set_ylim([yticks_ax1[0],yticks_ax1[1]])
ax2.set_ylim([yticks_ax2[0],yticks_ax2[1]])
ax1.yaxis.set_ticks(yticks_ax1)
ax1.set_yticklabels(yticks_ax1)
ax2.yaxis.set_ticks(yticks_ax2)
ax2.set_yticklabels(yticks_ax2)

ax1.set_ylabel('Moist adiabatic lapse rate',fontsize=label_fontsize)
ax2.set_ylabel('Fractional departure from dry adiabatic',fontsize=label_fontsize)
ax1.set_xlabel(xlabel,fontsize=label_fontsize)
save_name = imagedir + "moist_adiabatic_hw3_2017.png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()


