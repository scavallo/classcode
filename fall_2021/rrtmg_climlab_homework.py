# rrtmg_climlab_homework
#
# Steven Cavallo
# October 2021
#############################
import numpy as np
import matplotlib.pyplot as plt
import pylab as P
import climlab, pylab
import xarray as xr
import scipy.integrate as sp  #Gives access to the ODE integration package
from climlab.utils.thermo import pseudoadiabat

from mstats import *
import utilities_modules as um
#############################
# User options
#############################
moisture_option = 2 # 1 for Manabe and Wetherald JAS 1967, 
					# 2 for constant relative humidity.  Set relh below.
relh = 0.8 # If moisture_option = 1, this sets the fractional relative humidity at the surface
		   # If moisture_option = 2, this sets a constant (fractional) relative humidity value
ylims = [10,1000]
label_fontsize = 16
plot_logscale = False
datadir = '/Users/scavallo/Documents/data/'
#imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2021/homework/images/'

#############################
absorber_vmr = {'CO2':412/1e6,
                'CH4':1.82/1e6,
                'N2O':0.,
                'O2':0.,
                'CFC11':0.,
                'CFC12':0.,
                'CFC22':0.,
                'CCL4':0.,
                'O3':10/1e6}

#############################

def generate_idealized_temp_profile(SST, plevs, Tstrat=200):
    """
    Generates an idealized temperature profile with specified SST and Tstrat
    """
    solution = sp.odeint(pseudoadiabat, SST, np.flip(plevs))
    temp = solution.reshape(-1)
    temp[np.where(temp<Tstrat)] = Tstrat
    return np.flip(temp) # need to re-invert the pressure axis

def make_idealized_column(SST, num_lev=100, Tstrat=200):
    # Set up a column state
    state = climlab.column_state(num_lev=num_lev, num_lat=1)
    # Extract the pressure levels
    plevs = state['Tatm'].domain.axes['lev'].points
    # Set the SST
    state['Ts'][:] = SST
    # Set the atmospheric profile to be our idealized profile
    state['Tatm'][:] = generate_idealized_temp_profile(SST=SST, plevs=plevs, Tstrat=Tstrat)
    return state
    

#############################
# Make the temperature and moisture profiles
state = make_idealized_column(300)

# Water vapor mixing ratio profile following Manabe and Wetherald JAS 1967 
# Fixed surface relative humidity and a specified fractional profile.  Relative_humidity 
# is the specified surface RH qStrat is the minimum specific humidity, ensuring that there is some water vapor in the stratosphere.

if moisture_option == 1:
    h2o_out = climlab.radiation.water_vapor.ManabeWaterVapor(state=state,relative_humidity=relh)
else:
    qStrat=5e-07
    if relh == 0:
        qStrat = 0.
    h2o_out = climlab.radiation.water_vapor.FixedRelativeHumidity(state=state,relative_humidity=relh, qStrat=qStrat)

	

temperature = state['Tatm'].to_xarray()
pressure = state['Tatm'].domain.axes['lev'].points
h2o_q = h2o_out.q

es = climlab.utils.thermo.clausius_clapeyron(temperature)
e = climlab.utils.thermo.vapor_pressure_from_specific_humidity(pressure, h2o_out.q)
rhfrac = e/es

# Ozone profile
ozonepath = "http://thredds.atmos.albany.edu:8080/thredds/dodsC/CLIMLAB/ozone/apeozone_cam3_5_54.nc"
ozone = xr.open_dataset(ozonepath)
#  Dimensions of the ozone file
lat = ozone.lat
lon = ozone.lon
lev = ozone.lev
# Taking annual, zonal, and global averages of the ozone data
O3_zon = ozone.OZONE.mean(dim=("time","lon"))
#  Set the ozone mixing ratio
weight_ozone = np.cos(np.deg2rad(ozone.lat)) / np.cos(np.deg2rad(ozone.lat)).mean(dim='lat')
O3_global = (O3_zon * weight_ozone).mean(dim='lat')

levs_ozone = np.array(lev)
ozone_profile = np.array(O3_global.values)
ozone_interp=np.interp(pressure,levs_ozone,ozone_profile)
absorber_vmr['O3'] = ozone_interp

# RRTMG radiative parameterization.  This is the same parameterization used in both weather and climate models
rad = climlab.radiation.RRTMG_LW(state=state, specific_humidity=h2o_q,icld=1, irng=1, idrv=0, 
								 permuteseed=300, inflglw=2, iceflglw=1, liqflglw=1, tauc=0.0, 
								 tauaer=0.0, return_spectral_olr=False,absorber_vmr = absorber_vmr)
# Get the heating rates out
rad._compute_heating_rates()
lw_up = rad.LW_flux_up
lw_down = rad.LW_flux_down
pressure2 = rad.lev_bounds
lw_heating_rate = rad.TdotLW
lw_net = lw_up - lw_down

# Case for no water vapor
zero_moisture = np.zeros(len(h2o_q))
rad2 = climlab.radiation.RRTMG_LW(state=state, specific_humidity=zero_moisture,icld=1, irng=1, idrv=0, 
								 permuteseed=300, inflglw=2, iceflglw=1, liqflglw=1, tauc=0.0, 
								 tauaer=0.0, return_spectral_olr=False,absorber_vmr = absorber_vmr)
rad2._compute_heating_rates()
lw_heating_rate_nomoisture = rad2.TdotLW

# Case for no ozone
absorber_vmr['O3'] = 0.
rad3 = climlab.radiation.RRTMG_LW(state=state, specific_humidity=h2o_q,icld=1, irng=1, idrv=0, 
								 permuteseed=300, inflglw=2, iceflglw=1, liqflglw=1, tauc=0.0, 
								 tauaer=0.0, return_spectral_olr=False,absorber_vmr = absorber_vmr)
rad3._compute_heating_rates()
lw_heating_rate_noozone = rad3.TdotLW

# Plot the profiles
#############################
# Plot setup
#############################
# These are the levels I want to show on my plot explicitly   
if plot_logscale == False:
    yticks = [1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100.,50.]
else:
    yticks = [1000.,850.,700.,500.,300.,200.,100.,50.,1.]   
x0 = np.zeros_like(pressure).astype('f')

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

#############################
# Figure 1
#############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(ozone_interp,pressure,'r',linewidth=3.0,label='Ozone')
if plot_logscale == True:
    ax1.set_yscale('log')

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)    
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([ylims[1],ylims[0]])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Ozone',fontsize=label_fontsize)
save_name = imagedir + "rrtm_ozone_profile.png"
plt.savefig(save_name, bbox_inches='tight')

#############################
# Figure 2
#############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(temperature,pressure,'r',linewidth=3.0,label='Temperature')
if plot_logscale == True:
    ax1.set_yscale('log')

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)    
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([ylims[1],ylims[0]])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Temperature (K)',fontsize=label_fontsize)
save_name = imagedir + "rrtm_temperature_profile.png"
plt.savefig(save_name, bbox_inches='tight')

#############################
# Figure 3
#############################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(h2o_q,pressure,'g',linewidth=3.0,label='Specific humidity')
if plot_logscale == True:
    ax1.set_yscale('log')

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)    
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([ylims[1],ylims[0]])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel('Specific humidity (g/kg)',fontsize=label_fontsize)
save_name = imagedir + "rrtm_moisture_profile.png"
plt.savefig(save_name, bbox_inches='tight')

#############################
 # Figure 4
#############################   
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(lw_up,pressure2,'b',linewidth=3.0,label=r'Longwave$\uparrow$')
p2, = ax1.plot(lw_down,pressure2,'r',linewidth=3.0,label=r'Longwave$\downarrow$')
p3, = ax1.plot(lw_net,pressure2,'k',linewidth=3.0,label=r'Net Longwave')
if plot_logscale == True:
    ax1.set_yscale('log')

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)    
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([ylims[1],ylims[0]])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='upper right', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel(r'Flux (W m$^2$)',fontsize=label_fontsize)
save_name = imagedir + "rrtm_lw_profiles.png"
plt.savefig(save_name, bbox_inches='tight')

#############################
 # Figure 5
#############################   
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)
p1, = ax1.plot(lw_heating_rate,pressure,'b',linewidth=4.0,label=r'Longwave heating rate')
p2, = ax1.plot(lw_heating_rate_nomoisture,pressure,'g',linewidth=2.0,label=r' No water vapor')
p3, = ax1.plot(lw_heating_rate_noozone,pressure,'m',linewidth=2.0,label=r' No ozone')
if plot_logscale == True:
    ax1.set_yscale('log')

um.label_options(ax1,fontsize=label_fontsize,xaxis_opt=True,yaxis_opt=True,bold_opt=True)    
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_ylim([ylims[1],ylims[0]])
ax1.yaxis.set_ticks(yticks)
ax1.set_yticklabels(yticks)
ax1.grid(True, linestyle='-')
legend = ax1.legend(loc='lower left', shadow=True)
ax1.set_ylabel('Pressure (hPa)',fontsize=label_fontsize)
ax1.set_xlabel(r'Heating rate (K day$^{-1}$)',fontsize=label_fontsize)
save_name = imagedir + "rrtm_heatingrate_profiles.png"
plt.savefig(save_name, bbox_inches='tight')

plt.show()


