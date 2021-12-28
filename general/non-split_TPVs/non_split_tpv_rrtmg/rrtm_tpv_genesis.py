##############################################################################
#Imports
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
import climlab
import xarray as xr
import scipy.integrate as sp 
import scipy.ndimage
from climlab.utils.thermo import pseudoadiabat

##############################################################################
#Options and data paths 
##############################################################################
dataDir = '/Users/matthewbray/Desktop/Classes/Undergrad/Senior_Spring/Stats/project/data/ERA/'
saveLoc = '/Users/matthewbray/Desktop/Classes/Undergrad/Senior_Spring/Stats/project/output/'
numLev = 100 #number of climlab pressure levels
includeWV = True #include water vapor in RRTMG_LW?

eraProfilePressure_50 = [100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450,
                         500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875,
                         900, 925, 950, 975, 1000]
                        #pressure levels for 0.5 deg ERA-I data (RH)
eraProfilePressure_75 = np.array([200, 300, 500, 700, 850, 925, 1000])
                        #pressure levels for 0.75 deg ERA-I data (T, z, zeta, PV)
##############################################################################
#Define climlab functions to make idealized and real temperature profiles
##############################################################################
def generate_idealized_temp_profile(SST, plevs, Tstrat=200):
    """
    Generates an idealized temperature profile with specified SST and Tstrat
    """
    solution = sp.odeint(pseudoadiabat, SST, np.flip(plevs))
    temp = solution.reshape(-1)
    temp[np.where(temp<Tstrat)] = Tstrat
    return np.flip(temp) # need to re-invert the pressure axis

def make_idealized_column(SST, num_lev=200, Tstrat=200):
    # Set up a column state
    state = climlab.column_state(num_lev=num_lev, num_lat=1)
    # Extract the pressure levels
    plevs = state['Tatm'].domain.axes['lev'].points
    # Set the SST
    state['Ts'][:] = SST
    # Set the atmospheric profile to be our idealized profile
    state['Tatm'][:] = generate_idealized_temp_profile(SST=SST, plevs=plevs, Tstrat=Tstrat)
    return state, plevs

def make_real_column(tempArray,presArray,num_lev=100):
    # Set up a column state
    state = climlab.column_state(num_lev=num_lev, num_lat=1)
    # Extract the pressure levels
    plevs = state['Tatm'].domain.axes['lev'].points
    # Set the SST
    surfInd = np.argmax(presArray)
    state['Ts'][:] = tempArray[surfInd]
    # Set the atmospheric profile
    state['Tatm'][:] = np.interp(plevs,presArray,tempArray)
    return state, plevs

##############################################################################
#Read in ERA-I profiles in npy format taken from AAARG for the non-split cases
#Temp, RH, PV, vortex latitude, geopotential height, and vertical vorticity
#At 0 hour and 24 hour lags from genesis time, at genesis location
##############################################################################
nsGenTemp0 = np.load(dataDir + 'ns_t_profile_lag_0.npy')
nsGenRH0 = np.load(dataDir + 'ns_rh_profile_lag_0.npy')
nsGenPV0 = np.load(dataDir + 'ns_pv_profile_lag_0.npy')
nsGenLat0 =  np.load(dataDir + 'ns_lats_lag_0.npy')
nsGenZ0 = np.load(dataDir + 'ns_z_profile_lag_0.npy')
nsGenZeta0 = np.load(dataDir + 'ns_zeta_profile_lag_0.npy')

nsGenTemp24 = np.load(dataDir + 'ns_t_profile_lag_24.npy')
nsGenRH24 = np.load(dataDir + 'ns_rh_profile_lag_24.npy')
nsGenPV24 = np.load(dataDir + 'ns_pv_profile_lag_24.npy')
nsGenLat24 =  np.load(dataDir + 'ns_lats_lag_24.npy')
nsGenZ24 = np.load(dataDir + 'ns_z_profile_lag_24.npy')
nsGenZeta24 = np.load(dataDir + 'ns_zeta_profile_lag_24.npy')


#Average all non-split cases
nsGenTempMean0 = np.nanmean(nsGenTemp0,0)
nsGenRHMean0 = np.nanmean(nsGenRH0,0)
nsGenPVMean0 = np.nanmean(nsGenPV0,0)
nsGenLatMean0 =  np.nanmean(nsGenLat0)
nsGenZMean0 = np.nanmean(nsGenZ0,0)
nsGenZetaMean0 = np.nanmean(nsGenZeta0,0)

nsGenTempMean24 = np.nanmean(nsGenTemp24,0)
nsGenRHMean24 = np.nanmean(nsGenRH24,0)
nsGenPVMean24 = np.nanmean(nsGenPV24,0)
nsGenLatMean24 =  np.nanmean(nsGenLat24)
nsGenZMean24 = np.nanmean(nsGenZ24,0)
nsGenZetaMean24 = np.nanmean(nsGenZeta24,0)


##############################################################################
#Set up temperature and humidity profiles for climlab from 24 hr lag ERA-I data
#Interpolate everything to climlab regular plevs
##############################################################################
#state, plevs = make_idealized_column(300)
state, plevs = make_real_column(nsGenTempMean24,eraProfilePressure_75,num_lev=numLev)
#h2o = climlab.radiation.water_vapor.ManabeWaterVapor(state=state,relative_humidity=.8)
h2o = climlab.radiation.water_vapor.FixedRelativeHumidity(state=state)
h2o_rh = np.interp(plevs,eraProfilePressure_50,nsGenRHMean24) / 100
tatm = np.array(state['Tatm'])
h2o_e = 611 * np.exp(((2.5*(10**6))/(461))*((1/273)-(1/tatm))) * h2o_rh
h2o_w = (.622 * h2o_e) / (plevs*100 - h2o_e)
h2o_q = h2o_w / (1+h2o_w)
if includeWV == True:
    h2o.q = h2o_q
    h2o.RH_profile = h2o_rh
else: 
    h2o.q = np.zeros_like(h2o_q)
    h2o.RH_profile = np.zeros_like(h2o_rh)


nsGenPV0_interp = np.interp(plevs,eraProfilePressure_75,nsGenPVMean0)
nsGenZ0_interp = np.interp(plevs,eraProfilePressure_75,nsGenZMean0)
nsGenZeta0_interp = np.interp(plevs,eraProfilePressure_75,nsGenZetaMean0)

nsGenPV24_interp = np.interp(plevs,eraProfilePressure_75,nsGenPVMean24)
nsGenZ24_interp = np.interp(plevs,eraProfilePressure_75,nsGenZMean24)
nsGenZeta24_interp = np.interp(plevs,eraProfilePressure_75,nsGenZetaMean24)


##############################################################################
#Set absorber profiles with ozone from Cavallo
##############################################################################
absorber_vmr = {'CO2':412/1e6,
                'CH4':1.65/1e6,
                'N2O':3.06/1e7,
                'O2':0.21,
                'CFC11':0.,
                'CFC12':0.,
                'CFC22':0.,
                'CCL4':0.,
                'O3':0.}

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
ozone_interp=np.interp(plevs,levs_ozone,ozone_profile)
absorber_vmr['O3'] = ozone_interp


##############################################################################
#Run RRTMG_LW
#Calculate necessary terms for PV tendency equation (vertical only)
##############################################################################
rad24 = climlab.radiation.RRTMG_LW(state=state, specific_humidity=h2o.q, return_spectral_olr=False,absorber_vmr = absorber_vmr)
rad24.compute_diagnostics()

f = np.tile(2*(7.292e-05)*np.sin(np.deg2rad(nsGenLatMean24)),plevs.size)
rho = (plevs * 100) / (287 * np.array(h2o.Tatm))

heatRate = rad24.TdotLW
heatRate = scipy.ndimage.gaussian_filter(heatRate, 3)
dheatRatedz = -rho*9.81*np.gradient(heatRate,plevs*100)

dPVdt = ((nsGenZeta24_interp + f) / rho) * dheatRatedz * 10**6
PVDiff = (nsGenPV0_interp - nsGenPV24_interp) * 10 ** 6


##############################################################################
#Plot
##############################################################################
fig, [ax1,ax2,ax3,ax4,ax5] = plt.subplots(1,5,dpi=100,figsize=(12,4),sharey=True)
ax1.plot(h2o.q,h2o.lev,c='b')
ax2.plot(h2o.Tatm.to_xarray(),h2o.lev,c='r')
ax3.plot(heatRate,plevs,c='k')
ax4.plot(dheatRatedz * 1000,plevs,c='k')
ax5.plot(PVDiff,plevs,c='k',label = 'Actual PV Change')
ax5.plot(dPVdt,plevs,c='r',label='Calculated Tendency')
ax5.legend(loc='lower right',fontsize='small')
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax1.set_ylim(1000,50)
ax2.set_ylim(1000,50)
ax3.set_ylim(1000,50)
ax4.set_ylim(1000,50)
ax5.set_ylim(1000,50)
ax3.set_xlim(-2.5,0)
ax1.set_ylabel('Pressure (hPa)')
ax1.set_xlabel('Specific Humidity (g g$^{-1}$)')
ax2.set_xlabel('Temperature (K)')
ax3.set_xlabel('Heating Rate (K day$^{-1}$)')
ax4.set_xlabel('d(Heating Rate)/dz \n (K day$^{-1}$ km$^{-1}$)')
ax5.set_xlabel('Potential Vorticity \nChange (PVU day$^{-1}$)')
plt.tight_layout()
fig.subplots_adjust(top=0.92)
plt.suptitle('Non-likely Split TPV -24 Hour Lag RRTM PV Generation Simulation')
plt.savefig(saveLoc + 'ns_rrtm_PV_generation.png',bbox_inches='tight',dpi=200)






