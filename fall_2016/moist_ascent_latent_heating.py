# compute_latent_heating.py
# 
#
# Computes moist adiabatic parcel ascent and 
# latent heating rates
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
import os, datetime

import weather_modules as wm
from mstats import *
###########################
# User options
###########################
pres1 = 70000. # Pa
pres2 = 30000. # Pa
omega = -1.00 # Pa s-1
reference_latitude = 45.

read_temperature_fromfile = 'True'
label_fontsize = 16
datadir = '/Users/scavallo/data/'
imagedir = '' # Enter your path!
infile1 = 'nnrp_t700_vs_lat.dat'

###########################
# constants
###########################
cp = 1004.  # J K-1 kg-1
cv = 717.    # J K-1 kg-1
Rd = 287.    # J K-1 kg-1
Rv = 461.    # J K-1 kg-1
p0 = 100000. # Pa
epsi = Rd/Rv
Lv = 2.5*10**(6.0) # J kg-1

###########################
# Read file, compute stuff
###########################
if read_temperature_fromfile == 'True':
    ab = np.loadtxt(datadir+infile1, skiprows=0)       
    lats = ab[:,0]
    T1 = ab[:,1]
    [refind] = np.where(lats==reference_latitude)
    print refind
else:
    T1 = np.arange(-40,41,0.25)+273.15 # Virtual temperatures
    [refind] = np.where(T1 == 273.15)
# Check to make sure it read in correctly
mstats(T1)    
es1 = wm.claus_clap(T1) # saturation vapor pressure
ws1 = (epsi*es1)/(pres1-es1) # saturation mixing ratio

# Compute moist adiabatic lapse rate, first in z-coords
sat_lapse = (9.81/cp)*(( (1 + (Lv*ws1)/(Rd*T1))) / (1 + (ws1*(Lv**2)/(cp*Rv*T1**2.0)) ) )
Rho=pres1/(Rv*T1) # Density, moist air
sat_lapse_isobaric = -sat_lapse/(Rho*9.81) # Convert to p-coords

# Temperature at top of ascent layer (T2) is equal to what it was initially (T1) plus slope (-gamma_m) * deltap
T2 = T1 - (sat_lapse_isobaric)*(pres2-pres1) 

# Convert to Celsius (needed for plot if not reading from file)
T1C = T1 - 273.15
T2C = T2 - 273.15

# Convert to K/km and K/hundred hPa
sat_lapse = sat_lapse*1000. # per km
sat_lapse_isobaric = sat_lapse_isobaric*100.*100. # per 100 hPa
# Check values
mstats(sat_lapse)
###########################
# Calculate latent heating and plot below
###########################
