# read_sounding_data_hw2016.py
# 
# Example plotting script to read in sounding data and plot.
#
# Obtain sounding data from the University of Wyoming website:
# http://weather.uwyo.edu/upperair/sounding.html
# The function readsounding (located within weather_modules.py) will read data 
#    in the format given when type of plot is 'Text: List' from the dropdown menu
# 
# Steven Cavallo
# University of Oklahoma
# October 2016
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
# User settings
###########################
datadir = '/Users/scavallo/data/'
imagedir = '/Users/scavallo/scripts/python_scripts/images/'
infile1 = 'sounding_data_oun_2016101300.dat'
label_fontsize = 16
vertical_coordinate = 'height' # pressure or height
plot_logscale = 'false'
Lv = 2.5*10**6 # J kg-1
cp = 1004 # J K-1 kg-1
g = 9.81

###########################
# Read the sounding data
###########################
data,fieldnames = wm.readsounding(datadir+infile1)
# Print the plotting options
print fieldnames
 
 # Read in and separate out the fields into their own arrays  
pres=data["pres"][1::]
temp=data["temp"][1::]
dwpt=data["dwpt"][1::]
hght=data["hght"][1::]
theta=data["thta"][1::]
thetae=data["thte"][1::]
mixr=data["mixr"][1::]
thetav=data["thtv"][1::]
winddir=data["drct"][1::]
windspeed_kt=data["sknt"][1::]
Spd = windspeed_kt*0.514444

# If you want height above ground level, subtract height at first level...
hght = hght - hght[0]

# convert virtual theta to virtual temperature
tempv = wm.theta_to_temp(thetav,pres*100.)

# Some calcs that are used more than once
delta_z = (hght[1::]-hght[0:-1]) # meters
delta_theta = (theta[1::]-theta[0:-1]) 

# vertical gradient of theta
dtheta_dz_si = delta_theta / delta_z
dtheta_dz = dtheta_dz_si*10.0**3.0 # K km-1


temp_k = temp + 273.15 # Convert to Kelvin
temp_avg = ( tempv[1::] + tempv[0:-1] ) / 2.0;


u = -Spd * np.sin(winddir * (np.pi/180.))
v = -Spd * np.cos(winddir * (np.pi/180.))
delta_u = (u[1::]-u[0:-1]) 
delta_v = (v[1::]-v[0:-1]) 
u_avg = ( u[1::] + u[0:-1] ) / 2.0;
v_avg = ( v[1::] + v[0:-1] ) / 2.0;

Ri = ( ((9.81/temp_avg)*delta_theta)*delta_z ) / (delta_u**2 + delta_v**2);

du_dz_si = (u[1::]-u[0:-1]) / delta_z 
dv_dz_si = (v[1::]-v[0:-1]) / delta_z
du_dz = du_dz_si*10.0**3.0 # m s-1 km-1
dv_dz = dv_dz_si*10.0**3.0 # m s-1 km-1