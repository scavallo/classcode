
# coding: utf-8

# In[ ]:

# saturation_mixing_ratio
#
# Computes saturation mixing, w_s, given a temperature and pressure.
#
# Steven Cavallo
# University of Oklahoma
# September 2015


# In[5]:

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
# New import: Requires weather_modules.py and mstats.py.  
# The below says to import all functions inside weather_modules.py as wm.  So you must use "wm" to access them.
import weather_modules as wm


# In[6]:

T = 0 # degrees Celsius
p = 900 # hPa
Rd = 287. # The dot makes it a floating point value
Rv = 461.


# In[7]:

T = T + 273.15 # Convert to Kelvin
p = p * 100 # Convert to Pa
# Saturation mixing ratio is a function of saturation vapor pressure and total pressure.  Given temperature, we know we can compute saturation vapor pressure:
es = wm.claus_clap(T) # Pa
print 'Saturation vapor pressure is ', es/100, ' hPa' 
# Now we can compute saturation mixing ratio:
epsi = Rd/Rv
ws = epsi*(es/(p-es))*1000 # Convert to g/kg
print 'The saturation mixing ratio at ', T-273.15, ' degrees Celsius and ', p/100, ' hPa is ', ws, ' g/kg'


# In[8]:

# Just for fun, let's check to see what the error in saturation mixing ratio is by the approximation we made to the formula.
# ws = 0.622*(es/p).  
ws_approx = epsi*(es/p)*1000
ws_err = ws - ws_approx
perc_acc = ws_approx/ws*100
Tc = T-273.15
print "The accuracy in making the approximation at %5.2f degrees C and %5.2f hPa is %5.2f percent" %(Tc, p/100, perc_acc)
print 'The error value is ', ws_err, ' g/kg'


# In[ ]:




# In[ ]:



