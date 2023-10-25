# random_perturbation_timeseries
#
# Demonstration of Reynolds Averaging calculations
#
# Steven Cavallo
# November 2023
#####################################################
# Imports
#####################################################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
from mstats import *

#####################################################
# User settings
#####################################################
mean = 0
std = 1 
label_fontsize = 12
imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
#####################################################
# Calculations
#####################################################
t = np.linspace(0,120,100)


num_samples = len(t)
plotvar = np.random.normal(mean, std, size=num_samples)
plotvar2 = np.random.normal(mean, std, size=num_samples)

prod_pert = plotvar*plotvar
prod_pert_mean = np.nanmean(prod_pert)

prod_pert_covariance = plotvar*plotvar2
prod_pert_mean_covariance = np.nanmean(prod_pert_covariance)


#####################################################
# Figure setup
#####################################################
golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

titletext1 = r'$a^\prime$'
titletext2 = r'$a^\prime a^\prime$'
titletext3 = r'$a^\prime b^\prime$'

legendtext1 = r'$\overline{a^\prime}$'
legendtext2 = r'$\overline{a^\prime a^\prime}$'
legendtext3 = r'$\overline{a^\prime b^\prime}$'

#####################################################
# Figure 1
#####################################################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)


p1, = ax1.plot(t,plotvar,'b',linewidth=4.0)
ax1.grid(linestyle='-')
ax1.axhline(0,linewidth=4, color='k',label=legendtext1)

ax1.set_xlim([t[0],t[-1]])

ax1.set_ylabel('Value',fontsize=label_fontsize)
ax1.set_xlabel('Time (seconds)',fontsize=label_fontsize)
plt.title(titletext1,fontsize=28)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + 'random_perturbation_zeromean.png'
plt.savefig(save_name, bbox_inches='tight')

#####################################################
# Figure 2
#####################################################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)


p1, = ax1.plot(t,prod_pert,'b',linewidth=4.0)
ax1.grid(linestyle='-')
ax1.axhline(prod_pert_mean,linewidth=4, color='k',label=legendtext2)

ax1.set_xlim([t[0],t[-1]])


ax1.set_ylabel('Value',fontsize=label_fontsize)
ax1.set_xlabel('Time (seconds)',fontsize=label_fontsize)
plt.title(titletext2,fontsize=28)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + 'random_perturbation_products_zeromean.png'
plt.savefig(save_name, bbox_inches='tight')

#####################################################
# Figure 3
#####################################################
fig = pylab.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)

p1, = ax1.plot(t,prod_pert_covariance,'b',linewidth=4.0)
ax1.grid(linestyle='-')
ax1.axhline(prod_pert_mean_covariance,linewidth=4, color='k',label=legendtext3)

ax1.set_xlim([t[0],t[-1]])


ax1.set_ylabel('Value',fontsize=label_fontsize)
ax1.set_xlabel('Time (seconds)',fontsize=label_fontsize)
plt.title(titletext3,fontsize=28)

legend = ax1.legend(loc='upper left', shadow=True)
save_name = imagedir + 'random_perturbation_products_covariance_zeromean.png'
plt.savefig(save_name, bbox_inches='tight')

plt.show()

