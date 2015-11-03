
# coding: utf-8

# In[99]:

# mixed_layer_model
# 
#
# Steven Cavallo
# University of Oklahoma
# October 2015


# In[100]:

###########################
# imports
###########################
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os, datetime


import imp
wm = imp.load_source('weather_modules','/Users/scavallo/scripts/python_scripts/weather_modules.py')
from mstats import *


# In[101]:

datadir = '/Users/scavallo/data/'
imagedir = '/Users/scavallo/Documents/Work/Classes/METR5004/fall_2015/homework/images/'
label_fontsize = 16
label_fontsize_subplot = 8
figname = "abl_problem"

A1 = 0.2
A2 = -0.2
max_FHzi = 200.  # w m-2; If changing this, you will need to change coefficient 'a_used' below!
delta_theta_max = 10.0

# Constants
cp = 1004.
rho = 1.0
Ts = 300. 
Rd = 287.
g = 9.81


# In[102]:

t = np.arange(0,86401.0,1)


# In[103]:

T = 86400.
c = 2.0*np.pi/T
FHz0_traditional = max_FHzi*np.sin(c*t)
FHz0_kinematic = FHz0_traditional / (cp*rho)
delta_theta = -1.0*delta_theta_max*np.cos(c*t)+delta_theta_max

ratio = (FHz0_kinematic) / delta_theta

a = 1.0
ratio_guess = -a*(t-(T/2))**3
[ind] = np.where(t==T/2.5)
print ratio[ind], ratio_guess[ind]


#a_actual = ratio[ind]/ratio_guess[ind]
a_used = 1.0e-14

print a_actual, a_used
ratio_guess = -a_used*(t-(T/2))**3.0
delta_theta_calculate = (FHz0_kinematic)/ratio_guess

zi = np.zeros_like(t).astype('f')
c1 = (A1*a_used/4.0)*(-T/2.0)**4.0
print c1
zi = -(A1*a_used/4.0)*(t - (T/2.0))**4.0 + c1
[tinds] = np.where(t>T/2)
zi[tinds] = -(A2*a_used/4.0)*(t[tinds] - (T/2.0))**4.0 + c1
zi_km = zi/1000.

dtheta_dt = ((1+A1)*FHz0_kinematic)/zi
dtheta_dt[tinds] = ((1+A2)*FHz0_kinematic[tinds])/zi[tinds] 
dtheta_dt_perhour = dtheta_dt*60*60
dtheta_dt_perday = dtheta_dt*86400.0

w_e = A1*ratio_guess
w_e[tinds] = A2*ratio_guess[tinds] 
w_e_true = (A1*FHz0_kinematic) / delta_theta
w_e_true[tinds] = (A2*FHz0_kinematic[tinds]) / delta_theta[tinds]


# In[115]:

# deardorff velocity scale
[indnow] = np.where(t==T/4)
wstar = ((g*zi[indnow]/Ts)*FHz0_kinematic[indnow])**(1.0/3.0)
tstar = zi[indnow]/wstar
tstar_minutes = tstar/60.


print indnow, FHz0_kinematic[indnow],FHz0_traditional[indnow], zi[indnow], wstar, tstar_minutes


# In[105]:

thrs = t/(60.0*60.0)
xticks = np.arange(0,25,6)
zeroarr = np.zeros_like(thrs).astype('f')     


# In[106]:

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    # Figure properties for single and stacked plots
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 


# In[107]:

fig = plt.figure()   
ax1 = fig.add_subplot(1, 1, 1)

p1, = ax1.plot(thrs,delta_theta,'b',linewidth=3.0,label=r'$\Delta\theta$')
p2, = ax1.plot(thrs,delta_theta_calculate,'r',linewidth=3.0,label=r'$\Delta\theta$ calculated')
ax1.grid(True,linestyle='-')
ax1.set_ylabel('Cap strength (K)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax1.set_ylim([-20.0,20.0])
ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)


# In[108]:

fig = plt.figure()   
ax1 = fig.add_subplot(1, 1, 1)

p1, = ax1.plot(thrs,ratio,'r',linewidth=3.0,label='rh')
p2, = ax1.plot(thrs,ratio_guess,'r',linewidth=3.0,label='rh')
ax1.grid(True,linestyle='-')
ax1.set_ylabel(r'FHS/$\Delta\theta$',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax1.set_ylim([-1.0,1.0])
ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)

save_name = imagedir + "hw6_ratio.png"
plt.savefig(save_name, bbox_inches='tight')


# In[109]:

fig = plt.figure()   
ax1 = fig.add_subplot(1, 1, 1)


p1, = ax1.plot(thrs,zi_km,'k',linewidth=3.0,label='rh')
ax1.grid(True,linestyle='-')
ax1.set_ylabel(r'$z_i$ (km)',fontsize=label_fontsize_subplot,fontweight='bold')
#ax4.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax1.set_xlim([thrs[0],thrs[-1]])
ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticks)
ax1.set_xlabel('Time since sunrise (hrs)',fontsize=label_fontsize,fontweight='bold')
save_name = imagedir + "hw6_pblheight.png"
plt.savefig(save_name, bbox_inches='tight')


# In[110]:

fig = plt.figure()   
fig.subplots_adjust(hspace=.5)
ax2 = fig.add_subplot(5, 1, 1)
#ax1 = fig.add_subplot(2,1,2)
p1, = ax2.plot(thrs,FHz0_traditional,'r',linewidth=3.0)
p0, = ax2.plot(thrs,zeroarr,'k',linewidth=2.0)
ax2.grid(True,linestyle='-')
#ax2.set_ylabel(r'Surface heat flux (W m$^{-2}$)',fontsize=label_fontsize,fontweight='bold')
ax2.set_ylabel(r'$Q_{HS}$ (W m$^{-2}$)',fontsize=label_fontsize_subplot,fontweight='bold')
#ax2.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax2.set_ylim([-max_FHzi,max_FHzi])
cint = 100
yticks = np.arange(-max_FHzi,max_FHzi+(cint/2),cint)
#yticks = [-200,-100,0,100,200]
ax2.yaxis.set_ticks(yticks)
ax2.set_yticklabels(yticks)
ax2.set_xlim([thrs[0],thrs[-1]])
ax2.xaxis.set_ticks(xticks)
save_name = imagedir + "hw6_surfaceheatflux.png"
#plt.savefig(save_name, bbox_inches='tight')





# In[111]:

#fig = plt.figure()   
ax3 = fig.add_subplot(5, 1, 2)

p1, = ax3.plot(thrs,ratio_guess,'r',linewidth=3.0)
p0, = ax3.plot(thrs,zeroarr,'k',linewidth=2.0)
ax3.grid(True,linestyle='-')
ax3.set_ylabel(r'FHS/$\Delta\theta$ (m s$^{-1}$)',fontsize=label_fontsize_subplot,fontweight='bold')
#ax3.set_xlabel('Time since sunrise (hrs)',fontsize=label_fontsize,fontweight='bold')
ax3.set_xlim([thrs[0],thrs[-1]])
ax3.xaxis.set_ticks(xticks)
#ax3.set_xticklabels(xticks)

save_name = imagedir + "hw6_ratio.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[112]:

#fig = plt.figure()   
ax4 = fig.add_subplot(5, 1, 3)

p1, = ax4.plot(thrs,zi_km,'k',linewidth=3.0,label='rh')
ax4.grid(True,linestyle='-')
ax4.set_ylabel(r'$z_i$ (km)',fontsize=label_fontsize_subplot,fontweight='bold')
#ax4.set_xlabel('Time since sunrise (hours)',fontsize=label_fontsize,fontweight='bold')
ax4.set_xlim([thrs[0],thrs[-1]])
ax4.xaxis.set_ticks(xticks)
#ax4.set_ylim([0,np.max(zi_km)])
#cint = 0.5
#yticks = np.arange(0,np.max(zi_km)+(cint/2),cint)
ax4.set_ylim([0,4])
yticks = np.arange(0,4.1,1.0)
ax4.yaxis.set_ticks(yticks)
ax4.set_yticklabels(yticks)
#ax4.set_xticklabels(xticks)

save_name = imagedir + "hw6_pblheight.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[113]:

#fig = plt.figure()   
ax5 = fig.add_subplot(5, 1, 4)

#p1, = ax1.plot(thrs,dtheta_dt_perday,'r',linewidth=3.0,label='rh')
p2, = ax5.plot(thrs,dtheta_dt_perhour,'b',linewidth=3.0)
p0, = ax5.plot(thrs,zeroarr,'k',linewidth=2.0)             
ax5.grid(True,linestyle='-')
ax5.set_ylabel(r'$\frac{\partial{\theta}}{\partial{t}}$ (K hr$^{-1}$)',fontsize=label_fontsize_subplot,fontweight='bold')
#ax5.set_xlabel('Time since sunrise (hrs)',fontsize=label_fontsize,fontweight='bold')
ax5.set_xlim([thrs[0],thrs[-1]])
ax5.xaxis.set_ticks(xticks)

ax5.set_ylim([-0.6,0.6])
yticks = np.arange(-0.6,0.65,0.3)
ax5.yaxis.set_ticks(yticks)
ax5.set_yticklabels(yticks)

#ax5.set_xticklabels(xticks)

save_name = imagedir + "hw6_dthetadt.png"
#plt.savefig(save_name, bbox_inches='tight')


# In[114]:

#fig = plt.figure()   
ax6 = fig.add_subplot(5, 1, 5)

#p1, = ax1.plot(thrs,w_e_true,'r',linewidth=3.0,label='rh')
p2, = ax6.plot(thrs,w_e,'b',linewidth=3.0)
p0, = ax6.plot(thrs,zeroarr,'k',linewidth=2.0)
ax6.grid(True,linestyle='-')
ax6.set_ylabel(r'$w_e$ (m s$^{-1}$)',fontsize=label_fontsize_subplot,fontweight='bold')
ax6.set_xlabel('Time since sunrise (hrs)',fontsize=label_fontsize,fontweight='bold')
ax6.set_ylim([-0.5,0.5])
yticks = np.arange(-0.5,0.55,0.25)
ax6.yaxis.set_ticks(yticks)
ax6.set_yticklabels(yticks)
ax6.set_xlim([thrs[0],thrs[-1]])
ax6.xaxis.set_ticks(xticks)
ax6.set_xticklabels(xticks)

save_name = imagedir + "hw6_subplots.png"
plt.savefig(save_name, bbox_inches='tight')
plt.show()

