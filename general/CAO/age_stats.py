#########################################################
# age_stats.py by Dylan Lusk                            #
# 07/09/15                                              #
#                                                       #
# Purpose: Get stats on TPV age vs. var from wrf files  #
#                                                       #
#########################################################

import math
import numpy as np
import netCDF4 as nc
import scipy as sp
import scipy.stats
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pylab as pl
from scipy.interpolate import griddata
from scipy.signal import savgol_filter

######Data Load Location##################################
my_20c_dir = '/data01/Research/track_files/20c/new_tracks/'
my_20c_file = my_20c_dir + '2PVUall2day20cfiltered.dat'
my_21c_dir = '/data01/Research/track_files/21c/new_tracks/'
my_21c_file = my_21c_dir + '2PVUall2day21cfiltered.dat'

error = True #plot with error bars
##########################################################

def read_file(my_file):

    myfile = np.loadtxt(my_file)
    vorttemp = myfile[:,3]
    lenTrack = []
    lenPos = []
    lenTrack.append(int(vorttemp[0]))
    lenPos.append(0)
    #find the length of each track
    for i,val in enumerate(vorttemp):
        if i == np.sum(lenTrack) + len(lenTrack)*2:
            lenTrack.append(int(val))
            lenPos.append(i)

    return myfile,lenTrack,lenPos

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t.ppf((1+confidence)/2., n-1)
    return h, m-h, m+h
    
def bootstrap(data,boot_num):
    data = 1.0*np.array(data)
    n = len(data)
    boot_means = np.empty([boot_num])
    for i in np.arange(0,boot_num,1):
        current_boot = np.random.choice(data,n)
        boot_means[i] = np.mean(current_boot)
    return boot_means

def lat_v_age():
    
    my20cfile,len20cTrack,len20cPos = read_file(my_20c_file)
    my21cfile,len21cTrack,len21cPos = read_file(my_21c_file)
    
    tpv_20c_lats = my20cfile[:,0]
    tpv_21c_lats = my21cfile[:,0]
    
    lat_20c_array = np.zeros([21])
    lat_21c_array = np.zeros([21])
    
    lat_20c_var = []
    lat_21c_var = []
    
    counter_20c = 0
    for x,num in enumerate(len20cPos):
        if len20cTrack[x] < 21:
            continue
        counter_20c += 1
        lat_20c_var_row = []
        for y in np.arange(0,len(lat_20c_array),1):
            lat_20c_array[y] += tpv_20c_lats[num+y+1]
            lat_20c_var_row.append(tpv_20c_lats[num+y+1])
        lat_20c_var.append(lat_20c_var_row)
    
    counter_21c = 0
    for x,num in enumerate(len21cPos):
        if len21cTrack[x] < 21:
            continue
        counter_21c += 1
        lat_21c_var_row = []
        for y in np.arange(0,len(lat_21c_array),1):
            lat_21c_array[y] += tpv_21c_lats[num+y+1]
            lat_21c_var_row.append(tpv_21c_lats[num+y+1])
        lat_21c_var.append(lat_21c_var_row)
    
    lat_20c_array = savgol_filter((lat_20c_array/counter_20c),9,2)
    lat_21c_array = savgol_filter((lat_21c_array/counter_21c),9,2)
    
    lat_20c_var = np.array(lat_20c_var)
    lat_21c_var = np.array(lat_21c_var)
    
    if False:
    
        ci_20c = np.zeros([21])
        ci_21c = np.zeros([21])
    
        for i in np.arange(0,len(ci_20c),1):
            ci_20c[i] = mean_confidence_interval(lat_20c_var[:,i])
            ci_21c[i] = mean_confidence_interval(lat_21c_var[:,i])
            
    if True:
        ci_20c_u = np.zeros([21])
        ci_21c_u = np.zeros([21])
        ci_20c_l = np.zeros([21])
        ci_21c_l = np.zeros([21])
        for i in np.arange(0,len(ci_20c_u),1):
            boot_20c = bootstrap(lat_20c_var[:,i],10000)
            boot_21c = bootstrap(lat_21c_var[:,i],10000)
            ci_20c_u[i] = np.percentile(boot_20c - lat_20c_array[i],97.5)
            ci_21c_u[i] = np.percentile(boot_21c - lat_21c_array[i],97.5)
            ci_20c_l[i] = abs(np.percentile(boot_20c - lat_20c_array[i],2.5))
            ci_21c_l[i] = abs(np.percentile(boot_21c - lat_21c_array[i],2.5))
    
    days = np.arange(0,5.25,.25)
    
    if error:
        fig1 = plt.figure(1)
        (p1,caps1,_1) = plt.errorbar(days,lat_20c_array,yerr=[ci_20c_l,ci_20c_u],label='20c',c='black',linewidth=3.0)
        (p2,caps2,_2) = plt.errorbar(days,lat_21c_array,yerr=[ci_21c_l,ci_21c_u],label='21c',ls='--',c='green',linewidth=3.0)
        plt.ylabel('Latitude($^\circ$N)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([66,76])
        plt.legend(loc=4)
        for cap1 in caps1:
            cap1.set_linewidth(3.0)
            cap1.set_markeredgewidth(3.0)
        for cap2 in caps2:
            cap2.set_linewidth(3.0)
            cap2.set_markeredgewidth(3.0)
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/lat_vs_age_bootstrap.png')
        plt.show()
    
    else:
        fig1 = plt.figure(1)
        p1 = plt.plot(days,lat_20c_array,label='20c',c='black',linewidth=3.0)
        p2 = plt.plot(days,lat_21c_array,label='21c',ls='--',c='green',linewidth=3.0)
        plt.ylabel('Latitude($^\circ$N)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([68,76])
        plt.legend(loc=1)
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/lat_vs_age.png')
        plt.show()
        
    
def pres_amp_v_age():
    
    my20cfile,len20cTrack,len20cPos = read_file(my_20c_file)
    my21cfile,len21cTrack,len21cPos = read_file(my_21c_file)
    
    tpv_20c_pres = my20cfile[:,6]
    tpv_21c_pres = my21cfile[:,6]
    
    lat_20c_array = np.zeros([21])
    lat_21c_array = np.zeros([21])
    
    lat_20c_var = []
    lat_21c_var = []
    
    counter_20c = 0
    for x,num in enumerate(len20cPos):
        if len20cTrack[x] < 21:
            continue
        counter_20c += 1
        lat_20c_var_row=[]
        for y in np.arange(0,len(lat_20c_array),1):
            lat_20c_array[y] += tpv_20c_pres[num+y+1]
            lat_20c_var_row.append(tpv_20c_pres[num+y+1])
        lat_20c_var.append(lat_20c_var_row)
    
    counter_21c = 0
    for x,num in enumerate(len21cPos):
        if len21cTrack[x] < 21:
            continue
        counter_21c += 1
        lat_21c_var_row=[]
        for y in np.arange(0,len(lat_21c_array),1):
            lat_21c_array[y] += tpv_21c_pres[num+y+1]
            lat_21c_var_row.append(tpv_21c_pres[num+y+1])
        lat_21c_var.append(lat_21c_var_row)
    
    lat_20c_array = savgol_filter((lat_20c_array/counter_20c),9,2)
    lat_21c_array = savgol_filter((lat_21c_array/counter_21c),9,2)
    
    lat_20c_var = np.array(lat_20c_var)
    lat_21c_var = np.array(lat_21c_var)
    
    ci_20c = np.zeros([21])
    ci_21c = np.zeros([21])
    for i in np.arange(0,len(ci_20c),1):
        ci_20c[i] = mean_confidence_interval(lat_20c_var[:,i])
        ci_21c[i] = mean_confidence_interval(lat_21c_var[:,i])
    
    days = np.arange(0,5.25,.25)
    if error:
        fig1 = plt.figure(1)
        (p1,caps1,_1) = plt.errorbar(days,lat_20c_array,yerr=ci_20c,label='20c',c='black',linewidth=3.0)
        (p2,caps2,_2) = plt.errorbar(days,lat_21c_array,yerr=ci_21c,label='21c',ls='--',c='black',linewidth=3.0)
        plt.ylabel('Pressure Amplitude (hPa)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([20,70])
        plt.legend()
        for cap1 in caps1:
            cap1.set_linewidth(3.0)
            cap1.set_markeredgewidth(3.0)
        for cap2 in caps2:
            cap2.set_linewidth(3.0)
            cap2.set_markeredgewidth(3.0)
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/pres_amp_vs_age_ci.png')
        plt.show()
        
    else:
        fig1 = plt.figure(1)
        p1 = plt.plot(days,lat_20c_array,label='20c',c='black',linewidth=3.0)
        p2 = plt.plot(days,lat_21c_array,label='21c',ls='--',c='black',linewidth=3.0)
        plt.ylabel('Pressure Amplitude (hPa)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([30,70])
        plt.legend()
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/pres_amp_vs_age.png')
        plt.show()
    
def temp_amp_v_age():

    my20cfile,len20cTrack,len20cPos = read_file(my_20c_file)
    my21cfile,len21cTrack,len21cPos = read_file(my_21c_file)
    
    tpv_20c_temp = my20cfile[:,4]
    tpv_21c_temp = my21cfile[:,4]
    
    lat_20c_array = np.zeros([21])
    lat_21c_array = np.zeros([21])
    
    lat_20c_var = []
    lat_21c_var = []
    
    counter_20c = 0
    for x,num in enumerate(len20cPos):
        if len20cTrack[x] < 21:
            continue
        counter_20c += 1
        lat_20c_var_row=[]
        for y in np.arange(0,len(lat_20c_array),1):
            lat_20c_array[y] += tpv_20c_temp[num+y+1]
            lat_20c_var_row.append(tpv_20c_temp[num+y+1])
        lat_20c_var.append(lat_20c_var_row)
    
    counter_21c = 0
    for x,num in enumerate(len21cPos):
        if len21cTrack[x] < 21:
            continue
        counter_21c += 1
        lat_21c_var_row=[]
        for y in np.arange(0,len(lat_21c_array),1):
            lat_21c_array[y] += tpv_21c_temp[num+y+1]
            lat_21c_var_row.append(tpv_21c_temp[num+y+1])
        lat_21c_var.append(lat_21c_var_row)
        
    lat_20c_array = savgol_filter((lat_20c_array/counter_20c),9,2)
    lat_21c_array = savgol_filter((lat_21c_array/counter_21c),9,2)
    
    lat_20c_var = np.array(lat_20c_var)
    lat_21c_var = np.array(lat_21c_var)
    
    ci_20c = np.zeros([21])
    ci_21c = np.zeros([21])
    for i in np.arange(0,len(ci_20c),1):
        ci_20c[i] = mean_confidence_interval(lat_20c_var[:,i])
        ci_21c[i] = mean_confidence_interval(lat_21c_var[:,i])
    
    days = np.arange(0,5.25,.25)

    if error:
        fig1 = plt.figure(1)
        (p1,caps1,_1) = plt.errorbar(days,lat_20c_array,yerr=ci_20c,label='20c',c='black',linewidth=3.0)
        (p2,caps2,_2) = plt.errorbar(days,lat_21c_array,yerr=ci_21c,label='21c',ls='--',c='black',linewidth=3.0)
        plt.ylabel('Temperature Amplitude (K)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([4,15])
        plt.legend()
        for cap1 in caps1:
            cap1.set_linewidth(3.0)
            cap1.set_markeredgewidth(3.0)
        for cap2 in caps2:
            cap2.set_linewidth(3.0)
            cap2.set_markeredgewidth(3.0)
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/temp_amp_vs_age_ci.png')
        plt.show()

    else:    
        fig1 = plt.figure(1)
        p1 = plt.plot(days,lat_20c_array,label='20c',c='black',linewidth=2.0)
        p2 = plt.plot(days,lat_21c_array,label='21c',ls='--',c='black',linewidth=2.0)
        plt.ylabel('Temperature Amplitude (K)', fontsize = 12)
        plt.xlabel('Vortex Age (Days)', fontsize = 12)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.xlim([0,5])
        plt.ylim([4,15])
        plt.legend()
        plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/temp_amp_vs_age.png')
        plt.show()
    
def vortex_lifespan():

    my20cfile,len20cTrack,len20cPos = read_file(my_20c_file)
    my21cfile,len21cTrack,len21cPos = read_file(my_21c_file)
    
    vort_20c_life = np.zeros([max(len20cTrack)])
    for x in np.arange(0,max(len20cTrack),1):
        for num in len20cTrack:
            if x == num:
                vort_20c_life[x] += 1
    
    vort_21c_life = np.zeros([max(len21cTrack)])
    for x in np.arange(0,max(len21cTrack),1):
        for num in len21cTrack:
            if x == num:
                vort_21c_life[x] += 1
                            
    #need to convert the vort_20c_life to days here
    vort_20c_days = np.zeros([(len(vort_20c_life)//4)+1])
    for x in np.arange(0,len(vort_20c_life),1):
        vort_20c_days[x//4] = vort_20c_life[x]
        
    vort_21c_days = np.zeros([(len(vort_21c_life)//4)+1])
    for x in np.arange(0,len(vort_21c_life),1):
        vort_21c_days[x//4] = vort_21c_life[x]
    
    #need to add in line and t-test between them - look at ecmwf stuff
    
    
    days_20c = np.arange(0,len(vort_20c_days),1)
    days_21c = np.arange(0,len(vort_21c_days),1)
    fig1 = plt.figure(1)
    p1 = plt.scatter(days_20c,vort_20c_days,label='20c',marker='o',color='blue',s=30.0)
    p2 = plt.scatter(days_21c,vort_21c_days,label='21c',marker='x',color='green',s=30.0)
    
    #Create a truth array for the nonzero values - otherwise log plot fit fails
    mask_20c = np.nonzero(vort_20c_days)
    mask_21c = np.nonzero(vort_21c_days)
    
    days_m_20c = days_20c[mask_20c]
    days_m_21c = days_21c[mask_21c]
    vort_20c_m = vort_20c_days[mask_20c]
    vort_21c_m = vort_21c_days[mask_21c]

    slope_20c, intercept_20c, r_value_20c, p_value_20c, std_err_20c = scipy.stats.linregress(days_m_20c,np.log10(vort_20c_m))
    bfl_20c = np.zeros([len(days_m_20c)])
    bfl_20c = 10**((slope_20c*days_m_20c)+intercept_20c)
    
    slope_21c, intercept_21c, r_value_21c, p_value_21c, std_err_21c = scipy.stats.linregress(days_m_21c,np.log10(vort_21c_m))
    bfl_21c = np.zeros([len(days_m_21c)])
    bfl_21c = 10**((slope_21c*days_m_21c)+intercept_21c)
    
    p3 = plt.plot(days_m_20c, bfl_20c, label='20c Fit, R^2 = ' + '%.2f' % r_value_20c**2 + ', Slope = ' + '%.2f' % slope_20c, linewidth='2.0')
    p4 = plt.plot(days_m_21c, bfl_21c, label='21c Fit, R^2 = ' + '%.2f' % r_value_21c**2 + ', Slope = ' + '%.2f' % slope_21c, linewidth='2.0')
    plt.ylabel('Number', fontsize = 12)
    plt.xlabel('Vortex Lifespan (Days)', fontsize = 12)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.yscale('log')
    plt.xlim([0,30])
    plt.ylim([1,10**3])
    plt.legend()
    plt.savefig('/home/dylanl/Dropbox/Research/Lusk_Cavallo_2016/vortex_lifespan_fits.png')
    plt.show()
        
##############################################################

lat_v_age()
#pres_amp_v_age()
#temp_amp_v_age()
#vortex_lifespan()
