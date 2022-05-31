#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from scipy import ndimage
import scipy.io as sio
from scipy import stats

import warnings
warnings.filterwarnings("ignore")

###################################
# Set user options
###################################
date_firstrecord = '1979010100' # only if record_num > -1
date_lastrecord = '2018123118' # only if record_num > -1
map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
zoom = 'true'
proj_latlon = [90. , 270.]
hinc = 6 # number of hours between records;  only used if record_num > -1

plot_option = 2 # 1 for histogram, 2 for probability density function
plot_variable = 3 # 1 for max amplitude, 2 for median amplitude, 3 for minimum theta, 4 for minimum latitude, 5 for mean latitude, 6 for maximum latitude, 7 for radius, 8 for median year, 9 for total year points, 10 for genesis latitude, 11 for min trop theta anomaly, 12 for lifetime
stat_sig_opt = 0 # 2 for bootstrapping
climatology_file = 0 # 0 if input files do not contain climatology values, 1 otherwise
num_files = 2
legendtext = ['Description of File 1','Description of File 2','Description of File 3']
figname_description = 'file1_vs_file2'
make_tpv = 'False'
min_lifetime_days = 2
lat_origin_range = [-90, 90.]
label_fontsize = 18
npasses_smooth = 2
#confidence_interval_limits = [5,95]
confidence_interval_limits = [2.5,97.5]


fdir0 = '/Users/scavallo/Documents/data/track_files/'
fdir1 = '/Users/scavallo/Documents/data/track_files/'
textfile0 = '' # Enter name of track file 1
textfile1 = '' # Enter name of track file 2

imagedir = '/Users/scavallo/Documents/scripts/python_scripts/images/'
###################################
# END user options
###################################
R_earth = 6371200
R_earth = R_earth / 1000
pid = np.pi/180


plotvar_cntl = []
plotvar_exp1 = []
plotvar_exp2 = []
plotvar_exp3 = []
for ii in range(0,num_files):
    
    if ii == 0:
        textfile = textfile0
        fdir = fdir0
    elif ii == 1:
        textfile = textfile1
        fdir = fdir1
    elif ii == 2:
        fdir = fdir2
        textfile = textfile2
    elif ii == 3:
        fdir = fdir3
        textfile = textfile3
    
    fpath = fdir + textfile
    print(fpath)
    aa = np.loadtxt(fpath, skiprows=1)       
    datelist = aa[:,0]
    lat = aa[:,1]	
    lon = aa[:,2]
    thetamin = aa[:,3]
    if climatology_file == 1:
        thetaclim = aa[:,4]	    
        thetaamp = aa[:,5]
        vort_circ = aa[:,6]	    
        vort_radius = aa[:,7]   
    else:
        thetaamp = aa[:,4]
        vort_circ = aa[:,5]	    
        vort_radius = aa[:,6]

    if plot_variable == 8:
	    yyyy = []
	    for tt in range(0,len(datelist)):
		    datelist_strnow = str(int(datelist[tt]))
		    yyyy.append(int(datelist_strnow[0:4]))
	
    nrows, ncols = np.shape(aa)

    tinds = np.where(datelist>10000)
    lats_all = lat[tinds]
    lons_all = lon[tinds]
    Y = lats_all[0:-1]
    X = lons_all[0:-1]
    latf = lats_all[1:]
    lonf = lons_all[1:]

    dist_all = R_earth * np.arccos( np.sin(Y*pid) * np.sin(latf*pid) + np.cos(Y*pid) * np.cos(latf*pid) * np.cos((lonf - X)*pid));
    speed_all = (dist_all*1000.) / (hinc*3600.)
    speed_all = um.filter_numeric_nans(speed_all,120,float('NaN'),'high')

    ntracks_min = (24./hinc)*min_lifetime_days
    tt = 0
    tracknum = -1
    while tt < nrows:
        #print(tt)
        if datelist[tt] < 10000:    	        
            ntracks = datelist[tt]
            nhours = ntracks*hinc
            sind = tt+1
            eind = np.int(sind+ntracks-1)

            meanind = np.int(np.floor((sind+eind)/2))
            datesnow = datelist[meanind]
            datestrnow = str(datesnow)
            #print(sind,eind)	
            yyyynow = int(datestrnow[0:4])

            if make_tpv == 'True':
                latnow = lat[sind:eind]
                nlats_polar = np.where(np.abs(latnow)>=65.0)
                perc = (np.float(np.size(nlats_polar)) / np.float(np.size(latnow)))*100.
                lat_origin_range = [0.,90.]
                if min_lifetime_days < 2:
                    min_lifetime_days = 2.0
                ntracks_min = (24./hinc)*min_lifetime_days
                perc_thresh = 60.
				#print  perc
            else:
                perc = 0.0    	
                perc_thresh = -1.0        

            if ( (lat[sind] >= lat_origin_range[0]) and (lat[sind] <= lat_origin_range[1]) and (ntracks >= ntracks_min) and (perc >= perc_thresh) ):
                tracknum += 1
                latnow = lat[sind:eind]
                lonnow = lon[sind:eind]
                trthnow = thetamin[sind:eind]
                if climatology_file == 1:
                    trth_climnow = thetaclim[sind:eind]
                    trthanom = trthnow - trth_climnow
                ampsnow = thetaamp[sind:eind]	    
                radsnow = vort_radius[sind:eind]	    
                circsnow = vort_circ[sind:eind]	   
                
                #trthanom = trthnow - trth_climnow 

                if ii == 0:
                    if plot_variable == 1 :
                        plotvar_cntl.append(np.max(ampsnow))
                    elif plot_variable == 2 :
                        plotvar_cntl.append(np.median(ampsnow))
                    elif plot_variable == 3 :
                        plotvar_cntl.append(np.min(trthnow))		    	    
                    elif plot_variable == 4 :
                        plotvar_cntl.append(np.min(np.abs(latnow)))		    
                    elif plot_variable == 5 :
                        plotvar_cntl.append(np.mean(np.abs(latnow)))			
                    elif plot_variable == 6 :
                        plotvar_cntl.append(np.max(np.abs(latnow)))			
                    elif plot_variable == 7 :
                        plotvar_cntl.append(np.max(radsnow))			
                    elif plot_variable == 8 :
                        plotvar_cntl.append(np.max(yyyynow))		
                    elif plot_variable == 9 :
                        plotvar_cntl.append(yyyynow)			
                    elif plot_variable == 10 :
                        plotvar_cntl.append(np.abs(latnow[0]))
                    elif plot_variable == 11:
                        plotvar_cntl.append(np.nanmin(trthanom))
                    elif plot_variable == 12:
                        plotvar_cntl.append(nhours/24.)                        
                elif ii == 1:
                    if plot_variable == 1 :
                        plotvar_exp1.append(np.max(ampsnow))
                    elif plot_variable == 2 :
                        plotvar_exp1.append(np.median(ampsnow))
                    elif plot_variable == 3 :
                        plotvar_exp1.append(np.min(trthnow))		    	    
                    elif plot_variable == 4 :
                        plotvar_exp1.append(np.min(np.abs(latnow)))		    
                    elif plot_variable == 5 :
                        plotvar_exp1.append(np.mean(np.abs(latnow)))			
                    elif plot_variable == 6 :
                        plotvar_exp1.append(np.max(np.abs(latnow)))			
                    elif plot_variable == 7 :
                        plotvar_exp1.append(np.max(radsnow))			
                    elif plot_variable == 8 :
                        plotvar_exp1.append(np.max(yyyynow))		
                    elif plot_variable == 9 :
                        plotvar_exp1.append(yyyynow)			
                    elif plot_variable == 10 :
                        plotvar_exp1.append(np.abs(latnow[0]))
                    elif plot_variable == 11:
                        plotvar_exp1.append(np.nanmin(trthanom))
                    elif plot_variable == 12:
                        plotvar_exp1.append(nhours/24.)       
                elif ii == 2:
                    if plot_variable == 1 :
                        plotvar_exp2.append(np.max(ampsnow))
                    elif plot_variable == 2 :
                        plotvar_exp2.append(np.median(ampsnow))
                    elif plot_variable == 3 :
                        plotvar_exp2.append(np.min(trthnow))		    	    
                    elif plot_variable == 4 :
                        plotvar_exp2.append(np.min(np.abs(latnow)))		    
                    elif plot_variable == 5 :
                        plotvar_exp2.append(np.mean(np.abs(latnow)))			
                    elif plot_variable == 6 :
                        plotvar_exp2.append(np.max(np.abs(latnow)))			
                    elif plot_variable == 7 :
                        plotvar_exp2.append(np.max(radsnow))			
                    elif plot_variable == 8 :
                        plotvar_exp2.append(np.max(yyyynow))		
                    elif plot_variable == 9 :
                        plotvar_exp2.append(yyyynow)			
                    elif plot_variable == 10 :
                        plotvar_exp2.append(np.abs(latnow[0]))
                    elif plot_variable == 11:
                        plotvar_exp2.append(np.nanmin(trthanom))
                    elif plot_variable == 12:
                        plotvar_exp2.append(nhours/24.) 
                elif ii == 3:
                    if plot_variable == 1 :
                        plotvar_exp3.append(np.max(ampsnow))
                    elif plot_variable == 2 :
                        plotvar_exp3.append(np.median(ampsnow))
                    elif plot_variable == 3 :
                        plotvar_exp3.append(np.min(trthnow))		    	    
                    elif plot_variable == 4 :
                        plotvar_exp3.append(np.min(np.abs(latnow)))		    
                    elif plot_variable == 5 :
                        plotvar_exp3.append(np.mean(np.abs(latnow)))			
                    elif plot_variable == 6 :
                        plotvar_exp3.append(np.max(np.abs(latnow)))			
                    elif plot_variable == 7 :
                        plotvar_exp3.append(np.max(radsnow))			
                    elif plot_variable == 8 :
                        plotvar_exp3.append(np.max(yyyynow))		
                    elif plot_variable == 9 :
                        plotvar_exp3.append(yyyynow)			
                    elif plot_variable == 10 :
                        plotvar_exp3.append(np.abs(latnow[0]))
                    elif plot_variable == 11:
                        plotvar_exp3.append(np.nanmin(trthanom))
                    elif plot_variable == 12:
                        plotvar_exp3.append(nhours/24.) 

        tt+=1


mstats(plotvar_cntl)
mstats(plotvar_exp1)

golden = (pylab.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 8./ golden ), dpi=128)    
adjustprops = dict(left=0.15, bottom=0.1, right=0.90, top=0.93, wspace=0.2, hspace=0.2) 

if ( (plot_variable == 1) or (plot_variable == 2) ):
    binvals = np.arange(0,46,1)
    binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
    bincenters = binvals[0:-1]
    xlabeltext = 'Amplitude (K)'
    if plot_variable == 1:
        figname_prefix = 'tpv_maxamplitude_comparisons'
        ylims = [0,0.14]
        #ylims = [0,0.1]
    elif plot_variable == 2:
        figname_prefix = 'tpv_medianamplitude_comparisons'
        ylims = [0,0.22]
elif plot_variable == 3:
    binvals = np.arange(250,320,1)
    binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
    bincenters = binvals[0:-1]

    xlabeltext = 'Kelvin'
    figname_prefix = 'tpv_mintheta_comparisons'
    ylims = [0,0.08]
    #ylims = [0,0.06]
elif ( (plot_variable == 4) or (plot_variable == 5) or (plot_variable == 6) or (plot_variable == 10)):
    binvals = np.arange(22,91,1)
    binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
    bincenters = binvals[0:-1]    
    
    xlabeltext = r'Latitude ($^{\circ}$S)'
    ylims = [0,0.08]
    #ylims = [0,0.06]
    if plot_variable == 4:
        figname_prefix = 'tpv_minlat_comparisons'
    elif plot_variable == 5:
        figname_prefix = 'tpv_meanlat_comparisons'
    elif plot_variable == 6:
        figname_prefix = 'tpv_maxlat_comparisons'
    elif plot_variable == 10:
        figname_prefix = 'tpv_genesislat_comparisons'	
elif plot_variable == 7:
    binvals = np.arange(100,1301,50)
    binvalsh = np.arange(binvals[0]-25,binvals[-1]+25,50)
    bincenters = binvals[0:-1]    
    xlabeltext = 'Radius (km)'    
    figname_prefix = 'tpv_radius_comparisons'
    ylims = [0,0.003]
elif ( (plot_variable == 8) or (plot_variable == 9) ):
    binvals = np.arange(1981,2016,1)
    xlabeltext = 'Year'    
    ylims = [0,0.04]

    binvalsh = np.arange(1980.5,2016.5,1)
    #del bincenters
    bincenters = binvals
    if plot_variable == 8:
        figname_prefix = 'tpv_medianyear_comparisons'
    elif plot_variable == 9:
        figname_prefix = 'tpv_numyearpoint_comparisons'
elif ( plot_variable == 11) :
    binvals = np.arange(-60,1,1)
    binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,1)
    bincenters = binvals[0:-1]
    xlabeltext = 'Anomaly (K)'
    figname_prefix = 'tpv_maxanomaly_comparisons'
    ylims = [0,0.07]
    #ylims = [0,0.06]
elif ( plot_variable == 12) :
    binvals = np.arange(0,91,5)
    binvalsh = np.arange(binvals[0]-0.5,binvals[-1]+0.5,5)
    bincenters = binvals[0:-1]
    xlabeltext = 'Lifetime (Number of days)'
    figname_prefix = 'tpv_lifetime_comparisons'
    ylims = [0,0.14]

fig = plt.figure(**figprops)   # New figure
ax1 = fig.add_subplot(1, 1, 1)


if plot_option == 1:
   alphaval = 0.75
   normval = 0
elif plot_option == 2:
   alphaval = 0.0
   normval = 1
###################################
# Figure 1
###################################
#print(binvals)
plotvar_cntl = np.array(plotvar_cntl)
if num_files >= 2:
    plotvar_exp1 = np.array(plotvar_exp1)
if num_files >= 3:
    plotvar_exp2 = np.array(plotvar_exp2)
if num_files >= 4:
    plotvar_exp3 = np.array(plotvar_exp3)




if plot_option == 1:
    if num_files == 1:
        n, bins, patches = ax1.hist(plotvar_cntl, binvals, density=False, histtype='bar',color='b', alpha=alphaval, rwidth = 1.00, align='mid')
    else:      
        n, bins, patches = ax1.hist(plotvar_cntl,binvals, density=False, histtype='bar', color='b', alpha=alphaval, rwidth = 1.00, align='mid',label=legendtext[0])
        n2, bins2, patches2 = ax1.hist(plotvar_exp1, binvals, density=False, histtype='bar', facecolor='r', alpha=alphaval, rwidth = 0.50, align='mid',label=legendtext[1])
    
    
elif plot_option == 2:
    if num_files == 1:
        n, bins  = np.histogram(plotvar_cntl, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
    elif num_files == 2:
        n, bins  = np.histogram(plotvar_cntl, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
        n2, bins2 = np.histogram(plotvar_exp1, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)     
    elif num_files == 3:
        n, bins  = np.histogram(plotvar_cntl, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
        n2, bins2 = np.histogram(plotvar_exp1, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)     
        n3, bins3 = np.histogram(plotvar_exp2, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)     
    else:
        n, bins  = np.histogram(plotvar_cntl, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
        n2, bins2 = np.histogram(plotvar_exp1, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)     
        n3, bins3 = np.histogram(plotvar_exp2, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None) 
        n4, bins4 = np.histogram(plotvar_exp3, bins=binvalsh, range=(binvals[0],binvals[-1]), normed=normval,weights=None, density=None)
    
    n[1:] = um.smooth_onedim(n[1:],npasses_smooth) 
    if num_files >= 2:    
        n2[1:] = um.smooth_onedim(n2[1:],npasses_smooth)     
    if num_files >= 3:
        n3[1:] = um.smooth_onedim(n3[1:],npasses_smooth) 
    if num_files >= 4:
        n4[1:] = um.smooth_onedim(n4[1:],npasses_smooth) 
    
    ndiff = n2 - n
    nzeros = np.zeros_like(ndiff).astype('f')   

    if stat_sig_opt == 1:
        ul =  np.percentile(ndiff, confidence_interval_limits[1])
        ll =  np.percentile(ndiff, confidence_interval_limits[0])
        err1 = ndiff-ll
        err2 = ndiff-ul
    elif stat_sig_opt == 2:
        alpha = (confidence_interval_limits[0]*2.0)/100.
        #boot_means = wm.bootstrap(ndiff,10000)
        boot_means,ci_limits_out = wm.bootstrap_twosamps(n,n2,10000,alpha)
        ci_limits = [100.0*alpha/2.0,100.0*(1.0-alpha/2.0)]
        print(ci_limits)

        ci_upper = np.percentile(boot_means,ci_limits[0])
        ci_lower = np.percentile(boot_means,ci_limits[1])
        print(ci_limits_out)
      
        #if ci_upper > ci_lower:
        #    ci_limits_out = [ci_lower,ci_upper]    
        #else:
        #    ci_limits_out = [ci_upper,ci_lower]      
        #err1 = ndiff-ci_limits_out[0]
        #err2 = ndiff-ci_limits_out[1]
        err1 = n-ci_limits_out[0]
        err2 = n-ci_limits_out[1]
        
        err2a = n2-np.abs(ci_limits_out[0])
        err2b = n2+ci_limits_out[1]
        
        if num_files >= 3:
            err3a = n3-np.abs(ci_limits_out[0])
            err3b = n3+ci_limits_out[1]            
        if num_files >= 4:
            err4a = n4-np.abs(ci_limits_out[0])
            err4b = n4+ci_limits_out[1]                         
           
    p1, = ax1.plot(bincenters,n,'k',linewidth=6.0,label=legendtext[0])  
    if num_files >= 2: 
        p2, = ax1.plot(bincenters,n2,'g',linewidth=3.0,label=legendtext[1])
    if num_files >= 3:
        p3, = ax1.plot(bincenters,n3,'orange',linewidth=3.0,label=legendtext[2])   
    if num_files >= 4:
        p4, = ax1.plot(bincenters,n4,'0.8',linewidth=3.0,label=legendtext[3])   
    
    #ksstat,pval = stats.ks_2samp(n,n2)
    #print(ksstat,pval)
    ksstat1,pval1 = stats.ks_2samp(plotvar_cntl,plotvar_exp1)
    print("K-S test between control and sample 1: ", ksstat1,pval1) 
    #print("K-S test between control and sample 1 is %7.2f " % (ksstat1,pval1)) 
    if num_files >= 3:
        ksstat2,pval2 = stats.ks_2samp(plotvar_cntl,plotvar_exp2)
        print("K-S test between control and sample 2: ", ksstat2,pval2) 
        
        ksstat2b,pval2b = stats.ks_2samp(plotvar_exp1,plotvar_exp2)
        print("K-S test between sample 1 and sample 2: ", ksstat2b,pval2b)         
    #ksstat3,pval3 = stats.ks_2samp(plotvar_cntl,plotvar_exp3)
    #print("K-S test between control and sample 1: ", ksstat1,pval1)    
    #print("K-S test between control and sample 2: ", ksstat2,pval2)  
    #print("K-S test between control and sample 3: ", ksstat3,pval3)  
    
    if ( stat_sig_opt > 0):
    #if 1 == 1:
        ax1.fill_between(bincenters, err1, err2, alpha=0.5, edgecolor='0.5', facecolor='0.5')
        ax1.fill_between(bincenters, err2a, err2b, alpha=0.5, edgecolor='b', facecolor='b')	
        if num_files >= 3:
            ax1.fill_between(bincenters, err3a, err3b, alpha=0.5, edgecolor='r', facecolor='r')	
        if num_files >= 4:
            ax1.fill_between(bincenters, err4a, err4b, alpha=0.5, edgecolor='g', facecolor='g')	

        #upper_err = ci_limits_out[0] - np.zeros(len(n2))
        #lower_err = np.zeros(len(n2)) - ci_limits_out[1]
        #ax1.errorbar(bincenters,n2,yerr=[lower_err,upper_err])
        
        
    #else:
        #p1, = ax1.plot(bincenters,n,'r',linewidth=3.0)   	


ax1.grid(True, linestyle='-')
if ( (plot_variable == 3) ):
    legend = ax1.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
elif (  (plot_variable == 1)  or (plot_variable == 2) or (plot_variable == 3)  or (plot_variable == 7) ):
    legend = ax1.legend(loc='upper right', shadow=True, fontsize=label_fontsize)
elif  ( (plot_variable == 4) or (plot_variable == 5) or (plot_variable == 6) ):
    legend = ax1.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
elif  ( (plot_variable == 8) or (plot_variable == 9) ):
    legend = ax1.legend(loc='lower left', shadow=True, fontsize=label_fontsize)
elif plot_variable == 11:
   legend = ax1.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
elif plot_variable == 12:
   plt.ylim([0,0.08]) 
   legend = ax1.legend(loc='upper right', shadow=True, fontsize=label_fontsize)
else:
    legend = ax1.legend(loc='upper left', shadow=True, fontsize=label_fontsize)
#plt.title(titletext1,fontsize=label_fontsize)

plt.xlim([binvals[0],binvals[-1]])
if plot_option == 1:
    plt.xlim([binvals[0],binvals[-1]]) 
    xticks_all = binvals[:]+0.5     
    xtickcint = 2
    xticks = xticks_all[::xtickcint]    
        
    xticklabels = binvals[::xtickcint]
      
if plot_option == 2:        
    
    xtickcint = 2
    if plot_variable == 7:    
        xtickcint = 100 
    if plot_variable == 12:
        xtickcint = 5
    plt.xlim([bincenters[0],bincenters[-1]])
    xticks = np.arange(binvals[0],binvals[-1],xtickcint)  
    xticklabels = np.arange(binvals[0],binvals[-1],xtickcint)  



ax1.xaxis.set_ticks(xticks)
ax1.set_xticklabels(xticklabels, rotation = 90)
ax1.set_ylim(ylims[0],ylims[-1])
if plot_option == 1:
    plt.ylabel('Number',fontsize=label_fontsize)
    figname_suffix = 'histogram'
elif plot_option == 2:
    plt.ylabel('Probability',fontsize=label_fontsize)
    figname_suffix = 'pdf'
plt.xlabel(xlabeltext,fontsize=label_fontsize)
save_name = figname_prefix + '_' + figname_description + '_'  + figname_suffix + '_' + date_firstrecord +  '_' + date_lastrecord + '.png'
plt.savefig(imagedir + save_name, bbox_inches='tight')

plt.show()