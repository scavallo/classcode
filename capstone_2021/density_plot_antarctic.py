##################################################################
# density_plot_antarctic                        
# Created: 12/16/2020                                              
#                                                                
# Steven Cavallo                                                  
##################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import netCDF4
from scipy import ndimage
from netCDF4 import Dataset
from scipy.stats import gaussian_kde
import utilities_modules as um
from mstats import *
##################################################################
# User Settings
##################################################################

fSave = '/home/scavallo/scripts/python_scripts/images/' #Directory to save plot
fPlot = 'density_tpv_antarctic.png'
titletext = ''
fLoad = '/data2/scavallo/era_interim/sh_tracks/density_array_sh_60N65.npy'
fLoad2 = '/data2/scavallo/era_interim/sh_tracks/density_array_sh_60N65.npy'

nyears_tpv_count = 2020-1979+1
nyears_tpv = [nyears_tpv_count,nyears_tpv_count]
save = True #True = save fig without showing on screen, #False = show fig

plot_probability = False
plot_tpvsperyear = False
plot_logspace = False

norm = False #True = normalize by largest value
addcolm = False #True = fix my dumb arrays with more dumb cause I'm dumb
diff = False # If True, takes fLoad2 minus fLoad
subtract_mean = False # True to subtract the mean
subtract_mean_standardize = False # Same as subtract_mean, but standardizes instead
plot_atmos = False # True to plot fLoad_atm fields 
num_smooth_iterations = 1


hemisphere = 'southern'
map_projection = 'spstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
proj_latlon = [-50. , 180.]
zoom = 'false'
cen_long = 180.
label_fontsize = 22

##################################################################
# Create the tracks
##################################################################
if (plot_atmos == True):
    infile =  Dataset( fLoad_atm, 'r')
    lat_atm = infile.variables['lat'][:]
    lon_atm = infile.variables['lon'][:]
    ghgt_cntl = infile.variables['ghgt_cntl'][:]
    ghgt_exp = infile.variables['ghgt_exp'][:]
    slp_cntl = infile.variables['slp_cntl'][:]
    slp_exp = infile.variables['slp_exp'][:]     
    trth_cntl = infile.variables['trth_cntl'][:]
    trth_exp = infile.variables['trth_exp'][:]   
    infile.close()
    
    ghgt_diff = ghgt_exp - ghgt_cntl
    slp_diff = slp_exp - slp_cntl
    trth_diff = trth_exp - trth_cntl
    
    lonin = lon_atm
    ghgt_diff, lons_atm = um.addcyclic(ghgt_diff, lonin)
    trth_diff, lons_atm = um.addcyclic(trth_diff, lonin)
    slp_diff, lons_atm = um.addcyclic(slp_diff, lonin)
    ghgt_exp, lons_atm = um.addcyclic(ghgt_exp, lonin)
    trth_exp, lons_atm = um.addcyclic(trth_exp, lonin)
    slp_exp, lons_atm = um.addcyclic(slp_exp, lonin)    
    ghgt_cntl, lons_atm = um.addcyclic(ghgt_cntl, lonin)
    trth_cntl, lons_atm = um.addcyclic(trth_cntl, lonin)
    slp_cntl, lons_atm = um.addcyclic(slp_cntl, lonin)        
    Xatm, Yatm = np.meshgrid(lons_atm, lat_atm)    
            
    ghgt_plot = ghgt_exp
    trth_plot = trth_diff
    slp_plot = slp_diff	    
	    
    cint_ghgt = 60 #5
    cbar_min_ghgt = 5400 - (cint_ghgt*20)#-100
    cbar_max_ghgt = 5400 + (cint_ghgt*20)#100
    cmap_opt = plt.cm.RdBu_r
    cflevs_ghgt =  np.arange(cbar_min_ghgt, cbar_max_ghgt+(cint_ghgt/2), cint_ghgt)
    cflevs_ghgt_cntrs = cflevs_ghgt[::2]
    cflevs_ghgt_ticks = np.arange(cbar_min_ghgt,cbar_max_ghgt+(cint_ghgt/2),4*cint_ghgt)

    cint_trth = 1
    cbar_min_trth = -20
    cbar_max_trth = 20
    cmap_opt_trth = plt.cm.RdBu_r
    cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth+(cint_trth/2), cint_trth)
    cflevs_trth_cntrs = cflevs_ghgt[::2]
    cflevs_trthticks = np.arange(cbar_min_trth,cbar_max_trth+(cint_trth/2),4*cint_trth)


    cint_slp = 0.5#2
    cbar_min_slp = -30 #1004-(20*cint_slp)
    cbar_max_slp = 30 #1004+(20*cint_slp)
    cflevs_slp =  np.arange(cbar_min_slp, cbar_max_slp+(cint_slp/2), cint_slp)
    #cflevs_slp = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]    


if (subtract_mean_standardize == True):
    subtract_mean = False    


if (diff == False):   
   infile = np.load(fLoad)   
   den_arr = infile['den_arr']
   count_arr = infile['count_arr']
   hit_arr = infile['hit_arr']

   mstats(den_arr)

   tpvs_peryear = hit_arr / (nyears_tpv[0])
   if (plot_tpvsperyear == True) :
       den_arr = tpvs_peryear
   
   if (subtract_mean == True ):
       den_arr = den_arr - np.nanmean(den_arr)

   if (subtract_mean_standardize == True ):
       den_arr = (den_arr - np.nanmean(den_arr))/np.std(den_arr)
             
else:
   infile1 = np.load(fLoad)
   infile2 = np.load(fLoad2)
   den_arr1 = infile1['den_arr']
   den_arr2 = infile2['den_arr']
   count_arr = infile1['count_arr']
   hit_arr1 = infile1['hit_arr']   
   hit_arr2 = infile2['hit_arr']   
      
   if (subtract_mean == True):
       den_arr2 = den_arr2 - np.nanmean(den_arr2)
       den_arr1 = den_arr1 - np.nanmean(den_arr1)
   if (subtract_mean_standardize == True):      
       den_arr2 = (den_arr2 - np.nanmean(den_arr2))/np.std(den_arr2)
       den_arr1 = (den_arr1 - np.nanmean(den_arr1))/np.std(den_arr1)  
   
   try:
       den_arr = den_arr2 - den_arr1
   except:
       new_colm = den_arr1[-1,:]
       [ax,ay] = np.shape(den_arr2)
       den_arr1 = np.insert(den_arr1,ax-1,new_colm,axis=0)
       den_arr = den_arr2 - den_arr1
       
if (plot_probability == True):
    if (diff == False) :
        #den_arr = (hit_arr / count_arr)*100.
        den_arr = (tpvs_peryear/365.)*100
        #den_arr = (tpvs_peryear/90.)*100
	#den_arr = (tpvs_peryear/31.)*100
	#den_arr = ( (tpvs_peryear/90.) )
	    
    else:
        tpvs_peryear1 = hit_arr1 / (nyears_tpv[0])
	tpvs_peryear2 = hit_arr2 / (nyears_tpv[1])
	
	#den_arr = ((tpvs_peryear2/90.)*100) - ((tpvs_peryear1/90.)*100)
	#den_arr = ((tpvs_peryear2/90.)*100) - ((tpvs_peryear1/365.)*100)
	#den_arr = ((tpvs_peryear2/31.)*100) - ((tpvs_peryear1/365.)*100)
	den_arr = ((tpvs_peryear2/90.)) - ((tpvs_peryear1/90.))
	#den_arr = ((tpvs_peryear2/31.)*100) - ((tpvs_peryear1/31.)*100)
	   
if (norm == True):
   den_arr = den_arr / np.amax(den_arr)
   
if (addcolm == True):
    mstats(den_arr)
    #new_colm = den_arr[:,0]
    #den_arr = np.insert(den_arr,72,new_colm,axis=1)
    new_colm = den_arr[-1,:]
    [ax,ay] = np.shape(den_arr)
    den_arr = np.insert(den_arr,ax-1,new_colm,axis=0)

    latboxsize = 2
    lonboxsize = 5

    if hemisphere == 'northern':
        xlat = np.arange(30,90+(latboxsize/2),latboxsize)
        xlon = np.arange(0,365,lonboxsize)
    else:
        xlat = np.arange(-90,-30+(latboxsize/2),latboxsize)
        xlon = np.arange(0,365,lonboxsize)   

    xloni,xlati = np.meshgrid(xlon,xlat)
    
    
else:

##################################################################
# THIS MESH MUST ALIGN WITH MESH USED FOR NUMPY ARRAY WHEN CREATED
#
# Standard settings used: lat=2 lon=5
##################################################################
    
    latboxsize = 2
    lonboxsize = 5

    if hemisphere == 'northern':
        xlat = np.arange(30,90+(latboxsize/2),latboxsize)
        xlon = np.arange(0,365,lonboxsize)
    else:
        xlat = np.arange(-90,-30+(latboxsize/2),latboxsize)
        xlon = np.arange(0,365,lonboxsize)   
    xloni,xlati = np.meshgrid(xlon,xlat)    

if num_smooth_iterations > 0: 
    smooth_amp = 1.25
    for ii in range(0,num_smooth_iterations):
        den_arr[:,1:-1]  = ndimage.gaussian_filter(den_arr[:,1:-1],smooth_amp)
	if (plot_atmos == True):
	    ghgt_plot = ndimage.gaussian_filter(ghgt_plot,smooth_amp)
	    trth_plot = ndimage.gaussian_filter(trth_plot,smooth_amp)
            slp_plot = ndimage.gaussian_filter(slp_plot,smooth_amp)

if (diff == False):   
   min_den = 10
   max_den = 400
   cint_den = 10  
   cmap_opt = plt.cm.hot_r
   cbarlabel = 'TPVs within 555 km of point over record'
else:
   min_den = -1500
   max_den = 1500
   cint_den = 100
   cmap_opt = plt.cm.RdBu_r
   cbarlabel = 'TPV difference within 555 km of point over record'
   

   
if (plot_probability == True):   
   if (diff == False):
       min_den = 0
       max_den = 2
       cint_den = 0.1  
       cmap_opt = plt.cm.gist_heat_r
   else:
       min_den = -10
       max_den = 10
       cint_den = 0.5             
       cmap_opt = plt.cm.RdBu_r
   cbarlabel = 'TPVs per day'
if (plot_tpvsperyear == True) :
   min_den = 10
   max_den = 150
   cint_den = 5   
   cmap_opt = plt.cm.gist_heat_r
   cbarlabel = 'TPVs year-1'
if norm == True:
   min_den = 0
   max_den = 1
   cint_den = 0.05   
   cmap_opt = plt.cm.gist_heat_r
   cbarlabel = 'Normalized TPVs'
if (subtract_mean_standardize == True):
   min_den = -2.0
   max_den = 2.0
   cint_den = 0.1  
   cbarlabel = 'Number of standard deviations'        


if plot_logspace == True:
    cflevs_den =  np.logspace(0, 2, 30)
else:
    cflevs_den =  np.arange(min_den, max_den+(cint_den/2), cint_den)

cflevs_den_ticks = cflevs_den
golden = (np.sqrt(5)+1.)/2.
####################################################
# Figure 1
####################################################
fig = plt.figure(figsize=(8., 16./golden), dpi=128)  
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
if map_projection == 'ortho':
   if zoom == 'false':   
      m = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
               resolution = 'l', area_thresh = 1000.,ax=ax1)
   else:
      m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
               resolution = 'l', area_thresh = 1000.,ax=ax1)      		           

      width = m1.urcrnrx - m1.llcrnrx
      height = m1.urcrnry - m1.llcrnry

      #coef = 0.5
      coef = 0.7
      width = width*coef
      height = height*coef
      m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
          llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)		  

elif map_projection == 'lcc':
    m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax1)	
elif ( (map_projection == 'npstere') or (map_projection == 'spstere') ) :            
    m = Basemap(projection=map_projection,boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
F = plt.gcf()  # Gets the current figure

m.drawstates(color='#444444', linewidth=1.25)
m.drawcoastlines(color='#444444')
m.drawcountries(color='#444444', linewidth=1.25)

x, y = m(xloni,xlati)
mstats(den_arr)
if ( plot_atmos == True ): 
    xatm, yatm = m(Xatm, Yatm)
CS1 = m.contourf(x,y, den_arr, cmap=cmap_opt,levels=cflevs_den, extend='both',zorder=1)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)   
labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
cbar.ax.set_xticklabels(labels, size=16)
cbar.set_label(cbarlabel,size=20)
if (plot_probability == True):

    if (diff == True):
	levels = [0.,1.,5.,10.,20.,30.,40.,50.]
	CS2 = m.contour(x, y, den_arr, levels, colors='0.75', linewidths=1.25)      
	plt.clabel(CS2, levels, fmt = '%i', inline=True, fontsize=10)    
        levels2 = [-1.,-5.,-10.,-20.,-30.,-40.,-50.]
        CS3 = m.contour(x, y, den_arr, levels2, colors='0.75', linestyle='--',linewidths=1.25)      
        plt.clabel(CS3, levels2, fmt = '%i', inline=True, fontsize=10)	
    else:
	#levels = [1.,5.,10.,15.,20.]
	levels = [0.,1.,5.,10.,20.,30.,40.,50.]
	CS2 = m.contour(x, y, den_arr, levels, colors='0.75', linewidths=1.25)      
	plt.clabel(CS2, levels, fmt = '%i', inline=True, fontsize=10)    
if (plot_tpvsperyear == True) :    
    levels = [1.,10.,50.,100.,150.,200.]
    CS2 = m.contour(x, y, den_arr, levels, colors='0.75', linewidths=1.25)      
    plt.clabel(CS2, levels, fmt = '%i', inline=True, fontsize=10)
#plt.title('Probability per day (ERA Interim: DJF - Full year)',fontsize=28)
plt.title(titletext,fontsize=28)
if (plot_atmos == True):
    CS2 = m.contour(xatm, yatm, ghgt_exp, cflevs_ghgt, colors='0.75', linestyle='--',linewidths=1.5)
    plt.clabel(CS2,inline=1,fmt='%d',fontsize=10)
    CS3 = m.contour(xatm, yatm, ghgt_cntl, cflevs_ghgt, colors='0.25', linewidths=1.5)
    plt.clabel(CS3,inline=1,fmt='%d',fontsize=10)    
    

if (save == True):
   plt.savefig(fSave + fPlot, bbox_inches='tight'); #plt.close()

#if (save == False):



plt.show()

