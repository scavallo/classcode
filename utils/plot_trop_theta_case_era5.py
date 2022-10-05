#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab,sys
import numpy as np
import matplotlib as mpl
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy import ndimage
from matplotlib.colors import LogNorm

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *


import warnings
warnings.filterwarnings("ignore")
###################################
# Set user options
###################################
print("Python version")
print(sys.version)

#dates2plot = ['2006021800','2006031718']
#dates2plot = ['2004030100','2004040218']
#dates2plot = ['2021121500','2022013118']
#dates2plot = ['2020010100','2020012618']
#dates2plot = ['2011081300','2011091618']
dates2plot = ['2010041500','2010050918']
date_firstrecord = '2010040100'
date_firstrecord_sfcfile = '2010040100'
hinc = 6


map_projection = 'spstere' # 'npstere' for northern hemisphere polar stereorgraphic, 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
#proj_latlon = [65. , 0.]
#proj_latlon = [55. , 270.]
#proj_latlon = [70. , 270.]
#proj_latlon = [20. , 270.]
proj_latlon = [-30. , -180.]

#proj_latlon = [90. , 80.]
#proj_latlon = [90. , 310.]
plot_field = 'trth_ivt' # 'trth', 'trpr', 'trtemp','pw', 'trwind', 'sfcwind', 'totqv', 'quivwind','ghgt500', 'trth_slp_ice', 'trth_totqv', 'cloud_sw', 'cloud_lw','cloud_net','trth_ivt', 'blank' for blank map
plot_slp = 'True' # If true, plots slp contours
convert_to_celsius = 'False'
plot_trth_contour = 'False'
plot_trth_contours = 'True'
plot_cloud_cover = 'False'
smooth_plotvar = 'False'
plot_radius = 'False' # If true, plots radial circle
plot_wind = 'False'
plot_wind_level = -1
plot_spc_reports = 'false'
#fpath_spc = '/data1/scavallo/data/1950_2014_onetor.csv'
fpath_spc = '/data1/scavallo/data/1950-2016_all_tornadoes.csv'
plot_cross_section_line = 'false'
plot_radiosondes = 'True'
cenLat = 44.0 # for plot_cross_section_line only
cenLon = 286. # for plot_cross_section_line only
dLat = 0.     # for plot_cross_section_line only
dLon = 50.    # for plot_cross_section_line only

#latlon_radius = [35.77, -65.0] # Atlantic Ocean point
#latlon_radius = [35.77, -78.64] # Raleigh, NC
#latlon_radius = [41.26, -95.94] # Omaha, NE
#latlon_radius = [47.68, -116.78] # Coeur d'Alene
latlon_radius = [78.22, 15.63] 
zoom = 'False'
coef = 0.7

figname = 'era5_reanalysis'

russia_fir_lats= [60.0,90.,60.0]
russia_fir_lons = [30.0,30.0,191.6]

mosaic_lats = [82.,88.,85.,80.,75.]
mosaic_lons = [120.,120.,0.,0.,350.]

#mission3_lats = [78.22,78.22]
#mission3_lons = [15.62,15.62]

#mission1_lats = [78.22, 82., 76., 76., 82.,80., 86.,88.,85., 82., 72., 77., 82., 82.,82., 80.,78.22]
#mission1_lons = [15.62,225.,219.,195.,225.,302.,15.,05.,200.,225.,242.,268.,225.,195.,225.,255.,15.62]
#mission1_lats = [78.22, 82., 76., 76., 82., 81., 75., 82., 85., 82., 82., 80.,78.22]
#mission1_lons = [15.62,225.,219.,195.,225.,295.,240.,225.,200.,195.,225.,250.,15.62]
mission1_lats = [78.4, 82., 76., 76., 82., 81., 75., 82., 85.,78.4]
mission1_lons = [15.78,225.,219.,195.,225.,295.,240.,225.,200.,15.78]

mission3_lats = [78.4,87.,78.4]
mission3_lons = [15.78,10.,15.78]

mission2_lats = [78.4, 89., 74., 74., 80., 80., 86., 82.,86.,78.4]
mission2_lons = [15.78,300.,195.,240.,255.,195.,210.,298.,24.,15.78]

if 1 == 0:
    lmower_numlegs = 6
    lmower_latlon_start = [74.,240.]
    lmower_nlats = 2.0
    lmower_londist_km = 1500.

    lmower_latdist_km = lmower_nlats*111.
    lmower_latinc = lmower_latdist_km / 111.00
    mission2_lats = np.zeros(lmower_numlegs*2+2)
    mission2_lons = np.zeros(lmower_numlegs*2+2)
    mission2_lats[0] = latlon_radius[0]
    mission2_lons[0] = latlon_radius[1]
    mission2_lats[1] = lmower_latlon_start[0]
    mission2_lons[1] = lmower_latlon_start[1]
    mission2_lats[-1] = latlon_radius[0]
    mission2_lons[-1] = latlon_radius[1]
    count = 0
    signnow = -1.0
    for ii in range(2,len(mission2_lats)-1):
        if np.mod(count,2) != 0:
            mission2_lats[ii] = mission2_lats[ii-1]+lmower_latinc
        else:
            mission2_lats[ii] = mission2_lats[ii-1]    
        if np.mod(count,2) == 0:
            lmower_loninc = signnow*lmower_londist_km/(111.0*np.cos(mission2_lats[ii]*np.pi/180.))    
            signnow = -1.0*signnow
            mission2_lons[ii] = mission2_lons[ii-1]+lmower_loninc               
            if ( (mission2_lons[ii] < russia_fir_lons[2]) and (mission2_lons[ii] > russia_fir_lons[1]) ) :
               closest_lon = min(russia_fir_lons,key=lambda x:abs(x-mission2_lons[ii]))
               mission2_lons[ii]  = closest_lon + 3.0 #russia_fir_lons[2]

        else:
            mission2_lons[ii] = mission2_lons[ii-1]
        count += 1

datenow = date_firstrecord
count = 0
while datenow <= dates2plot[1]:
    #print(datenow)
    if datenow == dates2plot[0]:
        rec1 = count
    if datenow == dates2plot[1]:
        rec2 = count    
    datenow = um.advance_time(datenow,hinc)
    count+=1
record_num = [rec1, rec2]
print(record_num)

if record_num[0] == record_num[1]:
    loop_option = 'false'
else:    
    loop_option = 'true'

#data_gridnum = 4
data_gridnum = 20 # 3 = ncep
                  # 4 = ncep
          # 11 = era interim
          # 12 = era interim, alternative (2 PVU surface only)
          # 13 = era interim, from arctic server
	  # 20 = ERA5

label_fontsize = 16

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
#fdir = '/arctic3/datasets/ERA5/'
fdir = '/data1/scavallo/data/cases/ant_aprmay2010/'
#fdir = '/data1/scavallo/data/cases/ant_augsept2011/'
#fdir = '/data1/scavallo/data/cases/ant_janfeb2020/'
#fdir = '/data1/scavallo/data/cases/ant_marapr2006/'
#fdir = '/data1/scavallo/data/cases/arc_marapr1983/'
force_filename = True
#fname = 'era5_pt_1979-2018.nc'
fname = 'era5_2pvu_2010040100_2010053118.nc'
fname_sfc = 'era5_sfclevel_2010040100_2010053118.nc'
imagedir = '/home/scavallo/scripts/python_scripts/images/'
fpath_igra = '/home/scavallo/scripts/template/igra_stationlist.txt'
###################################
# END user options
###################################

# Open the netcdf file and read select variables
#dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
#fpath = os.path.join(fpath, fname)
#print fpath

lat_range = [cenLat-dLat, cenLat+dLat]
lon_range = [cenLon-dLon,cenLon+dLon]

if plot_radiosondes == 'True':
    ###################################
    # Get IGRA data
    ###################################
    f = open(fpath_igra, 'r')

    rad_lats = []
    rad_lons = []
    rad_endyear = []
    for line in f:
        line = line.strip()
        columns = line.split()

        #rad_lats.append(columns[1])
        #rad_lons.append(columns[2])    
        #rad_endyear.append(float(columns[-2]))
	
        rad_lats = np.append(rad_lats,columns[1])
        rad_lons = np.append(rad_lons,columns[2])
        rad_endyear = np.append(rad_endyear,columns[-2])

        rad_lats=[float(x) for x in rad_lats]
        rad_lons=[float(x) for x in rad_lons]    
        rad_endyear=[float(x) for x in rad_endyear] 



        rad_lats = um.filter_numeric_nans(rad_lats,90.0,float('NaN'),'both')
        rad_lons = um.filter_numeric_nans(rad_lons,-180.0,float('NaN'),'low')
        rad_endyear = um.filter_numeric_nans(rad_endyear,2100,float('NaN'),'high')

    
    [inds] = np.where(rad_endyear > 2013)

    tmp = rad_lats[inds]
    del rad_lats
    rad_lats = tmp
    del tmp

    tmp = rad_lons[inds]
    del rad_lons
    rad_lons = tmp
    del tmp

#########################
nplots = record_num[1] - record_num[0]


if loop_option == 'true':
   record_start = record_num[0]
   record_end = record_num[1]
else:
   record_start = record_num[0]
   record_end = record_num[0]

rcount = 0 
record_numnow = record_num[0] 
count_anal = record_numnow
print(record_start, record_end, nplots)
while record_start <= record_end:   
   #record_numnow = record_start
   print("rcount is %i while record_numnow is %i" %(rcount, record_numnow))

   if nplots > 1:
      count = 1
      
      if rcount == 0:
         datenow = dates2plot[0]     
      #while count <= record_num:
      # datenow = um.advance_time(datenow,hinc)      
      # count += 1


      date_string = datenow
      print(date_string)
      dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
     
      print(fname)
      fpath = fdir + fname
      
      
      fpath_sfc = fdir + fname_sfc
      if force_filename == False:
          print(fname_sfc)
          fpath_sfc = fdir + fname_sfc
      
      #if rcount == 0:
      #   fpath = os.path.join(fpath, fname)
      print(fpath)

   else:
      datenow = date_firstrecord
      #date_string = datenow
      hrs_ahead = (record_num[0]*hinc) - hinc
      date_string = um.advance_time(date_firstrecord,hrs_ahead)  
      dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
      #fpath = os.path.join(fpath, fname)
      if ( (data_gridnum == 3) or (data_gridnum == 4) or (data_gridnum == 12) ):
          #fname = 'erainterim_tropopause_' + datenow + '.nc'
          fname = fname
      elif force_filename == True:
          fname = fname
      else:
          fname = 'erainterim_pv_' + datenow + '.nc'
      #fname_sfc = 'erainterim_pressure_' + datenow + '.nc'
      fpath = fdir + fname
      fpath_sfc = fdir + fname_sfc
      print(fpath)
      print(fpath_sfc)


   print(fpath)
   f = netCDF4.Dataset(fpath,'r')
   if data_gridnum == 4:
      trtempin = f.variables['TMP_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      trpresin = f.variables['PRES_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      #uin = f.variables['UGRD_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      #vin = f.variables['VGRD_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      uin = f.variables['UGRD_P0_L100_GLL0'][record_numnow,19,::-1,:].squeeze()
      vin = f.variables['VGRD_P0_L100_GLL0'][record_numnow,19,::-1,:].squeeze()
      slp = f.variables['PRMSL_P0_L101_GLL0'][record_numnow,::-1,:].squeeze()/10**2
      
      lons = f.variables['lon_0'][:]
      lats = f.variables['lat_0'][::-1] # Read in reverse direction
      levs = f.variables['lv_ISBL0'][:]  
      print(levs)    
    
   elif data_gridnum == 11:
      f2 = netCDF4.Dataset(fpath_sfc,'r') 
      

      trth = f.variables['th'][record_numnow,3,::-1,:].squeeze()      
      uin = f.variables['u'][record_numnow,3,::-1,:].squeeze()
      vin = f.variables['v'][record_numnow,3,::-1,:].squeeze()
      if plot_field == 'sfcwind':         
         uin_sfc = f2.variables['u10'][record_numnow,:,:].squeeze()
         vin_sfc = f2.variables['v10'][record_numnow,:,:].squeeze()               
      if plot_field == 'ghgt500':         
         levels =  f2.variables['level'][:]
         [ind500] = np.where(levels==500)
         ghgt500 = (f2.variables['z'][record_numnow,ind500,:,:].squeeze())/9.81
     
             
      slp = (f2.variables['msl'][record_numnow,:,:].squeeze())/100.     
      trpresin = f.variables['p'][record_numnow,3,::-1,:].squeeze()
      
      xlons = f.variables['xlong'][::-1,:]
      lons = f.variables['xlong'][0,:]
      lats = f.variables['xlat'][::-1,0] 
      xlats = f.variables['xlat'][::-1,:] 
      levs = f.variables['levels'][:]             
      
      inds = np.where(trth>400)
      uin[inds] = float('NaN');
      vin[inds] = float('NaN');  
      
      if plot_field == 'pw':
          spechum = (f2.variables['q'][record_numnow,:,:,:].squeeze()) 
      tempin = (f2.variables['t'][record_numnow,:,:,:].squeeze()) 
      hght = (f2.variables['z'][record_numnow,:,:,:].squeeze())/9.81 
      pin = f2.variables['level'][:]*100         
      qv = (spechum/(1-spechum))    

      pres = np.zeros_like(tempin).astype('f')   
      for kk in range(0,len(levs)):      
          pres[kk,:,:] = pin[kk]     
          
      qvtot = wm.total_col(qv[::-1,:,:], pres[::-1,:,:], tempin[::-1,:,:], hght[::-1,:,:])         
      if ( (plot_field == 'totqv') or (plot_field == 'trth_totqv') ):
          qvtot = (f2.variables['q'][record_numnow,::-1,:].squeeze()) 
      
      #f2.close
   elif data_gridnum == 12:
      lons = f.variables['longitude'][:]
      lats = f.variables['latitude'][::-1]          
      if( (plot_slp == 'True') or (plot_field == 'trth_totqv') ):
          f2 = netCDF4.Dataset(fpath_sfc,'r') 
          sfclats = f2.variables['latitude'][::-1]
          sfclons = f2.variables['longitude'][:]
          xlons_sfc,xlats_sfc = np.meshgrid(sfclons,sfclats)
      loninnow = sfclons
      xlons_sfc, dummy = um.addcyclic(xlons_sfc, loninnow)
      xlats_sfc, dummy = um.addcyclic(xlats_sfc, loninnow)
      
      trth = f.variables['pt'][record_numnow,::-1,:].squeeze() 
      trpr = (f.variables['pres'][record_numnow,::-1,:].squeeze()) /100.
      if plot_wind_level >= 1:     
          f2 = netCDF4.Dataset(fpath_sfc,'r') 
          levels = f2.variables['level'][:]
          indnow = np.where(levels==plot_wind_level)
          print(indnow)
          uin = f2.variables['u'][record_numnow,indnow[0],::-1,:].squeeze()
          vin = f2.variables['v'][record_numnow,indnow[0],::-1,:].squeeze()          
      else:
          uin = f.variables['u'][record_numnow,::-1,:].squeeze()
          vin = f.variables['v'][record_numnow,::-1,:].squeeze()
      if plot_field == 'sfcwind':         
         uin_sfc = f2.variables['u10'][record_numnow,:,:].squeeze()
         vin_sfc = f2.variables['v10'][record_numnow,:,:].squeeze()               
      if plot_slp == 'True':
          slp = (f2.variables['msl'][record_numnow,::-1,:].squeeze())/100.     
      if plot_field == 'trth_slp_ice':
          lonsf2 = f2.variables['longitude'][:]
          latsf2 = f2.variables['latitude'][::-1]               
          xlonsf2, xlatsf2 = np.meshgrid(lonsf2, latsf2)  
          if 1 == 0:
              conc = (f2.variables['ci'][record_numnow,::-1,:].squeeze())  
          else:
              aa = f2.variables['siconc'][:].squeeze()
              #aa = f2.variables['ci'][:].squeeze()
              dconc = np.zeros_like(aa).astype('f')
              del aa
              #dconc[1:,:,:] = (f2.variables['ci'][1:,::-1,:]) -  (f2.variables['ci'][0:-1,::-1,:])
              dconc[1:,:,:] = (f2.variables['siconc'][0:-1,::-1,:]) -  (f2.variables['siconc'][1:,::-1,:])
              if record_numnow >= 11:
                  conc = np.nanmax(dconc[record_numnow-11:record_numnow+1,:,:],0)
              else:
                  conc = np.nanmax(dconc[0:record_numnow+1,:,:],0)
      else:
          slp = trth
      trpresin = f.variables['pres'][record_numnow,::-1,:].squeeze()
      trte = wm.theta_to_temp(trth,trpresin)

      
      #xlons = f.variables['longitude'][::-1,:]
      #lons = f.variables['longitude'][:]
      #lats = f.variables['latitude'][::-1] 
      #xlats = f.variables['latitude'][::-1,:] 
      xlons, xlats = np.meshgrid(lons, lats)   
      
      levs = [2.0]             
      
      inds = np.where(trth>400)
      uin[inds] = float('NaN');
      vin[inds] = float('NaN');  
      
      if plot_field == 'pw':
          spechum = (f2.variables['q'][record_numnow,:,:,:].squeeze()) 
      tempin = (f2.variables['t'][record_numnow,:,:,:].squeeze()) 
      hght = (f2.variables['z'][record_numnow,:,:,:].squeeze())/9.81 
      pin = f2.variables['level'][:]*100         
      qv = (spechum/(1-spechum))    

      pres = np.zeros_like(tempin).astype('f')   
      for kk in range(0,len(levs)):      
          pres[kk,:,:] = pin[kk]     
          
      qvtot = wm.total_col(qv[::-1,:,:], pres[::-1,:,:], tempin[::-1,:,:], hght[::-1,:,:])         
      if ( (plot_field == 'totqv') or (plot_field == 'trth_totqv') ):
          qvtot = (f2.variables['tcwv'][record_numnow,::-1,:].squeeze()) 
          f2.close
      if plot_cloud_cover == 'True':
          tcc = (f2.variables['tcc'][record_numnow,:,:].squeeze())
      plot_slp = 'False'
      f2.close      
      
      if plot_slp == 'True':
          f2.close
 
   elif data_gridnum == 20:
      lons = f.variables['longitude'][:]
      lats = f.variables['latitude'][::-1]          
      if( (plot_slp == 'True') or (plot_field == 'trth_totqv') or (plot_field == 'cloud_sw') or (plot_field == 'cloud_lw') or (plot_field == 'cloud_net') or (plot_field == 'ivt') or (plot_field == 'trth_ivt') ):
          f2 = netCDF4.Dataset(fpath_sfc,'r') 
          sfclats = f2.variables['latitude'][::-1]
          sfclons = f2.variables['longitude'][:]
          xlons_sfc,xlats_sfc = np.meshgrid(sfclons,sfclats)
          loninnow = sfclons
          xlons_sfc, dummy = um.addcyclic(xlons_sfc, loninnow)
          xlats_sfc, dummy = um.addcyclic(xlats_sfc, loninnow)
	  
          slp = (f2.variables['msl'][record_numnow,::-1,:].squeeze())/100.
	  
          if ( (plot_field == 'cloud_sw') or (plot_field == 'cloud_lw') or (plot_field == 'cloud_net') ):
              cloud_sw = (f2.variables['ssr'][record_numnow,::-1,:].squeeze())/3600
              cloud_sw_clearsky = (f2.variables['ssrc'][record_numnow,::-1,:].squeeze())/3600
              cloud_lw = (f2.variables['str'][record_numnow,::-1,:].squeeze())/3600
              cloud_lw_clearsky = (f2.variables['strc'][record_numnow,::-1,:].squeeze())/3600
	      
              crf_sw = cloud_sw_clearsky - cloud_sw
              crf_lw = cloud_lw - cloud_lw_clearsky
              crf_net = crf_lw + crf_sw
	      
              xlons = xlons_sfc
              xlats = xlats_sfc
      
          if ( (plot_field == 'ivt') or (plot_field == 'trth_ivt') ):
              ivt = (f2.variables['p72.162'][record_numnow,::-1,:].squeeze())
              xlons = xlons_sfc
              xlats = xlats_sfc	  
	  
      trth = f.variables['pt'][record_numnow,::-1,:].squeeze()
      #trpr = (f.variables['pres'][record_numnow,::-1,:].squeeze()) /100.
      
      if plot_field == 'trwind':
          uin = f.variables['u'][record_numnow,::-1,:].squeeze()
          vin = f.variables['v'][record_numnow,::-1,:].squeeze()
     
          
   elif data_gridnum == -1:
      trth = f.variables['TROP_THETA'][record_numnow,::-1,:].squeeze()
      slp = (f.variables['SLP'][record_numnow,::-1,:].squeeze())/100.     
      trpresin = f.variables['TROP_PRESSURE'][record_numnow,::-1,:].squeeze()
      
      lons = f.variables['lon'][record_numnow,:].squeeze()
      lats = f.variables['lat'][record_numnow,::-1].squeeze() # Read in reverse direction
      levs = []         
   else:
      trtempin = f.variables['TMP_3_PVL_10'][record_numnow,0,::-1,:].squeeze()
      trpresin = f.variables['PRES_3_PVL_10'][record_numnow,0,::-1,:].squeeze()
      try:
         slp = f.variables['MSLET_3_MSL_10'][record_numnow,::-1,:].squeeze()/10**2   
      except:
         slp = f.variables['PRMSL_3_MSL_10'][record_numnow,::-1,:].squeeze()/10**2   
      lons = f.variables['lon_3'][:]
      lats = f.variables['lat_3'][::-1] # Read in reverse direction
      levs = f.variables['lv_ISBL3'][:]

   f.close
   
   if data_gridnum > -1:
      if ( (data_gridnum != 11) & (data_gridnum != 12) & (data_gridnum != 20) ):
         trth = wm.temp_to_theta(trtempin, trpresin)

   lonin = lons
   trth, lons = um.addcyclic(trth, lonin)
   
   mstats(trth)
   if plot_field == 'trtemp':
       del trth
       trth, dummy = um.addcyclic(trte, lonin)
   if( (plot_slp == 'True') or (plot_field == 'trth_totqv') or (plot_field == 'ivt') or (plot_field == 'trth_ivt')):
       lonin_sfc = sfclons
       slp, dummy = um.addcyclic(slp, lonin_sfc)
       if ( (plot_field == 'ivt') or (plot_field == 'trth_ivt') ):
           lonin_sfc1 = sfclons
           ivt, dummy = um.addcyclic(ivt,lonin_sfc1)     
       if plot_field == 'tr_totqv':
           lonin_sfc1 = sfclons
           ivt, dummy = um.addcyclic(ivt,lonin_sfc1)  

   if ( (plot_field == 'cloud_sw') or (plot_field == 'cloud_lw') or (plot_field == 'cloud_net') ): 
       lonin_sfc1 = sfclons
       crf_sw, dummy = um.addcyclic(crf_sw, lonin_sfc1)
       lonin_sfc2 = sfclons
       crf_lw, dummy = um.addcyclic(crf_lw, lonin_sfc2) 
       lonin_sfc3 = sfclons
       crf_net, dummy = um.addcyclic(crf_net, lonin_sfc3)           
   if plot_wind == 'True':
      u, dummy = um.addcyclic(uin, lonin)
      v, dummy = um.addcyclic(vin, lonin)
      uvmag = np.sqrt(u**2 + v**2)
      
      uq = np.percentile(uvmag, 75)
      #cflevs_wind =  np.arange(60, 201, 10)
      cflevs_wind =  np.arange(np.round(uq), 201, 5)

   
   if data_gridnum != 11 :
      X, Y = np.meshgrid(lons, lats)      
   else:
      xlons, dummy = um.addcyclic(xlons, lonin)
      xlats, dummy = um.addcyclic(xlats, lonin)
      X = xlons
      Y = xlats

   if map_projection == 'ortho' :
      trth_thresh = 380
      cbar_min_trth = 270   
   
   else:
      #trth_thresh = 500
      #cbar_min_trth = 305
      trth_thresh = 380
      cbar_min_trth = 270
   

   if ( (plot_field == 'trth') or (plot_field == 'trth_totqv') or (plot_field == 'trth_ivt') ) :
       if (plot_field == 'trth'):
           figname_desc = 'trth'
       if (plot_field == 'trth_totqv'):
           figname_desc = 'trth_totqv'
       if (plot_field == 'trth_ivt'):
           figname_desc = 'trth_ivt'
           ivtplot = -1*ivt
           cntr_max_ivt = 1800
           cntr_min_ivt = 200
           cint_ivt = 300
           cflevs_ivt = np.arange(cntr_min_ivt, cntr_max_ivt+(cint_ivt/2), cint_ivt) 
       if convert_to_celsius == 'True':
           cbarlabel = r'$^{\circ}$C'
       else:
           #cbarlabel = 'Potential temperature (K)'
           cbarlabel = 'Kelvin'
   
       #cbar_max_trth = trth_thresh
       #cbar_max_trth = 335
       #cbar_max_trth = 350
       #cbar_max_trth = 400
       
       if convert_to_celsius == 'True':
           cbar_max_trth = np.ceil(cbar_max_trth - 273.15)
           cbar_min_trth = np.floor(cbar_min_trth - 273.15)
           trth_thresh = np.ceil(trth_thresh - 273.15)
       
       if plot_trth_contours == 'True':          
           cbar_max_trth = 380
           cint_trth = 1
       else:
           cbar_min_trth = 275   
           cbar_max_trth = 345
           cint_trth = 1
       
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth+cint_trth, cint_trth)
       #cflevs_trth_cntrs = cflevs_trth[0:12]
       cflevs_trth_cntrs = cflevs_trth
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth+cint_trth,4*cint_trth)
       cmap_opt = plt.cm.jet
       #cmap_opt = plt.cm.Blues_r

       trth = ndimage.gaussian_filter(trth,0.75)
       
       if convert_to_celsius == 'True':
           trth[:] = trth[:] - 273.15    


       trth = um.filter_numeric_nans(trth,trth_thresh,trth_thresh,'high')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       
       if plot_trth_contours == 'False':
           trth = um.filter_numeric_nans(trth,cbar_max_trth,cbar_max_trth,'high')
           trth = um.filter_numeric_nans(trth,cbar_min_trth,cbar_min_trth,'low')
       
       plotvar = trth
       mstats(plotvar)
                  
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)                                
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       if (plot_field == 'trth'):
           titletext1 = 'Tropopause potential temperature %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))             
       if (plot_field == 'trth_totqv'):
           titletext1 = 'Tropopause potential temperature %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))       
       if (plot_field == 'trth_ivt'):
           titletext1 = 'Tropopause potential temperature, SLP, and IVT %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'trth_slp_ice' :
       
       figname_desc = 'trth_slp_ice'   
   
       if convert_to_celsius == 'True':
           cbarlabel = r'$^{\circ}$C'
       else:
           #cbarlabel = 'Potential temperature (K)'
           cbarlabel = 'Fractional decrease'
   
       
       loninf2 = lonsf2
       xlonsf2, dummy = um.addcyclic(xlonsf2, loninf2)
       xlatsf2, dummy = um.addcyclic(xlatsf2, loninf2)
       Xf2 = xlonsf2
       Yf2 = xlatsf2       
       #cbar_min_trth = 0.15
       #cbar_max_trth = 1
       #cint_trth = 0.05
       
       if 1 == 0:
           cbar_min_trth = 0.8
           cbar_max_trth = 1.0        
           cint_trth = 0.01       
       else:
           cbar_min_trth = 0.02
           cbar_max_trth = 0.2        
           cint_trth = 0.01
       
       
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth+cint_trth, cint_trth)       
       #cflevs_trth_cntrs = cflevs_trth
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth+cint_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r
       #cmap_opt = plt.cm.Blues
       #cmap_opt = plt.cm.seismic_r

       trth = ndimage.gaussian_filter(trth,0.75)
       
       if convert_to_celsius == 'True':
           trth[:] = trth[:] - 273.15    


       trth = um.filter_numeric_nans(trth,trth_thresh+cint_trth,trth_thresh+cint_trth,'high')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       
       
       conc, dummy = um.addcyclic(conc, lonin)
       plotvar = conc
       mstats(plotvar)
       
       #plotvar = um.filter_numeric_nans(plotvar,cbar_min_trth,float('NaN'),'low')        
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)                                
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       cflevs_trth_cntrs =  np.arange(250, 311, 5)       
       
       titletext1 = 'Tropopause potential temperature %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))        
   if plot_field == 'trpr' :
       
       figname_desc = 'trpr' 
       cbarlabel = 'Pressure (hPa)'    
       
       cint_trth = 20
       cbar_max_trth = 800.
       cbar_min_trth = 300.
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       #cflevs_trth_cntrs = cflevs_trth[0:12]
       cflevs_trth_cntrs = cflevs_trth
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r
       #cmap_opt = plt.cm.Blues_r

       trpr = ndimage.gaussian_filter(trpr,0.75)
       
 


       trpr = um.filter_numeric_nans(trpr,cbar_max_trth+cint_trth,cbar_max_trth+cint_trth,'high')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       
       trpr, dummy = um.addcyclic(trpr, lonin)
       plotvar = trpr
       mstats(plotvar)
                  
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)                                
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       titletext1 = 'Tropopause pressure %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))       
   if plot_field == 'trtemp' :
       figname_desc = 'trtemp'
       trth_thresh = 220
       cbar_min_trth = 220       
       if convert_to_celsius == 'True':
           cbarlabel = r'$^{\circ}$C'
       else:
           #cbarlabel = 'Potential temperature (K)'
           cbarlabel = 'Kelvin'
   
       #cbar_max_trth = trth_thresh
       #cbar_max_trth = 335
       cbar_max_trth = 250
       #cbar_max_trth = 400
       
       if convert_to_celsius == 'True':
           cbar_max_trth = np.ceil(cbar_max_trth - 273.15)
       cbar_min_trth = np.floor(cbar_min_trth - 273.15)
       trth_thresh = np.ceil(trth_thresh - 273.15)
       
       cint_trth = 2
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[-10:]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       #cmap_opt = plt.cm.jet
       cmap_opt = plt.cm.Blues


       
       if convert_to_celsius == 'True':
           trth[:] = trth[:] - 273.15    


       trth = um.filter_numeric_nans(trth,trth_thresh-cint_trth,float('NaN'),'low')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       

       plotvar = trth
       mstats(plotvar)

                  
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)                                
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       titletext1 = 'Tropopause temperature gradient %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))        
   if ( (plot_field == 'pw') or (plot_field == 'totqv') or (plot_field == 'trth_totqv') ):
       if (plot_field == 'pw'):
           figname_desc = 'pw'
       if (plot_field == 'totqv'):
           figname_desc = 'totqv'
       if (plot_field == 'trth_totqv'):
           figname_desc = 'trth_totqv'       
       qvtot, dummy = um.addcyclic(qvtot, lonin)   
       
       if plot_field == 'trth_totqv':
           plotvar = trth
           cbarlabel = 'Kelvin'

           plot_trth_contours ='False'
           cbar_max_trth_overlay = 50.
           cbar_min_trth_overlay = 10.0
           cint_trth_overlay = 0.5


           cflevs_trth_overlay = np.logspace(np.log10(cbar_min_trth_overlay), np.log10(cbar_max_trth_overlay), num=10, endpoint=True, base=10.0)
           cflevs_trth_cntrs_overlay = cflevs_trth_overlay[::2]       
       
       else:       
           plotvar = qvtot
           cbarlabel = 'Precipitable water (mm)'
           cbar_max_trth = 20.
           cbar_min_trth = 0.1
           cint_trth = 0.5


       cflevs_trth = np.logspace(np.log10(cbar_min_trth), np.log10(cbar_max_trth), num=50, endpoint=True, base=10.0)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)

       #cflevs_trth = np.logspace(0.4, 1.2, num=100, endpoint=True, base=10.0) # Range is from 10^-1 to 10^9 meters with 100 points
       #cflevs_trth_cntrs = cflevs_trth[::2]
       #cflevs_trth_ticks = cflevs_trth[::2]

       #cmap_opt = plt.cm.BuGn 
       cmap_opt = plt.cm.Greens
       
       titletext1 = 'Precipitable water vapor %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))        
       
       if plot_trth_contours == 'True':
           cbar_min_trth_overlay = 260
           cbar_max_trth_overlay = 380
           cint_trth_overlay = 10       
           cflevs_trth_overlay =  np.arange(cbar_min_trth_overlay, cbar_max_trth_overlay+cint_trth_overlay, cint_trth_overlay)
       
       
                 
       
   if plot_field == 'trwind' :       
       figname_desc = 'trwind'   
       trwind = np.sqrt(uin**2 + vin**2)
       trwind, dummy = um.addcyclic(trwind, lonin)          
       plotvar = trwind
       cbarlabel = 'Tropopause wind (m s-1)'
   
       cbar_max_trth = 100
       cbar_min_trth = 20
       cint_trth = 2
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.hot_r
       
       titletext1 = 'Tropopause wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))      
   if plot_field == 'quivwind' :
       figname_desc = 'quivwind'
       trwind = np.sqrt(uin**2 + vin**2)       
       trwind, dummy = um.addcyclic(trwind, lonin)   
       uwind, dummy = um.addcyclic(uin, lonin) 
       vwind, dummy = um.addcyclic(vin, lonin)   
       plotvar = trwind
       cbarlabel = 'Tropopause wind (m s-1)'
   
       cbar_max_trth = 50
       cbar_min_trth = 10
       cint_trth = 1
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       #cmap_opt = plt.cm.RdBu_r
       cmap_opt = plt.cm.gist_heat_r
       
       titletext1 = 'Tropopause wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))      

   if plot_field == 'sfcwind' :
       
       figname_desc = 'sfcwind'

       trwind = np.sqrt(uin**2 + vin**2)
       sfcwind = np.sqrt(uin_sfc**2 + vin_sfc**2)
       sfcwind, dummy = um.addcyclic(sfcwind, lonin)   
       trwind, dummy = um.addcyclic(trwind, lonin)
       plotvar = sfcwind
       cbarlabel = '2-meter AGL wind (m s-1)'
   
       cbar_max_trth = 35
       cbar_min_trth = 10
       cint_trth = 1
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r
       
       titletext1 = 'Surface wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'ghgt500' :
       
       figname_desc = 'ghgt500'       
       
       ghgt500, dummy = um.addcyclic(ghgt500, lonin)   
       plotvar = ghgt500
       cbarlabel = '500 hPa heights (m)'
   
       cbar_max_trth = 5910
       cbar_min_trth = 4910
       cint_trth = 30
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.RdBu_r       
       #cmap_opt = plt.cm.jet
       contour2plot = 5420
       
       titletext1 = '500 hPa height %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'trthgrad':       
       figname_desc = 'trthgrad'

       dthdy, dthdx = wm.gradient_sphere(trth, lats, lons)
       plotvar = np.sqrt(dthdx**2 + dthdy**2)*1000
       plotvar = ndimage.gaussian_filter(plotvar,0.75)
       
       cbar_max_trth = 0.05
       cbar_min_trth = 0.01
       cint_trth = 0.001
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r          
       cbarlabel = 'Gradient (K km-1)'
   
       titletext1 = 'Tropopause potential temperature gradient %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 

   if ( (plot_field == 'cloud_sw') or (plot_field == 'cloud_lw') or (plot_field == 'cloud_net') ): 
       if plot_field == 'cloud_sw':
           figname_desc = 'cloud_sw'
           plotvar = crf_sw
           titletext1 = 'Shortwave cloud radiative forcing %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))     
           cbar_max_trth = 200
           cbar_min_trth = 0
           cint_trth = 5       
       if plot_field == 'cloud_lw':
           figname_desc = 'cloud_lw'
           plotvar = crf_lw
           titletext1 = 'Longwave cloud radiative forcing %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))     
           cbar_max_trth = 150
           cbar_min_trth = 6
           cint_trth = 3      
       if plot_field == 'cloud_net':
           figname_desc = 'cloud_net'
           plotvar = crf_net
           titletext1 = 'Net cloud radiative forcing %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))     
           cbar_max_trth = 600
           cbar_min_trth = 10
           cint_trth = 10    
       
       mstats(plotvar)


       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r #plt.cm.RdBu_r               
       cbarlabel = 'Cloud radiative focing (W m-2)'
   
   if plot_field == 'ivt':
       figname_desc = 'ivt'
       plotvar = -1*ivt # negative is southward.
       #plotvar = ivt
       titletext1 = 'Vertical integral of poleward water vapor flux %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))     
       cbar_max_trth = 1200
       cbar_min_trth = 100
       cint_trth = 10              
       
       mstats(plotvar)


       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r #plt.cm.RdBu_r               
       cbarlabel = 'IVT (kg m-1 s-1)'        

      
   
   base_cntr_slp = 996
   nslpconts = 20
   cint_slp = 4
   cbar_min_slp = base_cntr_slp-nslpconts*cint_slp
   cbar_max_slp = base_cntr_slp+nslpconts*cint_slp
   #cbar_max_slp = cbar_max_slp + (cint_slp/2)   
   cbar_max_slp = 992 #1008
   cflevs_slp =  np.arange(cbar_min_slp, cbar_max_slp+cint_slp, cint_slp)
   cflevs_slp_ticks = np.arange(cbar_min_slp,cbar_max_slp+cint_slp,4*cint_slp)

   if smooth_plotvar == 'True':
       plotvar = ndimage.gaussian_filter(plotvar,0.75) 

   
   
   # Calculate distance from Brunswick, Maine
   lat1 = latlon_radius[0]
   lon1 = latlon_radius[1]
   ed = um.earth_distm(lat1,lon1,Y,X)
      
   

 
   #titletext1 = 'Tropopause potential temperature valid %s at %s UTC' % (
   #        dt.strftime('%d %b %Y'), dt.strftime('%H00'))  
   #titletext1 = 'Tropopause potential temperature and SLP valid %s' % (dt.strftime('%Y%m%d%H'))     
   #titletext1 = 'Tropopause potential temperature %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   
   # Set global figure properties
   golden = (np.sqrt(5)+1.)/2.   

   # Figure 1
   
   fig = plt.figure(figsize=(12.,12.), dpi=128)   # New figure
   ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
   if map_projection == 'ortho':
    if zoom == 'False':   
        m = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
              resolution = 'l', area_thresh = 1000.,ax=ax1)
    else:
        m1 = Basemap(projection='ortho',lon_0=proj_latlon[1],lat_0=proj_latlon[0],resolution=None)

        width = m1.urcrnrx - m1.llcrnrx
        height = m1.urcrnry - m1.llcrnry
    
        width = width*coef
        height = height*coef
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        m = Basemap(projection='ortho',lon_0=proj_latlon[1],lat_0=proj_latlon[0],resolution='l',llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)      
          
   elif map_projection == 'lcc':
   #    m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-40.,urcrnrlat=50.352,\
   #    m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
       m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
               rsphere=(6378137.00,6356752.3142),\
               resolution='l',area_thresh=1000.,projection='lcc',\
               lat_1=50.,lon_0=-107.,ax=ax1)           
   elif map_projection == 'npstere':
       if zoom == 'false':
          m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
       else:
          m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
   elif map_projection == 'spstere':
       if zoom == 'false':
          m = Basemap(projection='spstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
       else:
          m = Basemap(projection='spstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')

   # draw countries, states, and differentiate land from water areas.
   m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
   m.drawcountries(linewidth=1, color='#444444', zorder=5)
   m.drawstates(linewidth=0.66, color='#444444', zorder=4)
   m.drawmapboundary
   if plot_field == 'trth_slp_ice':
       m.fillcontinents(color='Wheat',lake_color='lightblue', zorder=1) 
   # draw lat/lon grid lines every 30 degrees.
   #m.drawmeridians(np.arange(0, 360, 30))
   m.drawparallels(np.arange(-90, 90, 10))
   m.drawmeridians(np.arange(0, 360, 30),labels=[True,True,True,True])
   
   if plot_field == 'blank' :
   
       if map_projection == 'lcc':
           m.drawmeridians(np.arange(0, 360, 15),labels=[False,True,True,False])
           m.drawparallels(np.arange(-90, 90, 30),labels=[False,True,True,False])          
       
       m.fillcontinents(color='tan',lake_color='lightskyblue')
       m.drawmapboundary(fill_color='lightskyblue')
       m.drawparallels([60,60],labels=[1,1,1,1],latmax=60,linewidth=6,color='b')
       save_name = imagedir + 'blank_map' + '_' + date_string + ".png"
       plt.savefig(save_name, bbox_inches='tight')
       exit()
   x, y = m(X, Y)
   if plot_field == 'trth_slp_ice' :
       xf2, yf2 = m(Xf2,Yf2)

   #trth = um.filter_numeric_nans(trth,340,float('NaN') ,'high')
   
   # Contour tropopause potential temperature
   if plot_field == 'trth_slp_ice':
       CS1 = m.contourf(xf2,yf2,plotvar,cmap=cmap_opt,levels=cflevs_trth, extend='both',zorder=1)
   else:
       CS1 = m.contourf(x,y,plotvar,cmap=cmap_opt,levels=cflevs_trth, extend='both',zorder=1)
   if plot_field == 'trtemp' :
       CS1.cmap.set_under('w')
   if plot_field == 'trth' :
       CS1.cmap.set_over('w')
   cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
   
   #clabs = ['%i' % f for f in cflevs_trth]
   
   labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
   cbar.ax.set_xticklabels(labels, size=20)
   cbar.set_label(cbarlabel,size=20)

   if plot_field == 'quivwind' :
       ugrid,newlons = shiftgrid(180.,uwind,lons,start=False)
       vgrid,newlons = shiftgrid(180.,vwind,lons,start=False)
       # transform vectors to projection grid.
       
       uproj,vproj,xx,yy = m.transform_vector(ugrid[::-1,:],vgrid[::-1,:],newlons,lats[::-1],61,61,returnxy=True,masked=True)
       # now plot.
       Q = m.quiver(xx,yy,uproj,vproj,scale=700,color='g',linewidths=(1,), edgecolors=('g'), headaxislength=2)
  
   if plot_spc_reports == 'true':
       data = np.recfromtxt(fpath_spc, unpack=True, dtype=None, names=True, delimiter=',')
       date = data['date']
       time = data['time']
       #rating = data['f']
       rating = data['mag']
       slat = data['slat']
       slon = data['slon']
       elon = data['elon']
       elat = data['elat']

       for i in range(date.shape[0]):
           if 1 == 0 :
               print(date[i], time[i])


       dts = []
       print(date.shape)
       for i in range(date.shape[0]):
           try:
               datetime_object = datetime.datetime.strptime(date[i] + ' ' + time[i], '%m/%d/%y %H:%M:%S')
           except:
               datetime_object = datetime.datetime.strptime(date[i] + ' ' + time[i], '%Y-%m-%d %H:%M:%S')
       dts.append(datetime_object)

       print(dts[0])
       print(dts[100])
       print(dts[1])
       # Check to see if a date/time is earlier than another date:
       print(dts[0] < dts[100])

       # Check to see if two dates are equal:
       print(dts[0] == dts[1])

       # What's the difference between two dates:
       print(dts[100] - dts[0])


       begin_date = datetime.datetime(2010, 5, 10, 12, 0, 0)    # format here is "year, month, date, hour, minute, second"
       end_date = datetime.datetime(2010, 5, 11, 12, 0, 0)

       valid_slats = []
       valid_slons = []
       valid_elats = []
       valid_elons = []
       valid_ratings = []
       valid_dates = []

       # Now we loop through our original data to pull out the points that are valid
       torcount = 0
       for i in range(len(dts)):
           if dts[i] >= begin_date and dts[i] <= end_date:
               if rating[i] >= 0: # Filter for only violent tornadoes
                   valid_dates.append(dts[i])
                   valid_slats.append(slat[i])
                   valid_slons.append(slon[i])
                   valid_elats.append(elat[i])
                   valid_elons.append(elon[i])
                   valid_ratings.append(rating[i])   
                   torcount += 1
       print("There were %3d tornado reports meeting the specified criteria" %(torcount))
       for i in range(len(valid_slats)):
           if valid_elons[i] == 0 or valid_elats[i] == 0:
               valid_elons[i] = valid_slons[i]
               valid_elats[i] = valid_slats[i]

           new_slons, new_slats = m(valid_slons, valid_slats)
           new_elons, new_elats = m(valid_elons, valid_elats)

       for i in range(len(valid_dates)):
           x1 = new_slons[i]
           y1 = new_slats[i]
           x2 = new_elons[i]
           y2 = new_elats[i]
           #m.plot([x1, x2], [y1, y2], 'rv', linestyle='solid', linewidth=3,markersize=16)   
       #m.plot([x1, x2], [y1, y2],linestyle='None', marker='v', color='0.5', linewidth=3,markersize=20)
           m.plot([x1, x2], [y1, y2],'rv',linestyle='solid',linewidth=3,markersize=30)
     
   if plot_trth_contours == 'True':
       if plot_field == 'trth_slp_ice':
           CS2 = m.contour(x, y, trth, cflevs_trth_cntrs, colors='b', linewidths=3)
       elif ( (plot_field == 'totqv') or (plot_field == 'pw') ):
           CS2 = m.contour(x, y, trth, cflevs_trth_overlay, linestyles='solid',colors='0.5', linewidths=1.5)
           plt.clabel(CS2, cflevs_trth_overlay, fmt = '%i', inline=True, fontsize=10) 
       else:           
           CS2 = m.contour(x, y, trth, cflevs_trth_cntrs, colors='k', linewidths=0.25)
       
   if plot_field == 'trth_totqv':
           #CS2 = m.contour(x, y, qvtot, cflevs_trth_overlay, linestyles='solid',colors='0.5', linewidths=1.0)
       #plt.clabel(CS2, cflevs_trth_overlay, fmt = '%i', inline=True, fontsize=10) 
       
       CS2 = m.contour(x, y, qvtot, [12.,14.], linestyles='solid',colors='0.5', linewidths=4.0)
       plt.clabel(CS2, [12,14], fmt = '%i', inline=True, fontsize=10) 
       
              
   if plot_slp == 'True':      
      #xslp, yslp = m(X, Y)
      xslp, yslp = m(xlons_sfc,xlats_sfc)
      CS3 = m.contour(xslp, yslp, slp, cflevs_slp, colors='k', linewidths=1.5)
      #plt.clabel(CS3, inline=1, fontsize=10, fmt='%i')   
      plt.clabel(CS3, cflevs_slp, fmt = '%i', inline=True, fontsize=10)

   if plot_field == 'trth_ivt':
      xivt, yivt = m(xlons_sfc,xlats_sfc)
      #CS4 = m.contour(xivt, yivt, ivtplot, cflevs_ivt, colors='violet', linewidths=3)  
      #plt.clabel(CS4, cflevs_ivt, fmt = '%i', inline=True, fontsize=10)
      CS4a = m.contour(xivt, yivt, ivtplot, cflevs_ivt, colors='black', linewidths=6)  
      CS4b = m.contour(xivt, yivt, ivtplot, cflevs_ivt, colors='white', linewidths=4)  

   if plot_cloud_cover == 'True':
      tcc, lons = um.addcyclic(tcc, lonin)
      #tcc = ndimage.gaussian_filter(tcc,0.5)      
      
      
      xx,yy=np.where(tcc >= 0.1)
      sig=np.copy(tcc)
      sig[:,:]=1.0 
      sig[xx,yy]=0.0
      
      CS3 = m.contour(x, y, tcc, levels=[0.1,0.5], colors='0.5', linewidths=1.5)
      #CS3 = m.contourf(xcloud, ycloud, tcc, levels=[0.1,1.1], cmap=plt.cm.Greys_r,alpha=0.1)
      
      #CS3 = m.contourf(xcloud,ycloud,sig,levels=[0.1,1.0],hatches=["//"],alpha=0)

      cflevs_cloud = [np.nanmin(tcc), 0.9, np.nanmax(tcc)]     
      #CS3 = m.contourf(xcloud,ycloud,tcc,levels=cflevs_cloud,hatches=["", "//"], alpha=0) 
      
      #plt.clabel(CS3, inline=1, fontsize=10, fmt='%i')   
      #plt.clabel(CS3, cflevs_cloud, fmt = '%i', inline=True, fontsize=10)
   if plot_wind == 'True':
      #CSW = m.contour(x, y, uvmag, cflevs_wind, colors='Mediumorchid', linewidths=2.5)
      CSW = m.contour(x, y, uvmag, cflevs_wind, colors='0.7', linewidths=2.5)
      #plt.clabel(CSW, inline=1, fontsize=10, fmt='%i')   
      plt.clabel(CSW, cflevs_wind, fmt = '%i', inline=True, fontsize=10)   

      
   if plot_field == 'ghgt500':
      CS3 = m.contour(x, y, plotvar, [contour2plot, contour2plot], colors='green', linewidths=6.0)

   if plot_cross_section_line == 'true':             
       xx, yy = um.xsection_inds(lon_range[0], lat_range[0], lon_range[1], lat_range[1], X, Y, m)
       ax1.plot(x[xx,yy], y[xx,yy], color='0.75', lw=8)
     
   if plot_radius == 'True':
      # Calculate distance from Brunswick, Maine
      lat1 = latlon_radius[0]
      lon1 = latlon_radius[1]
      ed = um.earth_distm(lat1,lon1,Y,X)
      #loiter
      rcruise = 4900 # nm
      vcruise = 510  # knots
      loiter_times = np.array([6.0,4.0,2.0])
      cflevs_loiter = (['6-h','4-h','2-h'])
      dbase = (rcruise/2) - ((loiter_times*vcruise)/(2*1.16))

      #range_km = (dbase*1.15)/0.621371
      range_km = dbase*1.852
      print(range_km)
      
      
      strs = [ '6-h', '4-h', '2-h']

      
      CS4 = m.contour(x,y,ed,range_km,colors='m',linewidths=3.0)
      fmt = {}
      for l,s in zip( CS4.levels, strs ):
         fmt[l] = s      
      plt.clabel(CS4,CS4.levels,inline=True,fmt=fmt,fontsize=16) 
      
   if plot_trth_contour == 'True':            
      #trthcnt = ndimage.gaussian_filter(trth,0.5)
      trthcnt = trth
      if plot_field == 'trth':
          CS4 = m.contour(x,y,trthcnt,[305.,310.,315.],colors='k',linewidths=5.0)   
      if plot_field == 'trtemp':
            CS4 = m.contour(x,y,trthcnt,[230.,240.,250.],colors='k',linewidths=5.0)     
   if plot_radiosondes == 'True':
       xpt,ypt = m(rad_lons,rad_lats)
       if zoom == 'false':
           CSx = m.scatter(xpt,ypt,c='m',s=10,zorder=10,marker='+')   
       else:
           CSx = m.scatter(xpt,ypt,c='m',s=20,zorder=10,marker='+')   
   #ax1.set_title(titletext1)
   #plt.title('%s' % (titletext1), fontsize=label_fontsize,bbox=dict(facecolor='white', alpha=0.65),x=0.5,y=.95,weight = 'demibold',style='oblique', \
   #        stretch='normal', family='sans-serif')

   if loop_option == 'false':      
      #plt.show()
      aa = 1
   #for ii in range(0,len(russia_fir_lons)-1):       
   #    m.drawgreatcircle(russia_fir_lons[ii],russia_fir_lats[ii],russia_fir_lons[ii+1],russia_fir_lats[ii+1],linewidth=6,color='0.3',linestyle=':')
       
   #for ii in range(0,len(mosaic_lons)-1):       
   #    m.drawgreatcircle(mosaic_lons[ii],mosaic_lats[ii],mosaic_lons[ii+1],mosaic_lats[ii+1],linewidth=6,color='g',linestyle='--')  
   if 1 == 0:
       dist_mission2 = 0.
       for ii in range(0,len(mission2_lons)-1):       
           m.drawgreatcircle(mission2_lons[ii],mission2_lats[ii],mission2_lons[ii+1],mission2_lats[ii+1],linewidth=6,color='k')
           distnow = um.earth_distm(mission2_lats[ii],mission2_lons[ii],mission2_lats[ii+1],mission2_lons[ii+1]) 
           dist_mission2 = dist_mission2 + distnow
       xm,ym = m(mission2_lons,mission2_lats)
       CSM = m.scatter(xm,ym,c='m',edgecolors='k',linewidths=3,s=120,zorder=10,marker='o')       
       print("The total distance is mission 2 is %8.2f km " %(dist_mission2))    
   if 1 == 0:
       dist_mission1 = 0.
       for ii in range(0,len(mission1_lons)-1):       
           m.drawgreatcircle(mission1_lons[ii],mission1_lats[ii],mission1_lons[ii+1],mission1_lats[ii+1],linewidth=6,color='k')
           distnow = um.earth_distm(mission1_lats[ii],mission1_lons[ii],mission1_lats[ii+1],mission1_lons[ii+1]) 
           dist_mission1 = dist_mission1 + distnow
       xm,ym = m(mission1_lons,mission1_lats)
       CSM = m.scatter(xm,ym,c='m',edgecolors='k',linewidths=3,s=120,zorder=10,marker='o')       
       print("The total distance is mission 1 is %8.2f km " %(dist_mission1)) 
   if 1 == 0:
       dist_mission3 = 0.
       for ii in range(0,len(mission3_lons)-1):       
           m.drawgreatcircle(mission3_lons[ii],mission3_lats[ii],mission3_lons[ii+1],mission3_lats[ii+1],linewidth=6,color='k')
           distnow = um.earth_distm(mission3_lats[ii],mission3_lons[ii],mission3_lats[ii+1],mission3_lons[ii+1]) 
           dist_mission3 = dist_mission3 + distnow
       xm,ym = m(mission3_lons,mission3_lats)
       CSM = m.scatter(xm,ym,c='m',edgecolors='k',linewidths=3,s=120,zorder=10,marker='o')
       print("The total distance is mission 3 is %8.2f km " %(dist_mission3))       
   
   save_name = imagedir + figname + '_' + figname_desc + '_' + date_string + ".png"
   plt.savefig(save_name, bbox_inches='tight')   
   record_start += 1
   rcount +=1
   record_numnow +=1

   datenow = um.advance_time(datenow,hinc)      

if loop_option == 'false':
    plt.show()
