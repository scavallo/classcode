#!/usr/bin/python

import numpy as np
#import matplotlib.pyplot as mpl
import matplotlib as plt
from pyproj import Geod
#from mpl_toolkits.basemap import Basemap, addcyclic 
import math

from mstats import *

def nan2zero(data):
    ''' Convert NaNs to zero '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    #temp = np.reshape(data,np.prod(np.size(data)), 1) 
    # SMC: New August 2021
    temp = np.reshape(data, (np.prod(np.size(data)),1));          
    
    # Find indices with NaNs
    inds = np.argwhere(np.isnan(temp))    
    
    # Replace NaNs with zero
    temp[inds] = 0.                 
    
    # Turn vector back into array
    #data = np.reshape(temp,dimens,order='F').copy()
    # SMC new August 2021
    data = np.reshape(temp,(dimens))  
 
    return data

def nan2replace(data,replace_value):
    ''' Convert NaNs to zero '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    #temp = np.reshape(data,np.prod(np.size(data)), 1)  
    # SMC: New August 2021
    temp = np.reshape(data, (np.prod(np.size(data)),1));             
    
    # Find indices with NaNs
    inds = np.argwhere(np.isnan(temp))    
    
    # Replace NaNs with zero
    temp[inds] = replace_value                 
    
    # Turn vector back into array
    #data = np.reshape(temp,dimens,order='F').copy()
    # SMC new August 2021
    data = np.reshape(temp,(dimens))     
 
    return data

def zero2nan(data):
    ''' Convert zeros to Nans '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    #temp = np.reshape(data,np.prod(np.size(data)), 1)  
    # SMC: New August 2021
    temp = np.reshape(data, (np.prod(np.size(data)),1));     
    
    # Find indices with NaNs
    inds = np.argwhere(temp==0)    
    
    # Replace zeros with NaNs
    temp[inds] = float('NaN')                 
    
    # Turn vector back into array
    #data = np.reshape(temp,dimens,order='F').copy()
    # SMC new August 2021
    data = np.reshape(temp,(dimens)) 
 
    return data

def FixNaNs(arr):
    ''' Convert NaNs to nearest neighbor value '''
    ''' '''
    ''' arr:  Input data array '''
    ''' data: Output data array '''
    dimens = np.shape(arr)
    
    # Temporarily collapse data array
    #temp = np.reshape(arr,np.prod(np.size(arr)), 1) 
    # SMC: New August 2021
    temp = np.reshape(arr, (np.prod(np.size(arr)),1));
              
    idxs=np.nonzero(temp==temp)[0]

    if len(idxs)==0:
        return None

    ret=temp

    for i in range(0, idxs[0]):
        ret[i]=ret[idxs[0]]

    for i in range(idxs[-1]+1, ret.size):
        ret[i]=ret[idxs[-1]]

    # Turn vector back into array
    #data = np.reshape(ret,dimens,order='F').copy()
    # SMC new August 2021
    data = np.reshape(ret,(dimens))    
    
    return data

def filter_numeric_nans(data,thresh,repl_val,high_or_low) :
    ''' Filter numerical nans above or below a specified value'''
    ''' '''
    ''' data:        (Input) array to filter '''
    ''' thresh:      (Input) threshold value to filter above or below '''
    ''' repl_val:    (Input) replacement value'''
    ''' high_or_low: (Input)''' 


    dimens = np.shape(data)    
    #temp = np.reshape(data,np.prod(np.size(data)), 1)     
    # SMC: New August 2021
    temp = np.reshape(data, (np.prod(np.size(data)),1));
    if high_or_low=='high':        	
	    inds = np.argwhere(temp>thresh) 	
	    temp[inds] = repl_val	 
    elif high_or_low=='low':    
        inds = np.argwhere(temp<thresh) 
        temp[inds] = repl_val	  
    elif high_or_low =='both':
       	inds = np.argwhere(temp>thresh) 	
        temp[inds] = repl_val
        del inds
       	inds = np.argwhere(temp<-thresh) 	
        #temp[inds] = -repl_val	 
        temp[inds[1:,0]] = -repl_val                
    else:
        inds = np.argwhere(temp>thresh) 
        temp[inds] = repl_val	  
    
    # Turn vector back into array
    #data = np.reshape(temp,dimens,order='F').copy()
    # SMC new August 2021
    data = np.reshape(temp,(dimens))
    
 
    return data    
def label_options(ax,fontsize=None,xaxis_opt=True,yaxis_opt=True,bold_opt=True):
    if fontsize is None:
        fontsize = 14
    if xaxis_opt == True:   
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            if bold_opt == True:
                tick.label1.set_fontweight('bold')
    if yaxis_opt == True:      
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            if bold_opt == 'True':
                tick.label1.set_fontweight('bold')
            
def bold_labels(ax,fontsize=None):
    if fontsize is None:
        fontsize = 14
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

#def draw_map_background(m, ax=mpl.gca()):
#    ''' Setup the map background '''
#    m.drawcoastlines(ax=ax, linewidth=2, color='#444444', zorder=6)
#    m.drawcountries(ax=ax, linewidth=1, color='#444444', zorder=5)
#    m.drawstates(ax=ax, linewidth=0.66, color='#444444', zorder=4)
#    m.drawmapboundary
    
def lonswap(d,subtract=0.):
    sh = np.shape(d)	
    midl = sh[1]/2	
    midl = np.round(midl)
    h=d[:,midl:].copy()-subtract
    d[:,midl:]=d[:,:midl].copy()
    d[:,:midl]=h
    return d
	    
def periodic(d,add=0.):
	return np.append( d, (d[:,0].copy()+add).reshape(-1,1) , 1)    
	
def advance_time(timestrin,timeinc):
    ''' Advances or reverses a time by timeinc'''
    ''' '''
    ''' timestrin: (Input) time string in yyyymmddhh format'''
    ''' timeinc:   (Input) number of hours to increment or decrement'''
    '''             Use a negative sign to decrement '''
    
    import datetime
         
    yyyy = timestrin[0:4]
    mm = timestrin[4:6]
    dd = timestrin[6:8]
    hh = timestrin[8:10]	
	
    date=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
    date += datetime.timedelta(hours=timeinc)
    tt = date.timetuple()
    
    yyyy = str(tt[0])
    mm = str(tt[1])
    dd = str(tt[2])
    hh = str(tt[3])
    
    if tt[0]<1000: yy = '0'+mm 
    if tt[1]<10: mm = '0'+mm 
    if tt[2]<10: dd = '0'+dd
    if tt[3]<10: hh = '0'+hh
    
    timestrin = yyyy+mm+dd+hh        
    
    return timestrin
def get_cmap_cust():
    ''' Setup a custom colortable. '''
    cdict = {'red': ((0.00, 240/255., 220/255.),
                         (0.25, 40/255., 20/255.),
                         (0.50, 225/255., 255/255.),
                         (0.75, 150/255., 150/255.),
                         (1.00, 255/255., 255/255.)),

             'green': ((0.00, 240/255., 220/255.),
                         (0.25, 0/255., 50/255.),
                         (0.50, 255/255., 255/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.)),

             'blue': ((0.00, 255/255., 255/255.),
                         (0.25, 160/255., 150/255.),
                         (0.50, 255/255., 170/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.))}
    return plt.colors.LinearSegmentedColormap('cool2warm', cdict, 256)
def cmap_discretize(cmap, N):
   """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)
   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_whitezero(cmap,N,Nwhites,pos):
   """Whites out middle index of a colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
	Nwhites: Number of whites
	pos: Position for white bar; if -1, will place in middle
	                             if  0, will place at bottom

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   if ( pos == -1 ):
      mid = np.round(N/2)
      mid = int(mid)
   else:
      mid = pos
   
   nadd = Nwhites - 1
   
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)

   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       colors[mid] = 1.
       isodd = 0
       if ( np.mod(N,2) == 0 ) :           
           colors[mid-1] = 1.
           isodd = 1
       kk=mid-nadd
       kend=mid+nadd 
       if ( kk < kend ) :      
          while (kk <= kend) :
             colors[kk] = 1.
             kk += 1
       if (isodd == 1 ): colors[kk] = 1.   	  
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors       
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)

def earth_distm(lat1,lon1,lat2,lon2):
    """

   Calculates the distances between a point and all other points given a latitude and longitude

   Input:    
       lat1, lon1 - Coordinate pair to calculate distance from (must be single values)
       lat2, lon2 - Coordinates to calculate distance to (can be vector or array)
   Output:
       ed - distance between pairs in km
       
    Steven Cavallo
    February 2013
    University of Oklahoma
    
    """
    
    latshape = np.shape(lat2)
    latsize = np.size(lat2)
    ndims = len(latshape)
    
    
    R_earth = 6371200
    R_earth = R_earth / 1000
    pid = np.pi/180

    if ndims == 1:
        if lon2 < 0:
            lon2 = lon2 + 360
        if lon1 < 0:
            lon1 = lon1 + 360
	    
        X = lon1
        Y = lat1	     
    else:

        if latsize > 1:
            [iy,ix] = np.shape(lat2)	   

            X = np.zeros_like(lat2).astype('f')   
            Y = np.zeros_like(lon2).astype('f')   
            for ii in range(0,ix):
                for jj in range(0,iy):   
                    X[jj,ii] = lon1
                    Y[jj,ii] = lat1      
        else:
            X = lon1
            Y = lat1   
    	            

    # calculate distance
    ed = R_earth * np.arccos( np.sin(Y*pid) * np.sin(lat2*pid) + np.cos(Y*pid) * np.cos(lat2*pid) * np.cos((lon2 - X)*pid));

    return ed

def findRotAngle(clon,clat,rlon,rlat):
    '''
    
    Find azimuth angle between two lat-lon points (between 0 and 360 degrees)
    Returns angle measured clockwise from due north
    E.g., a due north azimuth from the center point to the rotated point returns 0
          a due east azimuth from the center point to the rotated point returns 90
    
    '''
    g = Geod(ellps='WGS84')
    az12,az21,dist = g.inv(clon,clat,rlon,rlat)
    
    if az12 <= 0:
        az12 = 360.+az12
    return(az12,dist)

def destagger_hor_wind(ustag,vstag):
   """
   u,v = destagger_hor_wind(ustag,vstag)

   destaggers horizontal wind

   Steven Cavallo
   March 2013

   """
   
   nshape = np.shape(ustag)
   nel = len(nshape)   
   if nel == 2:
      Ny,Nx = np.shape(ustag)         
      
      u = np.zeros((Ny,Nx-1))
      v = np.zeros((Ny,Nx-1))
      
      for jj in range(0,Ny):
         for ii in range(0,Nx-1):	 	     
             u[jj,ii] = (ustag[jj,ii+1] + ustag[jj,ii])/2


      for jj in range(0,Ny):
         for ii in range(0,Nx-1):			
             v[jj,ii] = (vstag[jj+1,ii] + vstag[jj,ii])/2   
   
   else:
      Nz,Ny,Nx = np.shape(ustag)   
       
      u = np.zeros((Nz,Ny,Nx-1))
      v = np.zeros((Nz,Ny,Nx-1))
      for kk in range(0,Nz):
          for jj in range(0,Ny):
              for ii in range(0,Nx-1):	 		
                  u[kk,jj,ii] = (ustag[kk,jj,ii+1] + ustag[kk,jj,ii])/2

      for kk in range(0,Nz):
          for jj in range(0,Ny):
              for ii in range(0,Nx-1):			
                  v[kk,jj,ii] = (vstag[kk,jj+1,ii] + vstag[kk,jj,ii])/2   

   return u,v
def destagger_vertical(varin):
   """
   varout = destagger_hor_wind(varin)

   destaggers in the vertical

   Steven Cavallo
   June 2014

   """
   
   nshape = np.shape(varin)
   nel = len(nshape)   
   
   if nel == 4:
      Nt,Nz,Ny,Nx = np.shape(varin)         
      
      varout = np.zeros((Nt,Nz-1,Ny,Nx))      
      
      for kk in range(0,Nz-1):         	 	     
          varout[:,kk,:,:] = (varin[:,kk+1,:,:] + varin[:,kk,:,:])/2
   
   else:
      Nz,Ny,Nx = np.shape(varin)   
       
      varout = np.zeros((Nz-1,Ny,Nx))
      for kk in range(0,Nz-1):	             	 		
          varout[kk,:,:] = (varin[kk+1,:,:] + varin[kk,:,:])/2



   return varout 
def grid_to_true_wind(lon,ug,vg,truelat1,truelat2,stdlon,proj_type):
    """

   converts from grid relative wind to true direction wind.  Based on FSL 
   mapping module.

   Input:    
       lon - 2D longitude array on grid wind points
       ug, vg - 2D arrays of grid u and v winds
       truelat1,truelat2 - true latitudes
       stdlon - standard longitude
       proj_type - projection type
   Output:
       ut,vt - output u and v true winds
       
    Steven Cavallo
    March 2013
    University of Oklahoma
    
    """
        
    
    Ny,Nx = np.shape(ug)         
      
    ut = np.zeros((Ny,Nx))
    vt = np.zeros((Ny,Nx))    
    
    pid = np.pi/180
    
    if proj_type == 'lambert':
       cone = (np.log(np.cos(truelat1*pid))-np.log(np.cos(truelat2*pid))) / (np.log(np.tan((90. - abs(truelat1)) * pid * 0.5 )) - \
            np.log(np.tan((90. - abs(truelat2)) * pid * 0.5 )) )           
    if proj_type == 'polar':
       cone = 1
       
    if ( (proj_type == 'polar') or (proj_type == 'lambert') ):
       diff = lon - stdlon
       Ny,Nx = np.shape(diff)       
       for jj in range(0,Ny):
          for ii in range(0,Nx):	 	     
              diffnow = lon[jj,ii] - stdlon
              if diffnow > 180:
                  diffnow = diffnow - 360
              if diffnow < -180:
                  diffnow = diffnow + 360                 
       
              alpha = diffnow * cone * pid * 1 * np.sign(truelat1); 
              ut[jj,ii] = vg[jj,ii] * np.sin(alpha) + ug[jj,ii] * np.cos(alpha);
              vt[jj,ii] = vg[jj,ii] * np.cos(alpha) - ug[jj,ii] * np.sin(alpha);    
    
    else:
        ut = ug
        vt = vg

    return ut,vt

def autocorr(datain,endlag):
    '''
    autocorr(datain,endlag)
    
    Input: 
         datain[0:N] is a data time series of size N
	     endlag is the number of time steps to find autocorrelation
    Output:
    	 aut[0:endlag] is the autocorrelation of datain from lag 0 to time step endlag	 
    
    Steven Cavallo
    University of Oklahoma
    July 2016
    '''
    
    N = np.size(datain)
    aut = []
    for lag in range(0,endlag):
        data1 = datain[0:N-lag]
        data1m = data1 - np.nanmean(data1)
        data2 = datain[lag:]
        data2m = data2 - np.nanmean(data2)
        aut.append(np.sum(data1m*data2m)/np.sqrt(np.sum(data1m**2.0)*np.sum(data2m**2.0)))

    return aut

def red_noise_spectrum(datam,conf,chunk_length,rr):
    '''
	red_noise_spectrum(datam,conf,chunk_length,rr)
	
	Determines the 95 and 99 percent confidence level factors for a red noise spectrum
	
	Input:
	    datam: data time series 
	    conf: either 95 or 99 for 95% confidence or 99% confidence interval, respectively
	    chunk_length: Chunk length, in days
	    rr: one-step autocorrelation
	Output:
	    spectout_red: Red noise spectrum
		specout_red_sig: The red noise spectrum's confidence interval (specout_red*specin)
	    fstat: The factors specin is multiplied by    
	    
	Steven Cavallo
	July 2020
    '''
    
    N = len(datam)
    dof = (2*N)/chunk_length
    print('Degrees of freedom for red noise test is %5.2f' %(dof))
    specout_red = np.zeros(N)
    Nr = np.int(N/2)
    freq=np.arange(0,Nr,1)/N
    # This is the expected red noise spectrum: Equation (9.77) of Wilks 3rd edition (2011), using also equation (9.22)
    for ii in range(0,Nr,1):
        white_noise_var = (1.0-rr**2.0)*np.var(datam)
        #xr[ii] = ((4.0*white_noise_var)/Nr)/(1.0 + (rr**2.0) - (2.0*rr*np.cos(2.*np.pi*(freq[ii]))))
        specout_red[ii] = ((4.0*white_noise_var)/((Nr-1)/(Nr-2)))/(1.0 + (rr**2.0) - (2.0*rr*np.cos(2.*np.pi*(freq[ii]))))


    f95 = [3.84,3.00,2.60,2.37,2.21,2.10,2.01,1.94,1.88,1.83,1.75,1.67,1.57,1.52,1.46,1.39,1.38,1.32,1.22,1.00]
    n95 = [1,2,3,4,5,6,7,8,9,10,12,15,20,24,30,40,50,60,120,1000000]
    f99 = [6.63,4.61,3.78,3,32,3.02,2.80,2.64,2.51,2.41,2.31,2.18,2.04,1.88,1.79,1.70,1.59,1.56,1.47,1.32,1.00]
    n99 = n95
        
    if conf==95:
        fstat = f95[-1]
        NF = np.size(f95)
        for ii in range(NF):
            if (dof <= n95[ii]):       
                fstat = f95[ii-1]+(dof-n95[ii-1])*(f95[ii]-f95[ii-1])/(n95[ii]-n95[ii-1])
                break   
    elif conf == 99:  
        fstat = f99[-1]
        NF = np.size(f99)
        for ii in range(NF):
            if (dof <= n99[ii]):        
                fstat = f99[ii-1]+(dof-n99[ii-1])*(f99[ii]-f99[ii-1])/(n99[ii]-n99[ii-1])  
                break
  
    specout_red_sig = specout_red*fstat
    return specout_red, specout_red_sig, fstat
    
def xsection_inds(slon, slat, elon, elat, lons, lats, m):
    '''
    Returns the indicies for creating a cross section.
    Note: The indices returned are generally south-to-north

    Parameters
    ----------
    slon : scalar
        The longitude value of the starting point
    slat : scalar
        The latitude value of the starting point
    elon : scalar
        The longitude value of the ending point
    elat : scalar
        The latitude value of the ending point
    lons : 2D array_like
        The longitude values of the underlying dataset
    lats : 2D array_like
        The latitude values of the underlying dataset
    m : Basemap Instance
        The basemap instance used for plotting data on a map

    Returns
    -------
    xinds : numpy array
        The first dimension indices of the cross section
    yinds : numpy array
        The second dimension indices of the cross section

    '''
    
    import scipy.spatial as ss
    from PIL import Image, ImageDraw
    
    x, y = m(lons, lats)
    gpoints = zip(x.ravel(), y.ravel())
    gtree = ss.cKDTree(gpoints)
    sx, sy = m(slon, slat)
    ex, ey = m(elon, elat)
    pts = np.array([(sx, sy), (ex, ey)])
    dists, inds = gtree.query(pts, k=1, distance_upper_bound=100*1000.)
    xinds, yinds = np.unravel_index(inds, x.shape)
    pts = ((yinds[0], xinds[0]), (yinds[1], xinds[1]))
    grid = np.zeros_like(lons)
    img = Image.new('L', grid.shape[::-1], 0)
    ImageDraw.Draw(img).line(pts, fill=1, width=1)
    img = np.array(img)
    xinds, yinds = np.where(img > 0)
    if slat > elat:
        xinds = xinds[::-1]
        yinds = yinds[::-1]
    return xinds, yinds

def nbrInds8_ll(iLat, iLon, nLat, nLon):
  #return list of 8-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
  #always have east, west neighbors. not north/south at respective pole
  #at the north pole, all pole-1 latitude cells are nbrs.
  
  iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
  iEast = (iLon+1)%nLon
  iSouth = iLat+1
  iNorth = iLat-1
  
  haveSouth = iSouth<nLat
  haveNorth = iNorth>-1
  
  nbrLats = [iLat, iLat]; nbrLons = [iWest, iEast]
  if (haveSouth):
    nbrLats.extend([iSouth, iSouth, iSouth])
    nbrLons.extend([iWest, iLon, iEast])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iNorth]*nLon)
    nbrLons.extend(range(nLon))
    
  if (haveNorth):
    nbrLats.extend([iNorth, iNorth, iNorth])
    nbrLons.extend([iWest, iLon, iEast])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iSouth]*nLon)
    nbrLons.extend(range(nLon))
    
  #to get flat array indices, call index_2dTo1d(lats, lons, nlon)
  #nbrLats = np.array(nbrLats, dtype=int); nbrLons = np.array(nbrLons, dtype=int)
  return (nbrLats, nbrLons)

def nbrInds4_ll(iLat, iLon, nLat, nLon):
  #return list of 4-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
  #always have east, west neighbors. not north/south at respective pole
  #at the north pole, all pole-1 latitude cells are nbrs.
  
  iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
  iEast = (iLon+1)%nLon
  iSouth = iLat+1
  iNorth = iLat-1
  
  haveSouth = iSouth<nLat
  haveNorth = iNorth>-1
  
  nbrLats = [iLat, iLat]; nbrLons = [iWest, iEast]
  if (haveSouth):
    nbrLats.extend([iSouth])
    nbrLons.extend([iLon])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iNorth]*nLon)
    nbrLons.extend(range(nLon))
    
  if (haveNorth):
    nbrLats.extend([iNorth])
    nbrLons.extend([iLon])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iSouth]*nLon)
    nbrLons.extend(range(nLon))
    
  #to get flat array indices, call index_2dTo1d(lats, lons, nlon)
  #nbrLats = np.array(nbrLats, dtype=int); nbrLons = np.array(nbrLons, dtype=int)
  return (nbrLats, nbrLons)

def nbrInds_ll(iLat, iLon, nLat, nLon):
  #return nbrInds8_ll(iLat, iLon, nLat, nLon)
  return nbrInds4_ll(iLat, iLon, nLat, nLon)

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  '''
  #return the distance between 1 ll1 point and >=1 ll2 points.
  on a sphere.
  input lat/lon in degrees!
  '''
  
  d2r = np.pi/180.;
  lat1 *= d2r; lon1 *= d2r; lat2 *= d2r; lon2 *= d2r
  
  dlat = lat2-lat1
  dlon = lon2-lon1
  latTerm = np.sin(.5*dlat); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2);
  dAngle = np.sqrt(latTerm+lonTerm)
  dist = 2.*r*np.arcsin(dAngle)
  
  return(dist)

def wrapLongitudeRange(lon):
  #given longitude value, return within range [-180, 180)
  while (lon>180):
    lon = lon-360
  while (lon<-180):
    lon = lon+360
  return lon

def gatherCellsOnLine_coords(lat, lon, latCoord, lonCoord):
  #return a unique array of cells obtained from the coordinate cells
  
  nLat = len(lat); nLon=len(lon); nPts = len(latCoord)
  latCells = np.empty(nPts, dtype=int); lonCells = np.empty(nPts, dtype=int)
  for iPt in xrange(nPts):
    iLat = np.argmin(np.absolute(lat-latCoord[iPt]))
    iLon = np.argmin(np.absolute(lon-lonCoord[iPt]))
    latCells[iPt] = iLat; lonCells[iPt] = iLon
      
  flatInds = latCells*nLon+lonCells
  u, ind = np.unique(flatInds, return_index=True)
  flatInds = u[np.argsort(ind)]

  iLon = ind%nLon
  iLat = (ind-iLon)/nLon  
  return (flatInds, nLon)

def gatherCellsOnLine(iLat0, iLon0, lat1, lon1, lat, lon):
  #the "line" is defined by center cell and latLon of point to reach.
  #return a list of cells crossed to point..
  
  doDiskConn = False
  
  nLat = len(lat); nLon=len(lon)
  
  latInds = [iLat0]; lonInds = [iLon0]
  iLat = iLat0; iLon = iLon0
  dMin = calc_distSphere_multiple(1., lat[iLat0], lon[iLon0], lat1, lon1)
  while (True):
    if (doDiskConn):
      nbrLat, nbrLon = gatherCells_region(iLat, iLon, nLat, nLon, lat, lon, 
                         6371., 30.)
    else:
      nbrLat, nbrLon = nbrInds_ll(iLat, iLon, nLat, nLon)
    d = calc_distSphere_multiple(1., lat1, lon1, lat[nbrLat], lon[nbrLon])
    iNbr = np.argmin(d)
    if (d[iNbr]<dMin):
      iLat = nbrLat[iNbr]; iLon = nbrLon[iNbr]
      latInds.append(iLat); lonInds.append(iLon)
      dMin = d[iNbr]
    else:
      break
  
  return (latInds, lonInds)

def driver_xsec_old(lat, lon, cenLat, cenLon, dLat, dLon):
  #given tpv lat/lon, create a path ll-dll -> ll -> ll+dll
  
  #mesh info ------------------------------
  #lat = data.variables['latitude'][:];
  #lon = data.variables['longitude'][:];
  
  #id cross-section cells ------------------------------
  #if cross pole in latitude, have to wrap to other side of globe
  ll0 = np.array([cenLat-dLat, cenLon-dLon])
  if (ll0[0]<-90.):
    ll0[0] = -90.-(ll0[0]+90.) #need > -np.pi/2
    ll0[1] += 180.
  ll1 = np.array([cenLat, cenLon])
  ll2 = np.array([cenLat+dLat, cenLon+dLon])
  if (ll2[0]>90.):
    ll2[0] = 90.-(ll2[0]-90.) #need < np.pi/2
    ll2[1] += 180.
  #for lat/lon mesh here, longitudes also have to be bounded
  ll0[1] = wrapLongitudeRange(ll0[1])
  ll1[1] = wrapLongitudeRange(ll1[1])
  ll2[1] = wrapLongitudeRange(ll2[1])
  
  #note that the above will fail in purpose if you input a silly dLat,
  #eg, degrees vs. radians
  print (ll0, ll1, ll2)
  
  #go minus to center
  iLat0 = np.argmin(np.absolute(lat-ll0[0]))
  iLon0 = np.argmin(np.absolute(lon-ll0[1]))
  latCells0, lonCells0 = gatherCellsOnLine(iLat0, iLon0, ll1[0], ll1[1], lat, lon)
  
  #go center to plus
  iLat0 = np.argmin(np.absolute(lat-ll1[0]))
  iLon0 = np.argmin(np.absolute(lon-ll1[1]))
  latCells1, lonCells1 = gatherCellsOnLine(iLat0, iLon0, ll2[0], ll2[1], lat, lon)
  
  #don't double count the tpv cell
  cellsOnLine_lat = np.array(latCells0[0:-1]+latCells1)
  cellsOnLine_lon = np.array(lonCells0[0:-1]+lonCells1)
  
  return (cellsOnLine_lat, cellsOnLine_lon)

def genLatLonCoords(lat0, lon0, dLat, dLon, nPts):
  
  latCoords = np.linspace(lat0-dLat, lat0+dLat, nPts)
  lonCoords = np.linspace(lon0-dLon, lon0+dLon, nPts)

  print (latCoords)
  print (lonCoords)
  
  #recover lats that cross pole 1x
  crossedP = latCoords<-np.pi/2
  latCoords[crossedP] = -np.pi/2-(latCoords[crossedP]+np.pi/2)
  lonCoords[crossedP] += np.pi
  
  crossedP = latCoords>np.pi/2
  latCoords[crossedP] = np.pi/2-(latCoords[crossedP]-np.pi/2)
  lonCoords[crossedP] += np.pi
  
  
  lonCoords = lonCoords%(2*np.pi)
  
  return (latCoords, lonCoords)

def driver_xsec(lat, lon, cenLat, cenLon, dLat, dLon):
  
  
  #id cross-section cells ------------------------------
  r2d = 180./np.pi; d2r = np.pi/180.; nPts = 601
  latCoords, lonCoords = genLatLonCoords(cenLat*d2r, cenLon*d2r, dLat*d2r, dLon*d2r, nPts)
  latCoords *= r2d; lonCoords *= r2d
  cellsOnLine_lat, cellsOnLine_lon = gatherCellsOnLine_coords(lat, lon, latCoords, lonCoords)
  
  return (cellsOnLine_lat, cellsOnLine_lon)

def boundary_buffer(field, npts_x, npts_y):
    ''' Replace field with NaNs near boundary defined by npts_x and npts_y '''
    ''' '''
    ''' '''
    ''' field: Input data array '''
    ''' npts_x, npts_y: number of points from boundaries to replace with NaNs '''
    
    dimens = np.shape(field)
    
    if len(dimens)==2:
        iy,ix = np.shape(field)
        field[0:npts_y,:] = float('NaN')
        field[-(npts_y+1),:] = float('NaN')	
        field[:,0:npts_x] = float('NaN')
        field[:,-(npts_x+1)] = float('NaN')		
	
    if len(dimens)==3:
        iz,iy,ix = np.shape(field)    
        field[:,0:npts_y,:] = float('NaN')
        field[:,-(npts_y+1),:] = float('NaN')	
        field[:,:,0:npts_x] = float('NaN')
        field[:,:,-(npts_x+1)] = float('NaN')		    
    
    return field 

def find_matching_indices(a, b):   
    ''' Find the indices of array a that match the values in array b'''
    ''' '''
    ''' a: Input array to match values to '''
    ''' b: Input array with values that are a subset of array a '''
    ''' '''
    ''' Output: minds (indices of a that match b)'''
    ''' '''
    ''' Steven Cavallo '''
    ''' December 2014 '''
    
    ab = np.in1d(a.ravel(), b).reshape(a.shape)
    linds = np.where(ab)

    return linds
def date_to_jd(year,month,day):
    """
    Convert a date to Julian Day.
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
    4th ed., Duffet-Smith and Zwart, 2011.
    Parameters
    ----------
    year : int
    Year as integer. Years preceding 1 A.D. should be 0 or negative.
    The year before 1 A.D. is 0, 10 B.C. is year -9.
    month : int
    Month as integer, Jan = 1, Feb. = 2, etc.
    day : float
    Day, may contain fractional part.
    Returns
    -------
    jd : float
    Julian Day
    Examples
    --------
    Convert 6 a.m., February 17, 1985 to Julian Day
    >>> date_to_jd(1985,2,17.25)
    2446113.75
    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        D = math.trunc(30.6001 * (monthp + 1))
        jd = B + C + D + day + 1720994.5
    return jd
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
    4th ed., Duffet-Smith and Zwart, 2011.
    Parameters
    ----------
    jd : float
    Julian Day
    Returns
    -------
    year : int
    Year as integer. Years preceding 1 A.D. should be 0 or negative.
    The year before 1 A.D. is 0, 10 B.C. is year -9.
    month : int
    Month as integer, Jan = 1, Feb. = 2, etc.
    day : float
    Day, may contain fractional part.
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    """
    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    return year, month, day 
def smooth_onedim(x,npasses):
    """ data = smooth_onedim(x)

     Uses a 5-point moving average with filter on vector x 
     with coefficients equal to the reciprocal of the span 
     to obtain vector dat

     x = input vector
     npasses = number of smoothing passes
     
     Returns vector called data
     
     Steven Cavallo
     January 2015     
    """
    
    if npasses == 0:
       print ('Number of smoothing passes is set to zero.  Set to a value greater than zero to smooth your data.')
       data = x
       return data
    
    for tt in range(0,npasses):       
        if tt == 0:
            data = np.zeros_like(x).astype('f')
            npts = len(x)

        data[0] = x[1]
        data[1] = (x[0] + x[1] + x[2])/3    
        for ii in range(2,npts-2):	   
            data[ii] = (x[ii-2] + x[ii-1] + x[ii] + x[ii+1] + x[ii+2] ) /5

        data[npts-2] = (x[npts-3] + x[npts-2] + x[npts-1])/3
        data[npts-1] = x[npts-1]

        x = data
    
    return data

def mjd_to_jd(mjd):
    """
    Convert Modified Julian Day to Julian Day.
        
    Parameters
    ----------
    mjd : float
        Modified Julian Day
        
    Returns
    -------
    jd : float
        Julian Day
    
        
    """
    return mjd + 2400000.5

    
def jd_to_mjd(jd):
    """
    Convert Julian Day to Modified Julian Day
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    mjd : float
        Modified Julian Day
    
    """
    return jd - 2400000.5


def date_to_jd(year,month,day):
    """
    Convert a date to Julian Day.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
    
    Returns
    -------
    jd : float
        Julian Day
        
    Examples
    --------
    Convert 6 a.m., February 17, 1985 to Julian Day
    
    >>> date_to_jd(1985,2,17.25)
    2446113.75
    
    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        
    D = math.trunc(30.6001 * (monthp + 1))
    
    jd = B + C + D + day + 1720994.5
    
    return jd
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day
    
    
def hmsm_to_days(hour=0,min=0,sec=0,micro=0):
    """
    Convert hours, minutes, seconds, and microseconds to fractional days.
    
    Parameters
    ----------
    hour : int, optional
        Hour number. Defaults to 0.
    
    min : int, optional
        Minute number. Defaults to 0.
    
    sec : int, optional
        Second number. Defaults to 0.
    
    micro : int, optional
        Microsecond number. Defaults to 0.
        
    Returns
    -------
    days : float
        Fractional days.
        
    Examples
    --------
    >>> hmsm_to_days(hour=6)
    0.25
    
    """
    days = sec + (micro / 1.e6)
    
    days = min + (days / 60.)
    
    days = hour + (days / 60.)
    
    return days / 24.
    
    
def days_to_hmsm(days):
    """
    Convert fractional days to hours, minutes, seconds, and microseconds.
    Precision beyond microseconds is rounded to the nearest microsecond.
    
    Parameters
    ----------
    days : float
        A fractional number of days. Must be less than 1.
        
    Returns
    -------
    hour : int
        Hour number.
    
    min : int
        Minute number.
    
    sec : int
        Second number.
    
    micro : int
        Microsecond number.
        
    Raises
    ------
    ValueError
        If `days` is >= 1.
        
    Examples
    --------
    >>> days_to_hmsm(0.1)
    (2, 24, 0, 0)
    
    """
    hours = days * 24.
    hours, hour = math.modf(hours)
    
    mins = hours * 60.
    mins, min = math.modf(mins)
    
    secs = mins * 60.
    secs, sec = math.modf(secs)
    
    micro = round(secs * 1.e6)
    
    return int(hour), int(min), int(sec), int(micro)
    

def datetime_to_jd(date):
    """
    Convert a `datetime.datetime` object to Julian Day.
    
    Parameters
    ----------
    date : `datetime.datetime` instance
    
    Returns
    -------
    jd : float
        Julian day.
        
    Examples
    --------
    >>> d = datetime.datetime(1985,2,17,6)  
    >>> d
    datetime.datetime(1985, 2, 17, 6, 0)
    >>> jdutil.datetime_to_jd(d)
    2446113.75
    
    """
    days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
    
    return date_to_jd(date.year,date.month,days)
    
    
def jd_to_datetime(jd):
    """
    Convert a Julian Day to an `jdutil.datetime` object.
    
    Parameters
    ----------
    jd : float
        Julian day.
        
    Returns
    -------
    dt : `jdutil.datetime` object
        `jdutil.datetime` equivalent of Julian day.
    
    Examples
    --------
    >>> jd_to_datetime(2446113.75)
    datetime(1985, 2, 17, 6, 0)
    
    """
    year, month, day = jd_to_date(jd)
    
    frac_days,day = math.modf(day)
    day = int(day)
    
    hour,min,sec,micro = days_to_hmsm(frac_days)
    
    return datetime(year,month,day,hour,min,sec,micro)


def timedelta_to_days(td):
    """
    Convert a `datetime.timedelta` object to a total number of days.
    
    Parameters
    ----------
    td : `datetime.timedelta` instance
    
    Returns
    -------
    days : float
        Total number of days in the `datetime.timedelta` object.
        
    Examples
    --------
    >>> td = datetime.timedelta(4.5)
    >>> td
    datetime.timedelta(4, 43200)
    >>> timedelta_to_days(td)
    4.5
    
    """
    seconds_in_day = 24. * 3600.
    
    days = td.days + (td.seconds + (td.microseconds * 10.e6)) / seconds_in_day
    
    return days
def cartopy_plot_prep(proj,lonarr,latarr,plotvar):
    """
    Prepares a plot variable for plotting with Cartopy
    
    Steven Cavallo
    May 2023
    """
    import cartopy.crs as ccrs
    
    out = proj.transform_points(ccrs.PlateCarree(),lonarr,latarr,plotvar)
    x = out[...,0]
    y = out[...,1]
    z = out[...,2]    
    
    return x, y, z

def geodesic(crs, start, end, steps):
    '''
    Construct a geodesic path between two points.
    This function acts as a wrapper for the geodesic construction available in `pyproj`.
    Parameters
    ----------
    crs: `cartopy.crs`
        Cartopy Coordinate Reference System to use for the output
    start: (2, ) array_like
        A latitude-longitude pair designating the start point of the geodesic (units are
        degrees north and degrees east).
    end: (2, ) array_like
        A latitude-longitude pair designating the end point of the geodesic (units are degrees
        north and degrees east).
    steps: int, optional
        The number of points along the geodesic between the start and the end point
        (including the end points).
    Returns
    -------
    `numpy.ndarray`
        The list of x, y points in the given CRS of length `steps` along the geodesic.
    See Also
    --------
    cross_section
    '''
    import cartopy.crs as ccrs    

    g = Geod(crs.proj4_init)
    geodesic = np.concatenate([
        np.array(start[::-1])[None],
        np.array(g.npts(start[1], start[0], end[1], end[0], steps - 2)),
        np.array(end[::-1])[None]
    ]).transpose()
    points = crs.transform_points(ccrs.Geodetic(), *geodesic)[:, :2]
    return points
def earthDistance(lons1, lats1, lons2, lats2):
    rEarth = 6371
    lons1 = np.radians(lons1)
    lats1 = np.radians(lats1)
    lons2 = np.radians(lons2)
    lats2 = np.radians(lats2)
    dlon = lons2 - lons1
    dlat = lats2 - lats1
    a = np.sin(dlat / 2)**2 + np.cos(lats1) * np.cos(lats2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = rEarth * c
    return(distance)   
def random_color_generator(red_bounds,green_bounds,blue_bounds):
    import random
    from datetime import datetime

    random.seed(datetime.now())    
    
    r = random.randint(red_bounds[0], red_bounds[1])
    g = random.randint(green_bounds[0], green_bounds[1])
    b = random.randint(blue_bounds[0], blue_bounds[1])
    return (r, g, b)
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    [ind] = np.where(array == array[idx])
    
    return array[idx], ind
