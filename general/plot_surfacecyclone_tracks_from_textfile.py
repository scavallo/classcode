#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

import tpv_tracker_plotting_modules as tpvmod

from scipy import ndimage
import scipy.io as sio
import math
import matplotlib.colors as col

import warnings
warnings.filterwarnings("ignore")



###################################
# Set user options
###################################
date_firstrecord = '2001083100' # used only if plot_tpvs_within_timerange = 'True'
date_lastrecord = '2001093118' # used only if plot_tpvs_within_timerange = 'True'
plot_tpvs_within_timerange = 'True' # If 'True', only finds vortices within time range listed above in date_firstrecord and date_lastrecord 
plot_tpvs_at_specific_time = 'False' # If 'True', finds all vortices located at time date_firstrecord

plot_cyclones_inbox = 'True' # If 'True', only finds TPVs that pass through a box centered at box_center_latlon and a distance of box_thresh from that lat/lon coordinate pair

box_center_latlon = [48., -51.9] # ending TC location
box_thresh = 500 # km; used only if plot_cyclones_inbox = 'True'



min_lifetime_days = 2 # finds only vortices within minimum lifetime
lat_origin_range = [20.,90.] # Vortex origin range.  I.e., if you want everything in Northern Hemisphere, set to [0.,90.] or if you want Arctic, set to [60.,90.]
make_tpv = 'False' # Restricts vortices to satisify definition of TPV: Must spend 60% of lifetime poleward of 65N latitude with a 2 day minimum lifetime


map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
zoom = 'true'
highlight_state = 'False'
proj_latlon = [40. , 270.]
#proj_latlon = [90. , 310.]
hinc = 6 # number of hours between records;  only used if record_num > -1

plot_metric_option = 1 # 1 for minimum slp, 2 for amplitude, 3 for radius, 4 for circulation, 
                       # 10 for genesis/lysis points, 11 for locations of maximum amplitudes, 12 for vortex speed, 13 for vortex direction,  99 for tracks only
		       # 100 for SLP
overlay_genesis_lysis = 'False'
overlay_genesis_only = 'True'
overlay_lysis_only = 'False'
overlay_max_locations = 'False'

label_fontsize = 18
track_thickness = 2
num_passes_smooth = 1

#fdir = '/data1/scavallo/data/seaice/seaice_loss_events/5percentile/'
#fdir = '/data2/scavallo/nnrp/cyclones_nsidc/'
fdir = '/data2/scavallo/nnrp/cyclones_nsidc/'
# Change the above to point to the track files located currently in the path below
# fdir = '/raid3/datasets/tpv_tracks/0p5/'

textfile = 'ncepstorms_allarctic_1958_2016.txt'
textfile_out = 'ncep_sfccylones_2001083100_2001093118.dat'


imagedir = '/home/scavallo/scripts/python_scripts/images/'
figname_prefix = 'sfcccylone_tracks_case_20010901'
###################################
# END user options
###################################
rEarth = 6370.e3 #radius of spherical Earth (m)
def colorline(x, y, z=None, cmap=plt.get_cmap('Blues_r'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0, ax=plt.gca()):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    

    segments = make_segments(x, y)


    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    #ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments

def index_2dTo1d(iLat, iLon, nLon):
  return iLat*nLon+iLon
  
def index_1dTo2d(ind, nLon):
  iLon = ind%nLon
  iLat = (ind-iLon)/nLon
  return (iLat, iLon)
  
def get_latLon_inds(self, inds):
    iLats, iLons = helpers.index_1dTo2d(inds, self.nLon)
    return (self.lat[iLats], self.lon[iLons])   

def flatten_2dTo1d(vals, nLat, nLon):
  #vals[lat][lon] goes to vals[iLat*nLon+iLon]
  
  valsOut = np.ravel(vals)
  return valsOut

def unflatten_1dTo2d(vals, nLat, nLon):
  #vals[iLat*nLon+iLon] goes to vals[lat][lon]
  
  valsOut = np.reshape(vals, (nLat,nLon))
  return valsOut

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  '''
  #return the distance between 1 ll1 point and >=1 ll2 points.
  on a sphere.
  input lat/lon in radians!!
  '''
  
  dlat = lat2-lat1
  dlon = lon2-lon1
  latTerm = np.sin(.5*dlat); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2);
  dAngle = np.sqrt(latTerm+lonTerm)
  dist = 2.*r*np.arcsin(dAngle)
  
  return dist

def distance_on_unit_sphere(lat1, long1, lat2, long2):

    '''
    Shamelessly taken from: http://www.johndcook.com/blog/python_longitude_latitude/
    '''
 
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) =
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    return arc
def mkcmap(): 
    white = '#ffffff'
    black = '#000000'
    red = '#ff0000'
    blue = '#0000ff'
    anglemap = col.LinearSegmentedColormap.from_list(
        'anglemap', [black, red, white, blue, black], N=256, gamma=1)
    return anglemap

anglemap = mkcmap()



fpath = fdir + textfile
aa = np.loadtxt(fpath, skiprows=1)       
datelist = aa[:,0]
lat = aa[:,1]	
lon = aa[:,2]
thetamin = aa[:,3]	    
thetaamp = aa[:,4]
vort_circ = aa[:,5]	    
vort_radius = aa[:,6]

nrows, ncols = np.shape(aa)
print nrows, ncols




if plot_metric_option == 1:
   varMin = 960.; varMax = 1008.
   #cmap_opt = plt.cm.RdBu_r
   cmap_opt = plt.cm.jet
elif plot_metric_option == 2:
   varMin = 5.; varMax = 25.
   cmap_opt = plt.cm.Blues
elif plot_metric_option == 3:
   varMin = 200.; varMax = 1000.
   cmap_opt = plt.cm.gist_heat_r
elif plot_metric_option == 4:
   varMin = 10.; varMax = 100.
   cmap_opt = plt.cm.gist_heat_r
elif plot_metric_option == 10: 
   overlay_genesis_lysis = 'True'
elif plot_metric_option == 12: 
   overlay_genesis_lysis = 'False'
   varMin = 5.; varMax = 40.
   cmap_opt = plt.cm.gist_heat_r   
elif plot_metric_option == 13: 
   overlay_genesis_lysis = 'False'
   varMin = 0.; varMax = 360.
   cmap_opt = plt.cm.gist_rainbow_r

   
   
   
elif plot_metric_option == 99:
   varMin = 260.; varMax = 325.
   cmap_opt = plt.cm.Greys_r
elif plot_metric_option == 100: 
   overlay_genesis_lysis = 'False'
   varMin = 980.; varMax = 1000.
   #cmap_opt = plt.cm.RdBu_r   
   cmap_opt = plt.cm.gist_heat_r

 
# Set global figure properties
golden = (np.sqrt(5)+1.)/2.
figprops = dict(figsize=(8., 16./golden), dpi=128)
fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#m = Basemap(projection='ortho',lon_0=0,lat_0=89.5, resolution='l')
if map_projection == 'ortho':
   if zoom == 'false':   
       m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
          resolution = 'l', area_thresh = 1000.,ax=ax1)
   else:
       m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
        	resolution = 'l', area_thresh = 1000.,ax=ax1)      		           

       width = m1.urcrnrx - m1.llcrnrx
       height = m1.urcrnry - m1.llcrnry

       coef = 0.7
       #coef = 0.9
       width = width*coef
       height = height*coef
       m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
           llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)		  
elif map_projection == 'lcc':
#   m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
   m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
       rsphere=(6378137.00,6356752.3142),\
       resolution='l',area_thresh=1000.,projection='lcc',\
       lat_1=50.,lon_0=-107.,ax=ax1)	       
elif map_projection == 'npstere':
   if zoom == 'false':
      m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
   else:
      m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')

#ax = plt.figure()
ax = plt.gca()
#m.drawcoastlines()
m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
m.drawcountries(linewidth=1, color='#444444', zorder=5)
m.drawstates(linewidth=0.66, color='#444444', zorder=4)
if ( (plot_metric_option != 20) and (plot_metric_option != 100) ):
   m.fillcontinents(color='Wheat',lake_color='lightblue', zorder=1) 
   m.drawmapboundary(fill_color='lightblue')     


m.drawmeridians(np.arange(0, 360, 30),labels=[True,True,True,True])
m.drawparallels(np.arange(-90, 91, 30),labels=[True,True,False,False])



   

dumval = -999.00
fpath_out = fdir + textfile_out
outfile = open(fpath_out,'w')
outfile.write('date lat lon thetamin amplitude circulation radius')
outfile.write('\n')

ntracks_min = (24./hinc)*min_lifetime_days
tt = 0
while tt < nrows:
    if datelist[tt] < 10000:
        ntracks = datelist[tt]
	sind = tt+1
	eind = np.int(sind+ntracks-1)

	datesnow = datelist[sind:eind]
	
	
        datestrinit = str(datelist[sind])	
	yyyyinit = datestrinit[0:4]
	mminit = datestrinit[4:6]
	ddinit = datestrinit[6:8]
	hhinit = datestrinit[8:10]
        
	try:
	    datestrfin = str(datelist[eind-1])	
	except:
	    #eind = np.int(sind+ntracks-1)
	    datestrfin = str(datelist[eind])	
	    print 'hi'
	yyyyfin = datestrfin[0:4]
	mmfin = datestrfin[4:6]
	ddfin = datestrfin[6:8]
	hhfind = datestrfin[8:10]	
        
	
	latnow = lat[sind:eind]
	lonnow = lon[sind:eind]
	thetanow = thetamin[sind:eind]
	ampnow = thetaamp[sind:eind]
	circnow = vort_circ[sind:eind]
	radnow = vort_radius[sind:eind] 
	
	speednow = np.zeros(len(latnow))	
	speednow[:] = float('NaN')
	dirnow = np.zeros(len(latnow))	
	dirnow[:] = float('NaN')	
	       
	if ( (plot_tpvs_within_timerange == 'True') and (plot_tpvs_at_specific_time == 'False') ):
            plotTrack1 = False
	    #plotTrack = True
            if ( int(datestrinit>=date_firstrecord) and  ( int(datestrfin<=date_lastrecord) ) ):	
		plotTrack1 = True
		print datestrinit, datestrfin
        else:
	    plotTrack1 = True
	
	if ((plot_tpvs_at_specific_time == 'True') and (plot_tpvs_within_timerange == 'False')):
            plotTrack1 = False
	    #plotTrack = True
            if ( int(datestrinit<=date_firstrecord) and  ( int(datestrfin>=date_firstrecord) ) ):	
		plotTrack1 = True
		print datestrinit, datestrfin
        else:
	    #plotTrack1 = True
	    aa = 1

        latnow = lat[sind:eind]
	lonnow = lon[sind:eind]	    
	lonnow = um.filter_numeric_nans(lonnow,361.,float('NaN'),'high')
	lonnow = um.filter_numeric_nans(lonnow,-1.,float('NaN'),'low')
	latnow = um.filter_numeric_nans(latnow,91.,float('NaN'),'high')
	latnow = um.filter_numeric_nans(latnow,-91.,float('NaN'),'low')
	x,y = m(lonnow,latnow)
        
	if plot_cyclones_inbox == 'True':			
	    plotTrack2 = False
	    for ii in xrange(0,len(latnow)):
		if lonnow[ii]<0:
	            lonnow[ii] = lonnow[ii] + 360.			 

		distance = abs(tpvmod.rEarth*tpvmod.distance_on_unit_sphere(box_center_latlon[0],box_center_latlon[1],latnow[ii],lonnow[ii]))
		distance = distance/1000.
        	if distance <= box_thresh:
	            plotTrack2 = True		    
		    break; break
	else:
	    plotTrack2 = True    
        
	if plot_metric_option == 12:	
	    for ii in xrange(1,len(latnow)):
		if lonnow[ii]<0:
	            lonnow[ii] = lonnow[ii] + 360.			 

		try:
		    distnow = abs(tpvmod.rEarth*tpvmod.distance_on_unit_sphere(latnow[ii],lonnow[ii],latnow[ii-1],lonnow[ii-1]))
		except:
		    distnow = 0.
		#distance = distance/1000.	
		speednow[ii] = distnow / (6.0*3600.)
	if plot_metric_option == 13:	
	    for ii in xrange(1,len(latnow)):
		if lonnow[ii]<0:
	            lonnow[ii] = lonnow[ii] + 360.			 

				    
		dlatnow = latnow[ii] - latnow[ii-1]
		dlonnow = lonnow[ii] - lonnow[ii-1]
		#distnow = math.atan2(dlatnow*(np.pi/180.),dlonnow*(np.pi/180.))
		#distnow = distnow*(180./np.pi)+90.
		distnow = math.atan2(dlatnow*(np.pi/180.),dlonnow*(np.pi/180.))
		if distnow < 0:
		    distnow = 90.-distnow*(180./np.pi)
		else:
		    distnow = 90.-np.abs(distnow*(180./np.pi))
		    if distnow < 0:
		        distnow = 360.- distnow
			
		#print dlatnow,dlonnow,distnow
		
		if ( (distnow <= 360.) and (distnow >= 0.) ):
		    dirnow[ii] = distnow 		
	        else:
		    dirnow[ii] = float('NaN')
	 
	if make_tpv == 'True':
	    nlats_polar = np.where(latnow>=65.0)
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
	
        
        if ( (lat[sind] >= lat_origin_range[0]) and (lat[sind] <= lat_origin_range[1]) and (ntracks >= ntracks_min) and (perc >= perc_thresh) and (plotTrack1 == True) and (plotTrack2 == True)):
	    
	    if (plotTrack2 == True):
	        #print ntracks,ntracks_min
	        print datestrinit, datestrfin

	    #plot track, with color representing value
	    #m.plot(x,y, 'b-')
	    if plot_metric_option == 1:
	        vals = thetamin[sind:eind]
		xlabunits = 'Kelvin'
		figname_suff = 'core'
	    elif plot_metric_option == 2:
	         vals = thetaamp[sind:eind]
		 xlabunits = 'Kelvin'
		 figname_suff = 'amps'
	    elif plot_metric_option == 3:
	         vals = vort_radius[sind:eind]
	         xlabunits = 'km'
		 figname_suff = 'radii'
	    elif plot_metric_option == 4:
	         vals = vort_circ[sind:eind]
	         xlabunits = r'km$^2$ s$^{-1}$'	 
		 figname_suff = 'circ'  
	    elif plot_metric_option == 10:
	         figname_suff = 'genesis_lysis'
	    elif plot_metric_option == 11:
	    	 figname_suff = 'maxamp_locations'	 
	    elif plot_metric_option == 12:
	         vals = speednow[:]		 
	         xlabunits = r'm s$^{-1}$'	 
		 figname_suff = 'speed'  
	    elif plot_metric_option == 13:
	         vals = dirnow[:]	
	         xlabunits = r'$^{\circ}$'	 
		 figname_suff = 'direction_deg'  
	    elif plot_metric_option == 99:
	    	 vals = thetamin[sind:eind]
		 vals[:] = 300.
		 varMin = 300.
		 varMax = 300.
	    	 figname_suff = 'tracks'	 	    
	    elif plot_metric_option == 100:
	    	 vals = thetamin[sind:eind]
		 xlabunits = 'hPa'		 	 
		 #vals[:] = 300.
		 varMin = 972.
		 varMax = 1008.
	    	 figname_suff = 'sfccyclone_tracks'
		 	    
            #print latnow, lonnow, vals
	    if plot_metric_option < 10:	         
	         CS3 = colorline(x, y, z=vals, cmap=cmap_opt, norm=plt.Normalize(varMin, varMax), linewidth=track_thickness, alpha=1.0, ax=ax)
 	    if ( (plot_metric_option == 12) or (plot_metric_option == 33) or (plot_metric_option == 99) or (plot_metric_option == 100)):
	         CS3 = colorline(x, y, z=vals, cmap=cmap_opt, norm=plt.Normalize(varMin, varMax), linewidth=3, alpha=1.0, ax=ax)
 	    if ( (plot_metric_option == 13) ):
	         CS3 = colorline(x, y, z=vals, cmap=anglemap, norm=plt.Normalize(varMin, varMax), linewidth=track_thickness, alpha=1.0, ax=ax)		 	    
	    if overlay_genesis_lysis == 'True':
	        #mark beginning and ending of track
		if (  (plot_metric_option == 99) ):
  	            m.scatter(x[0],y[0], marker='o', linewidth=6, color='0.5', s=45, zorder=2)
		else:
		    m.scatter(x[0],y[0], marker='o', linewidth=6, color='m', s=45, zorder=2)
	        m.scatter(x[-1],y[-1], marker='s',linewidth=6, color='g', s=45, zorder=2)
	    if overlay_genesis_only == 'True':
	        #mark beginning and ending of track
		if (  (plot_metric_option == 99) ):
  	            m.scatter(x[0],y[0], marker='o', linewidth=6, color='0.5', s=45, zorder=2)
		else:
		    m.scatter(x[0],y[0], marker='o', linewidth=6, color='m', s=45, zorder=2)
	    if overlay_lysis_only == 'True':
	        #mark beginning and ending of track
		if (  (plot_metric_option == 99) ):
  	            m.scatter(x[-1],y[-1], marker='s', linewidth=6, color='0.5', s=45, zorder=2)
		else:
		    m.scatter(x[-1],y[-1], marker='s', linewidth=6, color='m', s=45, zorder=2)
		    	           
	    if ( (plot_metric_option == 11) or (plot_metric_option == 2) ):
	        vals = thetaamp[sind:eind]
	        indnow = np.where(vals==np.max(vals))
		if overlay_max_locations == 'True':
  	            m.scatter(x[indnow],y[indnow], marker='s', linewidth=6, color='m', s=45, zorder=2)	            
	    if ( (plot_metric_option == 3)):	        
	        indnow = np.where(vals==np.max(vals))
		if overlay_max_locations == 'True':
  	            m.scatter(x[indnow],y[indnow], marker='s', linewidth=6, color='m', s=45, zorder=2)	            
	    if ( (plot_metric_option == 100)):
	        cm = plt.get_cmap('RdBu_r')
	        #CS3 = plt.scatter(x, y, s=100, c=vals, vmin=varMin,vmax=varMax,cmap=cm,alpha=0.5)
		if num_passes_smooth > 0:
		    x = um.smooth_onedim(x,num_passes_smooth)
		    y = um.smooth_onedim(y,num_passes_smooth)
		CS3 = colorline(x, y, z=vals, cmap=cmap_opt, norm=plt.Normalize(varMin, varMax), linewidth=track_thickness, alpha=1.0, ax=ax)
	    
	    wcount = 0				
	    outfile.write('%-10d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n' % (ntracks, dumval, dumval, dumval, dumval, dumval, dumval))
	    while wcount < ntracks-1:
	        

		datestrnow = str(datesnow[wcount])	
		yyyynow = datestrnow[0:4]
		mmnow = datestrnow[4:6]
		ddnow = datestrnow[6:8]
		hhnow = datestrnow[8:10]

		outfile.write('%-10s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n' % (str(int(datesnow[wcount])), latnow[wcount], lonnow[wcount], thetanow[wcount], ampnow[wcount], circnow[wcount], radnow[wcount]))
        	#timenow = um.advance_time(timenow,6)
        	wcount += 1

    tt += 1

if (plot_cyclones_inbox == 'True'):
    x1,y1 = m(box_center_latlon[1],box_center_latlon[0])
    m.scatter(x1,y1, marker='*', linewidth=4, color='0.5', s=220, zorder=2)
if ( (plot_metric_option < 10) or (plot_metric_option == 12) or (plot_metric_option == 13) or (plot_metric_option == 100)):
    cbar = plt.colorbar(CS3, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
    labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
    cbar.ax.set_xticklabels(labels, size=label_fontsize)
    cbar.set_label(xlabunits,size=label_fontsize)
    

if highlight_state == 'True':
    from matplotlib.patches import Polygon
    m.readshapefile('st99_d00', name='states', drawbounds=True)

    # collect the state names from the shapefile attributes so we can
    # look up the shape obect for a state by it's name
    state_names = []
    for shape_dict in m.states_info:
        state_names.append(shape_dict['NAME'])

    ax = plt.gca() # get current axes instance

    # get OK and draw the filled polygon
    seg = m.states[state_names.index('Oklahoma')]
    poly = Polygon(seg, facecolor='green',edgecolor='green',zorder=3)
    ax.add_patch(poly)


    #xx,yy = zip(*seg)  
    #color = 'g'
    #fill(xx,yy,color,edgecolor=color)
fignamenow = figname_prefix + '_' + figname_suff + '_' + date_firstrecord + '_' + date_lastrecord + '.png'
save_name = imagedir + fignamenow 
plt.savefig(save_name, bbox_inches='tight')


plt.show()
