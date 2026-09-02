#!/usr/bin/python 

# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import math
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *

from scipy import ndimage
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

rEarth = 6370.e3 #radius of spherical Earth (m)
def read_tracks_metrics(fNameTracks, metricNames):
  #Input the name of the track file and list of metricNames strings.
  #return list of numpy arrays list[iTrack][iTime,iMetric] with the metric properties of the tracked TPVs.
  
  nMetrics = len(metricNames)
  data = netCDF4.Dataset(fNameTracks,'r')
  nTracks = len(data.dimensions['nTracks'])
  
  trackList = []
  timeStartList = []
  
  for iTrack in xrange(nTracks):
    nTimes = data.variables['lenTrack'][iTrack]
    iTimeStart = data.variables['iTimeStart'][iTrack]
    timeStartList.append(data.variables['timeStamp'][iTimeStart])
    
    trackVals = np.empty((nTimes,nMetrics),dtype=float)
    for iMetric in xrange(nMetrics):
        key = metricNames[iMetric]
        trackVals[:,iMetric] = data.variables[key][iTrack,0:nTimes]
    trackList.append(trackVals)
        
  data.close()
  return trackList, timeStartList
  
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
