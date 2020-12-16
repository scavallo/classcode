import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma

from mstats import *

fpath_in = '/arctic3/datasets/tpv_SH/tracks/tracks_low_horizPlusVert.nc'
# Variables in track file above:
# circ, vortMean, ampMaxMin, rEquiv, thetaVol, ampMean, thetaExtr, latExtr, lonExtr

# Load data in
data = netCDF4.Dataset(fpath_in,'r') # you will need to change the path to point to your data
trackLen = data.variables['lenTrack'][:]
Dates = data.variables['timeStamp'][:]
mons = np.array([int(str(x).replace('-','')[4:6]) for x in Dates])
years = np.array([int(str(x).replace('-','')[:4]) for x in Dates])
Dates = np.array([int(str(x).replace('-','')) for x in Dates])
trackStartDate = data.variables['iTimeStart'][:]
trackStartMon = mons[trackStartDate]
trackStartYear = years[trackStartDate]
trackStartDate = Dates[trackStartDate]


min_theta_winter = []
max_amp_winter = []

min_theta_summer = []
max_amp_summer = []

for x in range(0,len(trackLen)):
  if x%10000 == 0:
    print ("On track {0}/{1}".format(x,len(trackLen)))
  if trackLen[x] < 8: # checking to make sure TPV track was longer than two days
    continue

  lat = data.variables['latExtr'][x,:]
  if not ma.is_masked(lat):
    per_life_in_arctic = float(np.where(lat<=-60)[0].shape[0])/float(lat.shape[0]) # checking if TPV spent 60% of lifetime in Arctic
  else:
    per_life_in_arctic = float(np.where((lat.data<=-60)&(lat.mask!=True))[0].shape[0])/float(np.where((lat.mask!=True))[0].shape[0])
  if per_life_in_arctic*100 < 60.:
    continue

  if trackStartMon[x] == 12 or trackStartMon[x] <= 2: # only getting tpv tracks for winter months
    min_theta_summer.append(np.amin(data.variables['thetaExtr'][x,:]))
    max_amp_summer.append(np.amax(data.variables['ampMaxMin'][x,:]))
  elif trackStartMon[x] > 5 and trackStartMon[x] < 9:
    min_theta_winter.append(np.amin(data.variables['thetaExtr'][x,:]))# only getting tpv tracks for summer months
    max_amp_winter.append(np.amax(data.variables['ampMaxMin'][x,:]))
  else:
    continue

data.close()

min_theta_winter = np.array(min_theta_winter)
min_theta_summer = np.array(min_theta_summer)
max_amp_winter = np.array(max_amp_winter)
max_amp_summer = np.array(max_amp_summer)


#########################################################################
################ Plotting section #######################################

fig,ax0 = plt.subplots(1,1)
labels = ['Winter','Summer']
parts = ax0.violinplot([min_theta_winter,min_theta_summer],showmeans=False,showmedians=False,showextrema=False)
for pc in parts['bodies']:
    pc.set_facecolor('#1b9e77')
    pc.set_edgecolor('#7570b3')
    pc.set_alpha(1)

quartile1_winter, medians_winter, quartile2_winter = np.percentile(min_theta_winter, [25, 50, 75], axis=0)
quartile1_summer, medians_summer, quartile2_summer = np.percentile(min_theta_summer, [25, 50, 75], axis=0)

whiskermin_winter,whiskermax_winter = np.percentile(min_theta_winter, [5,95], axis=0)
whiskermin_summer,whiskermax_summer = np.percentile(min_theta_summer, [5,95], axis=0)

medians = [medians_winter,medians_summer]
quartile1 = [quartile1_winter,quartile1_summer]
quartile2 = [quartile2_winter,quartile2_summer]
whiskersmin = [whiskermin_winter,whiskermin_summer]
whiskersmax = [whiskermax_winter,whiskermax_summer]

inds = np.arange(1, len(medians) + 1)
ax0.scatter(inds, medians, marker='_', color='yellow', s=30, zorder=6)
ax0.vlines(inds, quartile1, quartile2, color='k', linestyle='-', lw=5,zorder=3)
ax0.vlines(inds, whiskersmin, whiskersmax, color='k', linestyle='-', lw=1,zorder=3)

ax0.set_xticks(inds)
ax0.set_xticklabels(labels,rotation='vertical')
ax0.set_title('Violin Plots - Min Theta of TPVs')
ax0.set_ylabel('K')
plt.show()
#plt.savefig('violin_min_theta.png',dpi=300,bbox_inches='tight')
###############################################################################
fig,ax0 = plt.subplots(1,1)
labels = ['Winter','Summer']
parts = ax0.violinplot([max_amp_winter,max_amp_summer],showmeans=False,showmedians=False,showextrema=False)
for pc in parts['bodies']:
    pc.set_facecolor('#1b9e77')
    pc.set_edgecolor('#7570b3')
    pc.set_alpha(1)

quartile1_winter, medians_winter, quartile2_winter = np.percentile(max_amp_winter, [25, 50, 75], axis=0)
quartile1_summer, medians_summer, quartile2_summer = np.percentile(max_amp_summer, [25, 50, 75], axis=0)

whiskermin_winter,whiskermax_winter = np.percentile(max_amp_winter, [5,95], axis=0)
whiskermin_summer,whiskermax_summer = np.percentile(max_amp_summer, [5,95], axis=0)

medians = [medians_winter,medians_summer]
quartile1 = [quartile1_winter,quartile1_summer]
quartile2 = [quartile2_winter,quartile2_summer]
whiskersmin = [whiskermin_winter,whiskermin_summer]
whiskersmax = [whiskermax_winter,whiskermax_summer]

inds = np.arange(1, len(medians) + 1)
ax0.scatter(inds, medians, marker='_', color='yellow', s=30, zorder=6)
ax0.vlines(inds, quartile1, quartile2, color='k', linestyle='-', lw=5,zorder=3)
ax0.vlines(inds, whiskersmin, whiskersmax, color='k', linestyle='-', lw=1,zorder=3)

ax0.set_xticks(inds)
ax0.set_xticklabels(labels,rotation='vertical')
ax0.set_title('Violin Plots - Max Amplitude of TPVs')
ax0.set_ylabel('K')
plt.show()
#plt.savefig('violin_max_amp.png',dpi=300,bbox_inches='tight')
##################################################################################






