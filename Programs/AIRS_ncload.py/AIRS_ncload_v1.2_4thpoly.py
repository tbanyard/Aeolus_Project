#!/usr/bin/env python3
"""
AIRS data load from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---vg55perts---Granule 55 Temperature Perturbations---------------------
----------[BRANCH]-This_is_a_working_branch_version_of_this_file--------
------------------------------------------------------------------------
========================================================================
Reads .nc files of AIRS data from UBPC-2027://media/GWR/AIRS/
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import matplotlib.text as text
import os
# os.putenv('CODA_DEFINITION',
# '/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
os.putenv('CODA_DEFINITION',
'/usr/local/share/coda/definitions/AEOLUS-20200227.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.io import savemat
from scipy.signal import savgol_filter, butter
import scipy.ndimage as ndimage
from itertools import groupby
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Change current working directory to parent directory
os.chdir('..')

# cmaps
qbocmap = LinearSegmentedColormap('QBOcustomcmap', segmentdata=customcolormaps('QBOcmap7'), N=265)

strdirectory = '/home/tpb38/PhD/Bath/AIRS/3d_airs_2019/207/'
filename = 'airs_2019_207_055.nc'
infile = strdirectory + str(filename) # Specifies data file
data = nc.Dataset(infile)

# Longitude
data_lon = data.variables['l1_lon'][:]
# Latitude
data_lat = data.variables['l1_lat'][:]
# Altitude
data_alt = data.variables['ret_z'][:]
# Temperature
data_temp = data.variables['ret_temp'][:]

# ~ print("data_alt: ", data_alt)
# ~ print("shape of data_lon: ", np.shape(data_lon))
# ~ print("shape of data_lat: ", np.shape(data_lat))
# ~ print("shape of data_alt: ", np.shape(data_alt))
# ~ print("shape of data_temp: ", np.shape(data_temp))

data.close()

# ~ """Calculating perturbations"""
# ~ data_temp2 = savgol_filter(data_temp, 55, 2, axis = 1)
# ~ data_temp = data_temp - data_temp2
# ~ data_temp = ndimage.gaussian_filter(data_temp, sigma=0.75, order=0)

"""Calculating 4th order polynomial perturbations"""
perturbations = np.zeros((135, 90))
print(np.shape(perturbations))
for irow in range(90):
	raw = data_temp[:,irow,10]
	p = np.polyfit(data_lat[:, irow], data_temp[:,irow,10], 4)
	BG = np.polyval(p, data_lat[:, irow])
	iperts = raw-BG
	perturbations[:, irow] = iperts
	
print(perturbations)

print(np.shape(data_temp))
print(np.shape(data_lat))
print(np.shape(data_lon))
print(np.shape(data_alt))

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')
os.chdir('AIRS')

# Initialise figure
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

# Plotting map
bmlowlat, bmupperlat, bmleftlon, bmrightlon = -70, -40, -90, -50
map = Basemap(projection='cyl',llcrnrlat=bmlowlat,urcrnrlat=bmupperlat,\
			llcrnrlon=bmleftlon,urcrnrlon=bmrightlon,resolution='i', ax=ax1)
# ~ map = Basemap(projection='cyl',llcrnrlat=-80,urcrnrlat=-40,\
			# ~ llcrnrlon=-80,urcrnrlon=-40,resolution='i', ax=ax2)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
map.drawcoastlines(linewidth=0.25, color='black', zorder=3)
map.drawmeridians([-80,-70, -60, -50], linewidth=0.3)
map.drawparallels([-70, -60, -50], linewidth=0.3)
map.drawmeridians([-70,-60, -40, -20], linewidth=0.3)
map.drawparallels([20, 30, 40, 50, 60], linewidth=0.3)

# Plotting
# ~ cs = plt.contourf(data_lon, data_lat, data_temp[:,:,10], cmap='RdBu_r', zorder=2, vmin = 180, vmax = 230, levels=np.linspace(180,230,51))
cs = plt.contourf(data_lon, data_lat, perturbations, cmap=qbocmap, zorder=2, vmin = -15, vmax = 15, levels=np.linspace(-15,15,13), extend='both')

# Fix axes
ax1.set_xticks([-90,-80, -70, -60, -50])
ax1.set_yticks([-70, -60, -50, -40])
ax1.set_xlabel('Longitude / deg')
ax1.set_ylabel('Latitude / deg')
ax1.set_aspect('auto') # Stretch map to fill subplot
for axis in ['top','bottom','left','right']: # Set axes thickness
	ax1.spines[axis].set_linewidth(0.75)

# Title
plt.title('AIRS granule 55 for 2019-07-26 at 30 km')

# Add colorbar to figure
fig.subplots_adjust(bottom=0.3, right=0.88, left=0.12)
cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.025])
# ~ fig.colorbar(cs, cmap='magma', orientation='horizontal',
	# ~ label='Temperature / K', cax=cbar_ax, ticks=np.linspace(180,230,11))
fig.colorbar(cs, cmap=qbocmap, orientation='horizontal',
	label='Temperature Perturbation / K', cax=cbar_ax, ticks=np.linspace(-15,15,7))



# Saving figure
plt.savefig('AIRSg55perts_4thpoly.png',dpi=300)
print(os.getcwd())

plt.close()
