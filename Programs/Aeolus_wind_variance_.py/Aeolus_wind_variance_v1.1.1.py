#!/usr/bin/env python3
"""
Aeolus data load for wind variance plots from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Fixing array operations and indexing issues-------------------
---v1.1.1-SP-stereographic projection
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
produces geographical plots of wind variance from the datasets
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
os.putenv('CODA_DEFINITION',
'/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.io import savemat
from itertools import groupby
from mpl_toolkits.basemap import Basemap

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import ncload

# Change current working directory to parent directory
os.chdir('..')

# Set desired global variables
res = 5 # Resolution of density plot in degrees
resrcp = 1/res # Resolution Reciprocal

# Check that the input resolution is acceptable
if 180 % res != 0:
	raise ValueError('Invalid resolution. Must be both a factor of 180 \
and non-residual in binary.')

# Calculate the number of latitude and longitude points
num_lat = int(180 * resrcp)+1
# ~ print('num_lat: ', num_lat)
num_lon = int(360 * resrcp)+1
# ~ print('num_lon: ', num_lon)

# Initialise all arrays
vals = np.zeros(100000)
vals[:] = np.nan
vals = np.zeros(1, dtype='object')
itrn = np.zeros(1)

# Bin into a 5deg grid
# Build lat/lon grid for variances to be placed in:
data_vars = np.asarray([[0 for _ in np.linspace(-180,180,num_lon)] \
for _ in np.linspace(-90,90,num_lat)])
# Build lat/lon grids
lat_grid = np.asarray([_ for _ in np.linspace(-90,90,num_lat)])
lon_grid = np.asarray([_ for _ in np.linspace(-180,180,num_lon)])
# Build array for all values for each week
data_values = np.asarray([[vals for _ in np.linspace(-180,180,num_lon)] \
for _ in np.linspace(-90,90,num_lat)]) # itrn_array[latitude][longitude]
# Build lat/lon grid for variance_iterations to be placed in:
data_itrn = np.asarray([[0 for _ in np.linspace(-180,180,num_lon)] \
for _ in np.linspace(-90,90,num_lat)])
np.set_printoptions(threshold=sys.maxsize)

# Here I need to iterate through all. nc files and plot all of them
# into jpgs to view one after another
"""Find directory and read netCDF data"""
"""strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC3/'
file_itrn = 0
directory = os.fsencode(strdirectory)
for file in sorted(os.listdir(directory)):
	
	print("\n===========================================================")
	
	# Program Timing (Time taken to get through one file)	
	fstartTime = datetime.now()
	
	file_itrn += 1
	
	# Setting filename for dataload
	filename = os.fsdecode(file)
	print(str(filename), '\n')

	infile = strdirectory + str(filename) # Specifies data file
	print('netCDF file:')
	print(infile, '\n')
	data = nc.Dataset(infile)"""
	
	# ~ """=============================================================="""
	# ~ """======================Download Variables======================"""
	# ~ """=============================================================="""
"""# Longitude
	data_lon = data.variables['lon'][:]
	# Latitude
	data_lat = data.variables['lat'][:]
	# Altitude
	data_alt = data.variables['alt'][:]
	# Horizontal Line of Sight Wind speed
	data_HLOS = data.variables['HLOS_wind_speed'][:]
	# Rayleigh_Grouping
	RG = data.variables['RG'][:]
	# Time
	rayleigh_times = data.variables['time'][:]
	
	# Initialise arrays
	alts = np.empty((0,1))
	lats = np.empty((0,1))
	lons = np.empty((0,1))
	hlos = np.empty((0,1))
	times = np.empty((0,1))
	
	# Trim dataset to below 12km
	for t in range(len(data_alt)):
		if data_alt[t] < 12000:
			alts = np.append(alts, data_alt[t])
			lats = np.append(lats, data_lat[t])
			if data_lon[t] < 180:
				lons = np.append(lons, data_lon[t])
			if data_lon[t] >= 180:
				lons = np.append(lons, data_lon[t]-360)
			hlos = np.append(hlos, data_HLOS[t])
			times = np.append(times, rayleigh_times[t])
	
	# Building arrays
	for t in range(len(hlos)):
		val = hlos[t]
		# Find lat/lon point that satellite is closest to
		lat = min(range(len(lat_grid)), \
			key=lambda i: abs(lat_grid[i]-lats[t]))
		lon = min(range(len(lon_grid)), \
			key=lambda i: abs(lon_grid[i]-lons[t]))
				
		# Add hlos value
		idx = int(data_itrn[lat][lon])
		# ~ data_values[lat][lon][idx] += val
		data_values[lat][lon][0] = np.append(data_values[lat][lon][0], val)
		data_itrn[lat][lon] += 1
		# ~ print(data_values[lat][lon])

# ~ print(data_values[31][31][0])

for i in range(len(data_values)):
	for j in range(len(data_values[i])):
		data_vars[i][j] = np.nanvar(data_values[i][j][0])

# ~ np.set_printoptions(threshold=sys.maxsize)
# ~ for i in range(len(data_values)):
	# ~ print(i)
	# ~ print(data_values[i][34])

# ~ print(data_vars)"""
"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')
# Initialise figure
fig = plt.figure()
ax = fig.add_subplot(111)

mydata = np.load("july19_sub12km_var.npy", allow_pickle=True)

# Plotting data
map = Basemap(projection='spstere',boundinglat=-30,lon_0=90,\
				resolution='l', ax=ax)
map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
z = mydata
print(z)
lons, lats = np.meshgrid(lon_grid, lat_grid)
x, y = map(lons,lats)
cs = map.scatter(x, y, z, cmap='RdBu_r', alpha=1, zorder=3, latlon=True, ax=ax)

map.drawcoastlines(linewidth=0.25, color='#666666')
meridians = np.linspace(-180, 135, 8)
parallels = np.linspace(-80, -30, 6)
map.drawmeridians(meridians, linewidth=0.3)
map.drawparallels(parallels, linewidth=0.3)
# ~ ax.set_xticks(meridians)
# ~ ax.set_yticks(parallels)

# ~ np.save("july19_sub12km_var.npy", data_vars, allow_pickle=True)

# ~ fig.colorbar(cs, cmap='RdBu_r', orientation='horizontal',
			# ~ label='HLOS Rayleigh Wind Speed variance / (ms$^{-1}$)$^{2}$')

# Set figure title
str_plt_title = 'AEOLUS HLOS Rayleigh Wind Variance / (ms$^{-1}$)$^2$'
plt.title(str_plt_title)

# Saving figure
os.chdir('variances')
plt.savefig('july19_2.5_stereo.png',dpi=300)
os.chdir('..')
