#!/usr/bin/env python3
"""
Aeolus data load for Corwin's QBO .mat files
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Testing, still with pcolormesh only---------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .mat files converted from .nc files from the Aeolus database and 
produces plots
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
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.io import savemat, loadmat
from itertools import groupby
from mpl_toolkits.basemap import Basemap

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Change current working directory to parent directory
os.chdir('..')

matfile = '/home/tpb38/PhD/Bath/QBO/aeolus_EQ_newbins.mat'
matData = loadmat(matfile)['Data']
matLat = loadmat(matfile)['Lat']
matLon = loadmat(matfile)['Lon']
matTime = loadmat(matfile)['Time']
matZ = loadmat(matfile)['Z']
# ~ print(matData)
np.set_printoptions(threshold=sys.maxsize)
print('Lats: ', matLat)
print('Lons: ', matLon)
print('Times: ', matTime)
print('Zs: ', matZ)

"""==========================Binning to 2km grid============================="""
newData = np.copy(matData)
newData.fill(0)
for j in range(len(matTime[0])):
	for k in range(len(matZ[0])):
		a = k-2
		b = k+3
		if k == 0 | k == 1:
			a = 0
		if k == 87 | k == 88:
			b = 89
		newData[j,:,:,k] = np.nanmean(matData[j,:,:,a:b], axis = 2)
		
print(newData[0,:,:,15])
matData = np.copy(newData)

# Testing
"""
mean = np.mean(matData[0,:,:,15:20], axis = 2)
array = matData[0,:,:,15:20]
print('Mean: ', mean)
print('Array: ', array)
print('Mean Shape; ', np.shape(mean))
print('Array Shape: ', np.shape(array))"""

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')
# Initialise figure
fig = plt.figure()
ax = fig.add_subplot(111)
gridspec.GridSpec(4,4) # Initialise gridspec

"""===================================AX1===================================="""

# Plot 1
ax1 = plt.subplot2grid((4,4), (2,0), colspan=4, rowspan=1)

# Plotting data
map = Basemap(projection='cyl',llcrnrlat=-35,urcrnrlat=35,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax1)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
z = matData
print(np.shape(matData))
print(len(matTime[0]))
print(np.shape(matData[0,:,:,3]))

z = matData[80,:,:,50]
strheight = '17 +/- 1 km'
z = np.transpose(z)/100
lons, lats = np.meshgrid(matLon[0], matLat[0])
x, y = map(lons,lats)
print(np.shape(z))
print(np.shape(x))
print(np.shape(y))
cs1 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-40, vmax=40, alpha=1, zorder=3)

map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
ax1.set_xticks(meridians)
ax1.set_yticks(parallels)
ax1.set_title(strheight)

"""===================================AX2===================================="""
# Plot 2
ax2 = plt.subplot2grid((4,4), (1,0), colspan=4, rowspan=1)

# Plotting data
map = Basemap(projection='cyl',llcrnrlat=-35,urcrnrlat=35,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax2)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
z = matData
print(np.shape(matData))
print(len(matTime[0]))
print(np.shape(matData[0,:,:,3]))

z = matData[80,:,:,58]
strheight = '19 +/- 1 km'
z = np.transpose(z)/100
lons, lats = np.meshgrid(matLon[0], matLat[0])
x, y = map(lons,lats)
print(np.shape(z))
print(np.shape(x))
print(np.shape(y))
cs2 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-40, vmax=40, alpha=1, zorder=3)

map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
# ~ ax2.set_xticks(meridians)
ax2.set_yticks(parallels)
ax2.set_title(strheight)

"""===================================AX3===================================="""
# Plot 3
ax3 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=1)

# Plotting data
map = Basemap(projection='cyl',llcrnrlat=-35,urcrnrlat=35,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax3)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
z = matData
print(np.shape(matData))
print(len(matTime[0]))
print(np.shape(matData[0,:,:,3]))

z = matData[80,:,:,62] #62 = 21km
strheight = '21 +/- 1 km'
z = np.transpose(z)/100
lons, lats = np.meshgrid(matLon[0], matLat[0])
x, y = map(lons,lats)
print(np.shape(z))
print(np.shape(x))
print(np.shape(y))
cs3 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-40, vmax=40, alpha=1, zorder=3)

map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
# ~ ax2.set_xticks(meridians)
ax3.set_yticks(parallels)
ax3.set_title(strheight)

# ==============================================================================

cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.05])
fig.colorbar(cs1, cmap='RdYlBu_r', orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms$^{-1}$', cax=cbar_ax)

# Set figure title
# ~ str_plt_title = 'AEOLUS QBO Plot'
# ~ plt.title(str_plt_title, y=15.5)

fig.tight_layout(pad=0.5)

# Saving figure
os.chdir('CW_QBOplots')
plt.savefig('CW_QBOplot2.png',dpi=300)

"""
# Test Figure
for i in range(len(matZ[0])):
	fignum = str(matZ[0][i])
	if (float(fignum)*2)%1 == 0:
		fignum += '0'
	figname = 'test_' + fignum + 'km.png'
	plt.close()
	fig2 = plt.figure()
	ax = fig.add_subplot(111)
	x2 = matTime[0]
	y2 = matData[:,20,20,i]
	plt.plot(x2, y2)
	plt.savefig(figname,dpi=300)"""

os.chdir('..')
