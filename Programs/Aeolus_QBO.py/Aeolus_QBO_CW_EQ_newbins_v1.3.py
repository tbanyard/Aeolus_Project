#!/usr/bin/env python3
"""
Aeolus data load for Corwin's QBO .mat files
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Testing, still with pcolormesh only---------------------------
---v1.2---Improving plot quality----------------------------------------
---v1.3---Creating 2D Test figure to see data gaps----------------------
----------[DEPRECATED]-There_is_a_newer_version_of_this_file------------
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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import colorbar

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Change current working directory to parent directory
os.chdir('..')

matfile = '/home/tpb38/PhD/Bath/QBO/aeolus_EQ_newbins.mat'
matfile = '/home/tpb38/PhD/Bath/QBO/aeolus_EQ_2km_p25rolling.mat'
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
"""N.B. Only use the code below for Corwin's original file, aeolus_EQ_newbins"""
"""=========================================================================="""
"""newData = np.copy(matData)
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
matData = np.copy(newData)"""

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
cs1 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)

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
cs2 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)

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

z = matData[80,:,:,66] #66 = 21km
strheight = '21 +/- 1 km'
z = np.transpose(z)/100
lons, lats = np.meshgrid(matLon[0], matLat[0])
x, y = map(lons,lats)
print(np.shape(z))
print(np.shape(x))
print(np.shape(y))
x_lims = [-180, 180]
y_lims = [-35, 35]
im_interp = 'none'
# ~ cs3 = plt.imshow(z, aspect='auto', cmap='RdYlBu_r', extent=[x_lims[0],
				# ~ x_lims[1], y_lims[0], y_lims[1]], vmin=-30, vmax=30,
				# ~ interpolation=im_interp)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow
cs3 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)

map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
# ~ ax2.set_xticks(meridians)
ax3.set_yticks(parallels)
ax3.set_title(strheight)

# ==============================================================================

cbar_ax = fig.add_axes([0.1, 0.15, 0.76, 0.05])
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
	
"""=========================================================================="""
"""=============================2D Test Figure==============================="""
"""=========================================================================="""
plt.close()
fig2 = plt.figure()
ax = fig2.add_subplot(111)

# Create datetime array using basedate and matTime array
basedate = datetime(2000, 1, 1)
dtTime = np.full(len(matTime[0]), basedate)
for t in range(len(matTime[0])):
	dtTime[t] += timedelta(days=int(matTime[0][t]))
mplTime = dates.date2num(dtTime)

x = matTime[0]
x = mplTime
y = matZ[0][7:]
x, y = np.meshgrid(x, y)

alpha = 1
alpha_itrn = 0
for fractionrequired in [0.1, 0.25, 0.5, 0.75, 0.9]:
	alpha_itrn += 1
	alpha -= 0.35 / alpha_itrn
	# ~ fractionrequired = 0.1 # Between 0-1 for 0-100%
	gapsbinaryarray = np.empty((86, 89), dtype=int)
	for j in range(len(matTime[0])):
		for k in range(len(matZ[0])):
			isnanarray = np.isnan(matData[j,:,:,k])
			unique, counts = np.unique(isnanarray, return_counts=True)
			countdict = dict(zip(unique, counts))
			try:
				countdict[True]
			except:
				fractionempty = 1
			else:
				fractionempty = countdict[True] / np.size(isnanarray)
			
			if 1-fractionempty<fractionrequired:
				gapsbinaryarray[j, k] = 0
			else:
				gapsbinaryarray[j, k] = 1
	z = np.transpose(gapsbinaryarray[:,7:])
	newcmp = LinearSegmentedColormap('testCmap', segmentdata=customcolormaps('wg'), N=265)
	cs2d = plt.pcolormesh(x, y, z, cmap=newcmp, edgecolor = 'white', linewidth = 0.25, alpha=alpha, zorder=3)

ax.set_yticks(np.linspace(6,24,10))
ax.xaxis_date() # Initialises date axis
date_form = dates.DateFormatter('%b\n%Y') # Sets date format
ax.xaxis.set_major_formatter(date_form)
for label in ax.get_xticklabels(): # Center align x tick labels
    label.set_ha("center")
ax.set_xlabel('Date')
ax.set_ylabel('Altitude / km')
plt.ylim(5.75, 24)
ax.yaxis.set_minor_locator(plt.MaxNLocator(19))
plt.subplots_adjust(bottom = 0.15, right=0.85, top=0.9)
plt.title('Aeolus Data Coverage for the 35$^{\circ}$N/S Equatorial Band')
# ~ fig.colorbar(cs2d, ax=ax, orientation='vertical')

plottestcmp = LinearSegmentedColormap('testCmap', segmentdata=customcolormaps('2dplottest'), N=265)
cbar_ax = fig2.add_axes([0.87, 0.15, 0.02, 0.75], zorder=5)
colorbar.ColorbarBase(cbar_ax, cmap = plottestcmp, orientation='vertical',
	label='Coverage / %', boundaries = [0,10,25,50,75,90,100])
# ~ fig2.colorbar(cs2d, cmap='RdBu', orientation='vertical',
	# ~ label='Coverage / %', cax=cbar_ax, boundaries = [0,10,25,50,75,90,100])
for axis in ['top','bottom','left','right']: # Set axes thickness
	ax.spines[axis].set_linewidth(1.5)
ax.grid(False)
plt.savefig('2dplottest.png',dpi=300)

os.chdir('..')
