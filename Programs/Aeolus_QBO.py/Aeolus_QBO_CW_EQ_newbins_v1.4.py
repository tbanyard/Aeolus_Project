#!/usr/bin/env python3
"""
Aeolus data load for Corwin's QBO .mat files
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Testing, still with pcolormesh only---------------------------
---v1.2---Improving plot quality----------------------------------------
---v1.3---Creating 2D Test figure to see data gaps----------------------
---v1.4---Improving quality of maps-------------------------------------
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
import matplotlib.ticker as ticker
import os
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata, interp2d
from scipy.io import savemat, loadmat
from scipy.signal import savgol_filter
import scipy.ndimage as ndimage
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
# Colormap
qbocmap = LinearSegmentedColormap('QBOcustomcmap', segmentdata=customcolormaps('QBOcmap2'), N=265)
grayhatchescmap = LinearSegmentedColormap('Grayhatchescmap', segmentdata=customcolormaps('grayhatches'), N=265)
blackhatchescmap = LinearSegmentedColormap('Blackhatchescmap', segmentdata=customcolormaps('blackhatches'), N=265)

"""===================================AX1===================================="""

# Plot 1
ax1 = plt.subplot2grid((4,4), (2,0), colspan=4, rowspan=1)

# Initialise basemap
map = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax1)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')

# Initialise z array and apply filters
z = matData[80,:,1:-2,50]
z[-1,:] = z[0,:]
strheight = '17 $\pm$ 1 km'
z = np.transpose(z)/100
z = np.nan_to_num(z) # This needs to be changed to the mean or interpolate
z = savgol_filter(z, 9, 2, axis = 0) # Meridional boxcar
# ~ z = savgol_filter(z, 5, 2, axis = 1) # Zonal boxcar
# ~ z = ndimage.gaussian_filter(z, sigma=0.5, order=0) # Have tried sigma=0.6

# 2D interpolation to round edges
f = interp2d(matLon[0], matLat[0][1:-2], z, kind='cubic') # Have tried linear and quintic
xi2 = np.linspace(matLon[0][0], matLon[0][-1], 500)
yi2 = np.linspace(matLat[0][1], matLat[0][-3], 100)
zi2 = f(xi2, yi2)

# Initialise x,y grids
lons, lats = np.meshgrid(matLon[0], matLat[0][1:-2])
x, y = map(lons,lats)

# Plotting Data
# ~ cs3 = plt.imshow(z, aspect='auto', cmap='RdYlBu_r', extent=[x_lims[0],
				# ~ x_lims[1], y_lims[0], y_lims[1]], vmin=-30, vmax=30,
				# ~ interpolation=im_interp)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow
# ~ cs3 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs3 = plt.contourf(x, y, z, cmap='RdYlBu_r', levels = 7, vmin=-30, vmax=30, alpha=1, zorder=3)
cs3 = plt.contourf(xi2, yi2, zi2, cmap=qbocmap, levels=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs4 = plt.contour(x, y, z, colors = 'k', linestyles = 'dotted', linewidths = 0.3, interpolation = 'bicubic', vmin=-30, vmax=30, alpha=1, zorder=3)
cs4 = plt.contour(xi2, yi2, zi2, levels=[-30,-20,-10,10,20,30], linewidths = 0.4, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cs5 = plt.contour(xi2, yi2, zi2, levels=[0], linewidths = 0.8, colors = 'k', alpha=1, zorder=4)
# ~ cs4 = plt.contour(xi2, yi2, zi2, linestyles = 'solid', levels = 11, linewidths = 0.3, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cls = plt.clabel(cs4, [-30,-20,-10,10,20,30], fmt = '%i', fontsize=4, inline_spacing=1)
cls2 = plt.clabel(cs5, [0], fmt = '%i', fontsize=4, inline_spacing=1)

# Complete basemap
map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
ax1.set_xticks(meridians)
ax1.set_yticks(parallels)
ax1.tick_params(labelsize=8)
ax1.set_title(strheight)

# The following section of code labels the meridian ticks correctly
f = lambda x,pos: "{}$^{{\circ}}${EW}".format(int(abs(x)), EW=lonhemi(x)) # Double braces {{}} to treat literally
xformatter = ticker.FuncFormatter(f)
ax1.xaxis.set_major_formatter(xformatter)

# The following section of code labels the parallel ticks correctly
f = lambda x,pos: "{}$^{{\circ}}${NS}".format(int(abs(x)), NS=lathemi(x)) # Double braces {{}} to treat literally
yformatter = ticker.FuncFormatter(f)
ax1.yaxis.set_major_formatter(yformatter)
ax1.set_title(strheight, {'fontsize': 11})

"""===================================AX2===================================="""
# Plot 2
ax2 = plt.subplot2grid((4,4), (1,0), colspan=4, rowspan=1)

# Initialise basemap
map = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax2)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')

# Initialise z array and apply filters
z = matData[80,:,1:-2,58]
z[-1,:] = z[0,:]
strheight = '19 $\pm$ 1 km'
z = np.transpose(z)/100
z = np.nan_to_num(z) # This needs to be changed to the mean or interpolate
z = savgol_filter(z, 9, 2, axis = 0) # Meridional boxcar
# ~ z = savgol_filter(z, 5, 2, axis = 1) # Zonal boxcar
# ~ z = ndimage.gaussian_filter(z, sigma=0.5, order=0) # Have tried sigma=0.6

# 2D interpolation to round edges
f = interp2d(matLon[0], matLat[0][1:-2], z, kind='cubic') # Have tried linear and quintic
xi2 = np.linspace(matLon[0][0], matLon[0][-1], 500)
yi2 = np.linspace(matLat[0][1], matLat[0][-3], 100)
zi2 = f(xi2, yi2)

# Initialise x,y grids
lons, lats = np.meshgrid(matLon[0], matLat[0][1:-2])
x, y = map(lons,lats)

# Plotting Data
# ~ cs3 = plt.imshow(z, aspect='auto', cmap='RdYlBu_r', extent=[x_lims[0],
				# ~ x_lims[1], y_lims[0], y_lims[1]], vmin=-30, vmax=30,
				# ~ interpolation=im_interp)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow
# ~ cs3 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs3 = plt.contourf(x, y, z, cmap='RdYlBu_r', levels = 7, vmin=-30, vmax=30, alpha=1, zorder=3)
cs3 = plt.contourf(xi2, yi2, zi2, cmap=qbocmap, levels=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs4 = plt.contour(x, y, z, colors = 'k', linestyles = 'dotted', linewidths = 0.3, interpolation = 'bicubic', vmin=-30, vmax=30, alpha=1, zorder=3)
cs4 = plt.contour(xi2, yi2, zi2, levels=[-30,-20,-10,10,20,30], linewidths = 0.4, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cs5 = plt.contour(xi2, yi2, zi2, levels=[0], linewidths = 0.8, colors = 'k', alpha=1, zorder=4)
# ~ cs4 = plt.contour(xi2, yi2, zi2, linestyles = 'solid', levels = 11, linewidths = 0.3, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cls = plt.clabel(cs4, [-30,-20,-10,10,20,30], fmt = '%i', fontsize=4, inline_spacing=1)
cls2 = plt.clabel(cs5, [0], fmt = '%i', fontsize=4, inline_spacing=1)

# Complete basemap
map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
# ~ ax2.set_xticks(meridians)
ax2.set_yticks(parallels)
ax2.tick_params(labelsize=8)
ax2.set_title(strheight)

# The following section of code labels the parallel ticks correctly
f = lambda x,pos: "{}$^{{\circ}}${NS}".format(int(abs(x)), NS=lathemi(x)) # Double braces {{}} to treat literally
yformatter = ticker.FuncFormatter(f)
ax2.yaxis.set_major_formatter(yformatter)
ax2.set_title(strheight, {'fontsize': 11})

"""===================================AX3===================================="""
# Plot 3
ax3 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=1)

# Initialise basemap
map = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax3)
# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')

# Initialise z array and apply filters
z = matData[80,:,1:-2,67] #66 = 21km
z[-1,:] = z[0,:]
strheight = '21 $\pm$ 1 km'
z = np.transpose(z)/100

# NAN stuff
isnanarray = np.isnan(z)
hatcharray = np.copy(z)
mean = np.nanmean(z)
print(mean)
print('isnanarray: ', isnanarray)
print(np.shape(isnanarray))
print(len(isnanarray))
for lonidx in range(len(isnanarray)):
	for latidx in range(len(isnanarray[lonidx])):
		if isnanarray[lonidx][latidx] == True:
			z[lonidx][latidx] = mean
			hatcharray[lonidx][latidx] = 1
		else:
			hatcharray[lonidx][latidx] = 0
		
z = np.nan_to_num(z) # This needs to be changed to the mean or interpolate
z = savgol_filter(z, 9, 2, axis = 0) # Meridional boxcar
# ~ z = savgol_filter(z, 5, 2, axis = 1) # Zonal boxcar
# ~ z = ndimage.gaussian_filter(z, sigma=0.5, order=0) # Have tried sigma=0.6

# Attempted griddata interpolation with imshow
# ~ extent = (min(matLon[0]), max(matLon[0]), min(matLat[0]), max(matLat[0]))
# ~ xs,ys = np.mgrid[extent[0]:extent[1]:30j, extent[2]:extent[3]:30j]
# ~ z = griddata((matLon[0],matLat[0]), z, (xs, ys), method='cubic')

# 2D interpolation to round edges
f = interp2d(matLon[0], matLat[0][1:-2], z, kind='cubic') # Have tried linear and quintic
xi2 = np.linspace(matLon[0][0], matLon[0][-1], 500)
yi2 = np.linspace(matLat[0][1], matLat[0][-3], 100)
zi2 = f(xi2, yi2)

# Initialise x,y grids
lons, lats = np.meshgrid(matLon[0], matLat[0][1:-2])
x, y = map(lons,lats)
# Declare x and y limits if imshow is required
# ~ x_lims = [-180, 180]
# ~ y_lims = [-30, 30]
# ~ im_interp = 'none'

# Plotting Data
# ~ cs3 = plt.imshow(z, aspect='auto', cmap='RdYlBu_r', extent=[x_lims[0],
				# ~ x_lims[1], y_lims[0], y_lims[1]], vmin=-30, vmax=30,
				# ~ interpolation=im_interp)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow
# ~ cs3 = plt.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs3 = plt.contourf(x, y, z, cmap='RdYlBu_r', levels = 7, vmin=-30, vmax=30, alpha=1, zorder=3)
cs3 = plt.contourf(xi2, yi2, zi2, cmap=qbocmap, levels=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], vmin=-30, vmax=30, alpha=1, zorder=3)
# ~ cs4 = plt.contour(x, y, z, colors = 'k', linestyles = 'dotted', linewidths = 0.3, interpolation = 'bicubic', vmin=-30, vmax=30, alpha=1, zorder=3)
cs4 = plt.contour(xi2, yi2, zi2, levels=[-30,-20,-10,10,20,30], linewidths = 0.4, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cs5 = plt.contour(xi2, yi2, zi2, levels=[0], linewidths = 0.8, colors = 'k', alpha=1, zorder=4)
# ~ cs4 = plt.contour(xi2, yi2, zi2, linestyles = 'solid', levels = 11, linewidths = 0.3, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
cls = plt.clabel(cs4, [-30,-20,-10,10,20,30], fmt = '%i', fontsize=4, inline_spacing=1)
cls2 = plt.clabel(cs5, [0], fmt = '%i', fontsize=4, inline_spacing=1)
hatchescontour = plt.contour(x, y, hatcharray, levels = [0.75], colors = 'k', linestyles = 'dashed', linewidths = 1, alpha = 1, zorder = 5)
hatches = plt.contourf(x, y, hatcharray, levels = 2, cmap = blackhatchescmap, alpha = 0.25, zorder = 5)

# Complete basemap
map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
meridians = np.linspace(-180, 180, 9)
parallels = np.linspace(-30, 30, 7)
map.drawmeridians(meridians, linewidth=0.3, zorder=5)
map.drawparallels(parallels, linewidth=0.3, zorder=5)
# ~ ax2.set_xticks(meridians)
ax3.set_yticks(parallels)
ax3.tick_params(labelsize=8)

# The following section of code labels the parallel ticks correctly
f = lambda x,pos: "{}$^{{\circ}}${NS}".format(int(abs(x)), NS=lathemi(x)) # Double braces {{}} to treat literally
yformatter = ticker.FuncFormatter(f)
ax3.yaxis.set_major_formatter(yformatter)
ax3.set_title(strheight, {'fontsize': 11})

# ==============================================================================
# Add colorbar to plot
cbar_ax = fig.add_axes([0.125, 0.125, 0.785, 0.025], zorder=5)
colorbar.ColorbarBase(cbar_ax, cmap = qbocmap, orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms$^{-1}$', boundaries = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], extend='both')
for label in cbar_ax.get_xticklabels(): # Center align colorbar tick labels
    label.set_ha("center")


# Make axes thicker
for axis in ['top','bottom','left','right']: # Set axes thickness
	ax1.spines[axis].set_linewidth(1.5)
	ax2.spines[axis].set_linewidth(1.5)
	ax3.spines[axis].set_linewidth(1.5)

# Set figure title
str_plt_title = 'AEOLUS QBO Plot'
plt.title(str_plt_title, {'fontsize': 16}, y=32)

# Use tight layout
fig.tight_layout(pad=0.5)
plt.subplots_adjust(left = 0.125, bottom = -0.03, right=0.91, top=0.875)

# Saving figure
os.chdir('CW_QBOplots')
plt.savefig('CW_QBOplot2.png',dpi=300)

sys.exit(0) # Do not continue onto 2D Test Figure? (Toggle on/off)
	
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
