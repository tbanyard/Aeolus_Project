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
---v1.5---Looping over height levels and days---------------------------
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
import matplotlib.ticker as ticker
import os
import errno
import warnings
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

# Edit warnings
warnings.filterwarnings("error", message="No contour levels were found within the data range.", category=UserWarning)

# Change current working directory to parent directory
os.chdir('..')

# Loading data and printing dimensions
matfile = '/home/tpb38/PhD/Bath/QBO/aeolus_EQ_newbins.mat'
matfile = '/home/tpb38/PhD/Bath/QBO/aeolus_EQ_2km_p25rolling_v3_10days.mat'
matData = loadmat(matfile)['Data']
matLat = loadmat(matfile)['Lat']
matLon = loadmat(matfile)['Lon']
matTime = loadmat(matfile)['Time']
matZ = loadmat(matfile)['Z']
np.set_printoptions(threshold=sys.maxsize)
print('Lats: ', matLat)
print('Lons: ', matLon)
print('Times: ', matTime)
print('Zs: ', matZ)

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
# Get into the correct directory for plotting
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')

# Set colormaps
qbocmap = LinearSegmentedColormap('QBOcustomcmap', segmentdata=customcolormaps('QBOcmap2'), N=265)
grayhatchescmap = LinearSegmentedColormap('Grayhatchescmap', segmentdata=customcolormaps('grayhatches'), N=265)
grayhatchescmap_r = grayhatchescmap.reversed()
blackhatchescmap = LinearSegmentedColormap('Blackhatchescmap', segmentdata=customcolormaps('blackhatches'), N=265)
blackhatchescmap_r = blackhatchescmap.reversed()

# Create datetime array using basedate and matTime array
basedate = datetime(2000, 1, 1)
dtTime = np.full(len(matTime[0]), basedate)
for t in range(len(matTime[0])):
	dtTime[t] += timedelta(days=int(matTime[0][t]))
mplTime = dates.date2num(dtTime)

for t in range(len(matTime[0])):
	# Initialise figure
	fig = plt.figure()
	outer = gridspec.GridSpec(4,4) # Initialise gridspec
	# ~ subplotparams = gridspec.GridSpec.get_subplot_params(outer)
	
	"""========Looping through each of the three subplots========"""
	ax1 = plt.subplot2grid((4,4), (2,0), colspan=4, rowspan=1) # 19km [58]
	ax2 = plt.subplot2grid((4,4), (1,0), colspan=4, rowspan=1) # 21km [66]
	ax3 = plt.subplot2grid((4,4), (0,0), colspan=4, rowspan=1) # 23km [74]

	for ax in [ax1, ax2, ax3]:
		# Create basemap
		map = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
				llcrnrlon=-180,urcrnrlon=180,resolution='l', lon_0=0, ax=ax)
		
		# Set height elements
		if ax == ax1:
			alt = 19
			altidx = 58
		elif ax == ax2:
			alt = 21
			altidx = 66
		elif ax == ax3:
			alt = 23
			altidx = 74
		
		# Initialise z array and apply filters
		z = matData[t,:,1:-2,altidx]
		z[-1,:] = z[0,:] # Wrap around the dateline
		strheight = str(alt) + ' $\pm$ 1 km'
		z = np.transpose(z)/100
		# ~ print(strheight)
		# ~ print(z)
		
		# Set NaNs to mean and create binary array of NaNs.
		isnanarray = np.isnan(z)
		hatcharray = np.copy(z)
		mean = np.nanmean(z)
		for lonidx in range(len(isnanarray)):
			for latidx in range(len(isnanarray[lonidx])):
				if isnanarray[lonidx][latidx] == True:
					z[lonidx][latidx] = mean
					hatcharray[lonidx][latidx] = 1
				else:
					hatcharray[lonidx][latidx] = 0
					
		if False in isnanarray:
					
			# ~ z = ndimage.uniform_filter(z, size=(3,2), mode = 'reflect')
			z = savgol_filter(z, 9, 2, axis = 0) # Meridional S-G filter
			# ~ z = savgol_filter(z, 5, 2, axis = 1) # Zonal S-G filter
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
			# ~ cs3 = ax.pcolormesh(x, y, z, cmap='RdYlBu_r', shading = 'bilinear', vmin=-30, vmax=30, alpha=1, zorder=3)
			# ~ cs3 = ax.contourf(x, y, z, cmap='RdYlBu_r', levels = 7, vmin=-30, vmax=30, alpha=1, zorder=3)
			cs3 = ax.contourf(xi2, yi2, zi2, cmap=qbocmap, levels=[-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50], vmin=-30, vmax=30, alpha=1, zorder=3)
			# ~ cs4 = ax.contour(x, y, z, colors = 'k', linestyles = 'dotted', linewidths = 0.3, interpolation = 'bicubic', vmin=-30, vmax=30, alpha=1, zorder=3)
			# ~ cs4 = ax.contour(xi2, yi2, zi2, linestyles = 'solid', levels = 11, linewidths = 0.3, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
			try:
				cs4 = ax.contour(xi2, yi2, zi2, levels=[-50,-40,-30,-20,-10,10,20,30,40,50], linewidths = 0.4, colors='k', vmin=-30, vmax=30, alpha=1, zorder=3)
			except:
				print("Unable to plot contours")
			else:
				cls = ax.clabel(cs4, [-40,-30,-20,-10,10,20,30,40], fmt = '%i', fontsize=4, inline_spacing=1)
			try:		
				cs5 = ax.contour(xi2, yi2, zi2, levels=[0], linewidths = 0.8, colors = 'k', alpha=1, zorder=4)
			except:
				print("Unable to plot zero contour")
			else:
				cls2 = ax.clabel(cs5, [0], fmt = '%i', fontsize=4, inline_spacing=1)
			mask = ax.pcolor(x, y, hatcharray, cmap = grayhatchescmap, facecolor = 'r', zorder = 5, edgecolor = 'none')
			# ~ hatchescontour = ax.contour(x, y, hatcharray, levels = [0.75], interpolation = 0, colors = 'k', linestyles = 'dashed', linewidths = 1, alpha = 1, zorder = 5)
			# ~ hatches = ax.contourf(x, y, hatcharray, levels = 2, cmap = blackhatchescmap, alpha = 0.25, zorder = 5)

		# Complete basemap
		map.drawcoastlines(linewidth=0.25, color='#666666', zorder=4)
		meridians = np.linspace(-180, 180, 9)
		parallels = np.linspace(-30, 30, 7)
		map.drawmeridians(meridians, linewidth=0.3, zorder=5)
		map.drawparallels(parallels, linewidth=0.3, zorder=5)
		ax.tick_params(labelsize=8)
		
		# The following section of code labels the meridian ticks correctly
		if ax == ax1:
			ax.set_xticks(meridians)
			f = lambda x,pos: "{}$^{{\circ}}${EW}".format(int(abs(x)), EW=lonhemi(x)) # Double braces {{}} to treat literally
			xformatter = ticker.FuncFormatter(f)
			ax.xaxis.set_major_formatter(xformatter)
		
		# The following section of code labels the parallel ticks correctly
		ax.set_yticks(parallels)
		f = lambda x,pos: "{}$^{{\circ}}${NS}".format(int(abs(x)), NS=lathemi(x)) # Double braces {{}} to treat literally
		yformatter = ticker.FuncFormatter(f)
		ax.yaxis.set_major_formatter(yformatter)		
		ax.set_title(strheight, {'fontsize': 11})

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

	# Find date for title and savename
	### Title
	curr_date = dtTime[t].strftime('%d-%b-%Y')
	next_date = dtTime[t] + timedelta(days=9)
	next_date = next_date.strftime('%d-%b-%Y')
	
	### Savename
	curr_date_yyyymmdd = dtTime[t].strftime('%Y-%m-%d')

	# Set figure title
	str_plt_title = curr_date + ' to ' + next_date
	plt.title(str_plt_title, {'fontsize': 16}, y=32)

	# Use tight layout
	fig.tight_layout(pad=0.5)
	plt.subplots_adjust(left = 0.125, bottom = -0.03, right=0.91, top=0.875)

	# Saving figure
	os.chdir('CW_QBOplots')
	str_fig_name = curr_date_yyyymmdd + '_Aeolus_QBO_map.png'
	plt.savefig(str_fig_name,dpi=300)
	print(str_fig_name, " saved.")
	os.chdir('..')
	plt.close()

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
