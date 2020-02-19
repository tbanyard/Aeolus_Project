#!/usr/bin/env python3
"""
Plot Aeolus data segment as pcolormesh or imagesc/imshow from npy file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Adding subplot with lat/lon grid showing satellite track------
---v1.2---Editing for running on ubpc-2027: 19.02.2020_work-------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Plots the Aeolus data from an npy file for an orbital segment as either
a pcolormesh plot or an imagesc/imshow plot.
========================================================================
"""

"""-----------------------------------------------------------------------------
----------------------------------Imports---------------------------------------
-----------------------------------------------------------------------------"""

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
from itertools import groupby
from mpl_toolkits.basemap import Basemap

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import ncload

# ~ enterdirectory("test")
os.chdir('..')
os.chdir('..')
os.chdir('Plots')
print(os.getcwd())

"""-----------------------------------------------------------------------------
----------------------------------Program---------------------------------------
-----------------------------------------------------------------------------"""
# Choose pcolor or imshow
pc_or_im = 'im'

# Load data
data = np.load('data.npy', allow_pickle=True)
date_time = np.load('date_time.npy', allow_pickle=True)
data_track = np.load('data_track.npy', allow_pickle=True)
rayleigh_times = np.load('rayleigh_times.npy', allow_pickle=True)
x, y, z = data[0], data[1], data[2]
z = np.array(z, dtype='float64')
data_lat, data_lon = data_track[0], data_track[1]

# Initialise figure
fig = plt.figure()
gridspec.GridSpec(5,5) # Initialise gridspec

# Main plot
ax1 = plt.subplot2grid((5,5), (0,0), colspan=3, rowspan=4)

# Choose pc or im
if pc_or_im == 'pc':
	cs = plt.pcolormesh(x, y/1000, z, cmap='RdBu_r', vmin=-200, vmax=200)
	fixmydateaxis(ax1, date_time) # Uses my own function 'fixmydateaxis'

elif pc_or_im == 'im':
	# Plot imshow plot
	# Flatten x_meshgrid and convert from datetime format to mpl_compatible nums
	x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
	x_lims = dates.date2num(x_lims)
	# Y limits
	y_lims = [20.5, -0.5]
	fixnanswithmean(z) # Uses my own function 'fixnanswithmean'
	# Plots using imshow
	cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
		x_lims[1], y_lims[0], y_lims[1]], vmin=-200, vmax=200,
		interpolation='spline36')
	ax1.xaxis_date() # Initialises date axis
	date_form = dates.DateFormatter('%H:%M') # Sets date format
	ax1.xaxis.set_major_formatter(date_form)
	plt.gca().invert_yaxis() # Invert axis for imshow

# z array check
print(z)

# Fix axes

### X axis
ax1.set_xlabel('Time')
# Ensure the number of date ticks is sensible
timerange = date_time[-1] - date_time[0]
date_intvl = int(np.ceil(timerange.seconds/(5*60)))
ax1.xaxis.set_major_locator(plt.MaxNLocator(5)) # Maximum number of date ticks
ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=date_intvl))
# ~ ax1.grid(color='k', linestyle = 'dashed', linewidth = 0.25, axis='x')

### Y axis
ax1.set_ylabel('Altitude / km')
ax1.set_yticks(np.arange(len(y)))
ax1.yaxis.set_major_locator(plt.MaxNLocator(11))
ax1.yaxis.set_minor_locator(plt.MaxNLocator(21))
# ~ ax1.tick_params(axis='y', which='minor', left=True) # Minor ticks
ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='y' ,which='both')

# Satellite track plot
ax2 = plt.subplot2grid((5,5), (0,3), colspan=2, rowspan=4)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
map = Basemap(projection='cyl',llcrnrlat=-80,urcrnrlat=-40,\
            llcrnrlon=-80,urcrnrlon=-40,resolution='i', ax=ax2)
map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
map.drawcoastlines(linewidth=0.25, color='#666666')
map.drawmeridians([-70, -60, -50], linewidth=0.3)
map.drawparallels([-70, -60, -50], linewidth=0.3)
map.scatter(data_lon - 360, data_lat, marker = 'x', color = 'red', s=0.3, zorder=2)
map.plot(data_lon - 360, data_lat, color = 'red', linewidth = '0.5')
midpointlonindex = int(np.floor(len(data_lon)/2))
midpointlatindex = int(np.floor(len(data_lat)/2))
# ~ map.quiver(data_lon[midpointlonindex]-360, data_lat[midpointlatindex], data_lon[midpointlonindex+100] - data_lon[midpointlonindex-100], data_lat[midpointlatindex+100] - data_lat[midpointlatindex-100], color='black', zorder=3, headaxislength=0, headlength=0, pivot='middle')
# ~ map.quiver(data_lon[-1]-360, data_lat[-1], data_lon[-2]-data_lon[-1], data_lat[-2]-data_lat[-1], color='black', zorder=4)

# Add lines at time intervals corresponding to date_time ticks
date_ticks = dates.num2date(ax1.xaxis.get_ticklocs())

for i in range(len(date_ticks)):
	# ~ print(date_ticks[i].minute)
	# Co-locate between ticks and rayleigh_times array
	nearest_date = find_nearest(rayleigh_times, date_ticks[i])
	idx = (np.where(dates.date2num(rayleigh_times) == nearest_date)[0])[0]
	str_time = date_ticks[i].strftime('%H:%M') # Sets string form of each tick
	# Location of each line
	lon_loc = data_lon[idx]-360
	lat_loc = data_lat[idx]
	# ~ map.scatter([lon_loc], [lat_loc], s=50, marker='+', zorder=111, color='black')
	# Plots line on map and annotates with text
	try: 
		data_lat[idx+30] == False
		data_lat[idx-30] == False
		data_lon[idx+30] == False
		data_lon[idx-30] == False
	except:
		continue
	map.quiver(lon_loc, lat_loc, data_lat[idx-30] - data_lat[idx+30], data_lon[idx+30] - data_lon[idx-30], angles='xy', color='black', zorder=3, headaxislength=0, headlength=0, pivot='middle', units='xy')
	ax2.annotate(str_time, (lon_loc+1, lat_loc+1), fontsize=6, zorder=5)
	
	# ~ a = ax1.xaxis.get_ticklabels()
	# ~ print(text.Text(agg_filter=a))

# ~ map.plot([-70, -65], [-70, -67], color='k')

# Fix axes
ax2.set_xticks([-80, -70, -60, -50, -40])
ax2.set_yticks([-80, -70, -60, -50, -40])
ax2.set_xlabel('Longitude / deg')
ax2.set_ylabel('Latitude / deg')
ax2.set_aspect('auto') # Stretch map to fill subplot
for axis in ['top','bottom','left','right']: # Set axes thickness
	ax1.spines[axis].set_linewidth(0.75)
	ax2.spines[axis].set_linewidth(0.75)

# Add colorbar to figure
fig.subplots_adjust(bottom=0.2, right=0.88, left=0.12)
cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.05])
fig.colorbar(cs, cmap='RdBu_r', orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms$^{-1}$', cax=cbar_ax)

# Set figure title
plt.title('Aeolus Orbit HLOS Rayleigh Wind Cross-section', y=15)

# Save figure
os.chdir('/home/tpb38/PhD/Bath/Aeolus_Project/Plots/dump')
plt.savefig('testbed.png',dpi=300)
