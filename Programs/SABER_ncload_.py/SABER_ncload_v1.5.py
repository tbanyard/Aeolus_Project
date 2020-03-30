#!/usr/bin/env python3
"""
Testbed
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Splitting by time, into sections------------------------------
---v1.2---40-80 Andes box: Preparing code for this experiment.----------
----------Solved issue with the time formatting and overlapping times.--
---v1.3---Complete 40-80 Andes box.-------------------------------------
---v1.4---Satellite track subplot.--------------------------------------
---v2.0---23.03.20-Updates: Applying S-G Filtering----------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Testbed
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import os
os.putenv('CODA_DEFINITION',
'/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import scipy.ndimage as ndimage
from itertools import groupby
from mpl_toolkits.basemap import Basemap

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import ncload

# ~ enterdirectory("test")
# Change current working directory to parent directory
os.chdir('..')
print(os.getcwd())

"""Find directory and read netCDF data"""
strdirectory = '/home/tpb38/PhD/Bath/SABER/'
strfile = 'SABER_Temp_O3_H2O_July2019_v2.0.nc'
infile = strdirectory + strfile
data = nc.Dataset(infile)

# Toggles
xaxis_type = 'time' # 'track'
pc_or_im = 'im'
im_interp = 'none'

# Print full arrays without truncation
np.set_printoptions(threshold=sys.maxsize)

"""Find out which indices correspond to each date in the month"""
# Date
data_all_dates = data.variables['date'][:] # Load entire dataset of dates
date = data_all_dates[0] # First date in yyyydoy format
day_starts = np.empty((0,1), dtype=int) # Initialise numpy array of day_starts

# Initialise numpy array of dates in datetime format
str_date = str(date)
YYYY = int(str_date[0:4])
doy = int(str_date[4:7])
dtdate = doy_to_yyyymmdd(YYYY, doy)
dtdates = np.array(dtdate)

# Create an array of the elements at which each day begins
for elmnt in range(1, len(data_all_dates)):
	curr_date = data_all_dates[elmnt]
	if curr_date != date:
		day_starts = np.append(day_starts, elmnt)
		date = curr_date
		
		# Build dtdates array
		str_date = str(date)
		YYYY = int(str_date[0:4])
		doy = int(str_date[4:7])
		dtdate = doy_to_yyyymmdd(YYYY, doy)
		dtdates = np.append(dtdates, dtdate)
	
# Date loop	
end = 0
for day in range(len(dtdates)):
	# Build units for datetime objects in data_time array
	str_date = str(dtdates[day])[0:10]
	str_units = 'milliseconds since ' + str_date + ' 00:00:00'
	start = end
	if day == len(day_starts):
		end = len(data_all_dates)
	else:
		end = day_starts[day]
	
	"""===================================================================="""
	"""=========================Download Variables========================="""
	"""===================================================================="""
	# Longitude
	data_lon = data.variables['tplongitude'][start:end]
	# Latitude
	data_lat = data.variables['tplatitude'][start:end]
	# Altitude
	data_alt = data.variables['tpaltitude'][start:end]
	# Geopotential Altitude
	data_gpalt = data.variables['tpgpaltitude'][start:end]
	# Pressure
	data_press = data.variables['pressure'][start:end]
	# Density
	data_dens = data.variables['density'][start:end]
	# Kinetic Temperature
	data_ktemp = data.variables['ktemp'][start:end]
	# Date
	data_date = data.variables['date'][start:end]
	# Time
	data_time = data.variables['time'][start:end,350]
	# Fix times which overlap past midnight
	overlap_times = np.where(data_time[:-10] > 86400000)
	for time_elmnt in overlap_times[0]:
		data_time[time_elmnt] -= 86400000
	# Converted time
	data_time = nc.num2date(data_time,\
	calendar = 'standard', units = str_units)

	# Print full arrays without truncation
	np.set_printoptions(threshold=sys.maxsize)

	"""=============================================================="""
	"""=====Test to see if orbit is sufficiently within Andes box===="""
	"""=============================================================="""
	mnopib = 6 # minimum_number_of_profiles_in_box
	# ~ print(np.where(data_lat<-80, 0, (np.where(data_lat>-40, 0, 1))))
	# Find where the satellite is within the Andes box
	box = np.where(data_lat[:,350]<-80, 0, (np.where(data_lat[:,350]>-40, 0,
		(np.where(data_lon[:,350]>320, 0, (np.where(data_lon[:,350]<280, 0, 1)))))))
	diffs = np.diff(box) # Array of diffs for box
	# Grouping the differences between the elements in diffs
	# (2nd derivative)
	grouped_diffs = [(k, sum(1 for i in g)) for k,g in groupby(diffs)]
	# Returns:[(0, 6206),(1, 1),...,(0, 1748),(-1, 1),...,(0, 8617)]
	# ~ print(box)
	# ~ print(grouped_diffs)
	# Finding the start and end elements of the desired section
	itrn = 0
	start_elmnt = 0
	end_elmnt = 0
	complete_boxes = 0
	for u in grouped_diffs:
		if u[0] != 0: # Bypass 1's and -1's
			itrn += u[1]
		elif u[0] == 0:
			# Is this section the Andes box?
			if data_lat[:,350][itrn] < -40 and data_lat[:,350][itrn] > -80 and \
			box[itrn] == 1:
				# Are there enough profiles in the box?
				if u[1] > mnopib:
					print(u[1])
					if start_elmnt == 0:
						# First profile in box
						start_elmnt = itrn
						# Last profile in box
						end_elmnt = itrn + u[1]
					itrn += u[1]
				else:
					itrn += u[1]
			else:
				itrn += u[1]

		if end_elmnt == 0:
			continue
		
		# Amend arrays
		data_lon_new = data_lon[start_elmnt:end_elmnt+1]
		data_lat_new = data_lat[start_elmnt:end_elmnt+1]
		data_alt_new = data_alt[start_elmnt:end_elmnt+1]
		data_gpalt_new = data_gpalt[start_elmnt:end_elmnt+1]
		data_press_new = data_press[start_elmnt:end_elmnt+1]
		data_dens_new = data_dens[start_elmnt:end_elmnt+1]
		data_ktemp_new = data_ktemp[start_elmnt:end_elmnt+1]
		data_time_new = data_time[start_elmnt:end_elmnt+1]
		
		"""================================================================"""
		"""==========================Plotting=============================="""
		"""================================================================"""
		os.chdir('..')
		os.chdir('Plots')
		
		"""==================Track distance along x-axis?=================="""
		if xaxis_type == 'track':
			lats = np.deg2rad(data_lat_new[:,350])
			lons = np.deg2rad(data_lon_new[:,350])
			 # Levels are inverted
			levels = np.linspace(len(data_lat_new[0])-1, 0, len(data_lat_new[0]))
			alts = np.linspace(0, 154, 78)
			track = np.zeros(len(data_lat_new))
			track[0] = 0
			trackdistance = 0
			R = 6371 # Earth's radius in km

			for i in range(1, len(lats)):
				# Haversine formula
				a = (np.sin((lats[i]-lats[i-1])/2) ** 2) + (np.cos(lats[i-1]) * \
					np.cos(lats[i]) * (np.sin((lons[i]-lons[i-1])/2) ** 2))
				c = 2 * np.arcsin(np.sqrt(a))
				d = R * c
				trackdistance += d
				track[i] = trackdistance
			xaxis = track
		
		"""=========================Time along x-axis?==========================="""
		if xaxis_type == 'time':
			alts = np.linspace(0, 154, 78)
			xaxis = data_time_new
		
		"""Plotting"""
		# Initialise arrays
		z = np.zeros((len(alts), len(xaxis)))
		z_itrn = np.zeros((len(alts), len(xaxis)))

		# Placing temperature values into bins of height 2km
		for event in range(len(data_ktemp_new)):
			for level in range(len(data_ktemp_new[0])):
				val = find_nearest(alts, data_alt_new[event][level])
				alt_elmnt = np.where(alts == val)[0][0]
				z[alt_elmnt][event] += data_ktemp_new[event][level]
				z_itrn[alt_elmnt][event] += 1
		
		# Build plot		
		z /= z_itrn
		x, y = np.meshgrid(xaxis, alts)
		x, y, z = x[10:51], y[10:51], z[10:51]
		fig = plt.figure()
		gridspec.GridSpec(5,5) # Initialise gridspec
		
		# Main plot
		ax1 = plt.subplot2grid((5,5), (0,0), colspan=3, rowspan=4)
		
		# Choose pc or im
		if pc_or_im == 'pc':
			cs = plt.pcolormesh(x, y, z, cmap='plasma', vmin=160, vmax=300) # vmin = 120, vmax = 300
		
		if pc_or_im == 'im':
			# Plot imshow plot
			# Flatten x_meshgrid and convert from datetime format to mpl_compatible nums
			x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
			if xaxis_type == 'time':
				x_lims = dates.date2num(x_lims)
			# Y limits
			y_lims = [101, 19]
			fixnanswithmean(z) # Uses my own function 'fixnanswithmean'
			
			z = ndimage.uniform_filter(z, size=(3,7), mode = 'reflect')
			z1 = np.copy(z)
			z2 = savgol_filter(z, 15, 2, axis = 0) # Vertical S-G filter
			# ~ z2 = savgol_filter(z2, 5, 2, axis = 1) # Vertical S-G filter
			z = z1-z2
			
			# Plots using imshow
			cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=-20, vmax=20,
				interpolation=im_interp)
			if xaxis_type == 'time':
				ax1.xaxis_date() # Initialises date axis
				date_form = dates.DateFormatter('%H:%M') # Sets date format
				ax1.xaxis.set_major_formatter(date_form)
			plt.gca().invert_yaxis() # Invert axis for imshow
		
		ax1.set_ylabel('Altitude / km')
		if xaxis_type == 'time':
			fixmydateaxis(ax1, xaxis) # Uses my own function 'fixmydateaxis'
			ax1.set_xlabel('Time of day / HH:MM')
			# Ensure the number of date ticks is sensible
			timerange = data_time_new[-1] - data_time_new[0]
			date_intvl = int(np.ceil(timerange.seconds/(5*60)))
			ax1.xaxis.set_major_locator(plt.MaxNLocator(5)) # Maximum number of date ticks
			ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=date_intvl))
		elif xaxis_type == 'track':
			ax1.set_xlabel('Along track distance / km')
			
		# Satellite track plot
		ax2 = plt.subplot2grid((5,5), (0,3), colspan=2, rowspan=4)	
		
		# =====
		
		ax2.yaxis.set_label_position("right")
		ax2.yaxis.tick_right()
		map = Basemap(projection='cyl',llcrnrlat=-80,urcrnrlat=-40,\
					llcrnrlon=-80,urcrnrlon=-40,resolution='i', ax=ax2)
		map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
		map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
		map.drawcoastlines(linewidth=0.25, color='#666666')
		map.drawmeridians([-70, -60, -50], linewidth=0.3)
		map.drawparallels([-70, -60, -50], linewidth=0.3)
		map.plot(data_lon_new[:,150:] - 360, data_lat_new[:,150:], color = 'red', linewidth = '0.5', alpha=0.05, zorder=1)
		map.scatter(data_lon_new[:,150:] - 360, data_lat_new[:,150:], marker = 'x', color = 'black', s=0.3, alpha=0.05, zorder=2)
		midpointlonindex = int(np.floor(len(data_lon_new)/2))
		midpointlatindex = int(np.floor(len(data_lat_new)/2))
		# ~ map.quiver(data_lon[midpointlonindex]-360, data_lat[midpointlatindex], data_lon[midpointlonindex+100] - data_lon[midpointlonindex-100], data_lat[midpointlatindex+100] - data_lat[midpointlatindex-100], color='black', zorder=3, headaxislength=0, headlength=0, pivot='middle')
		# ~ map.quiver(data_lon[-1]-360, data_lat[-1], data_lon[-2]-data_lon[-1], data_lat[-2]-data_lat[-1], color='black', zorder=4)
		
		if xaxis_type == 'time':
			# Add lines at time intervals corresponding to date_time ticks
			date_ticks = dates.num2date(ax1.xaxis.get_ticklocs())

			for i in range(len(date_ticks)):
				# ~ print(date_ticks[i].minute)
				# Co-locate between ticks and rayleigh_times array
				nearest_date = find_nearest(data_time_new, date_ticks[i])
				idx = (np.where(dates.date2num(data_time_new) == nearest_date)[0])[0]
				str_time = date_ticks[i].strftime('%H:%M') # Sets string form of each tick
				# Location of each line
				lon_loc = data_lon_new[idx,350]-360
				lat_loc = data_lat_new[idx,350]
				# ~ map.scatter([lon_loc], [lat_loc], s=50, marker='+', zorder=111, color='black')
				# Plots line on map and annotates with text
				try: 
					data_lat_new[idx+1,350] == False
					data_lat_new[idx-1,350] == False
					data_lon_new[idx+1,350] == False
					data_lon_new[idx-1,350] == False
				except:
					continue
				# ~ map.quiver(lon_loc, lat_loc, data_lat_new[idx-1,350] - data_lat_new[idx+1,350], data_lon_new[idx+1,350] - data_lon_new[idx-1,350], angles='xy', color='black', zorder=3, headaxislength=0, headlength=0, pivot='middle', units='xy')
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
		
		# =====
			
		# Add colorbar to figure
		fig.subplots_adjust(bottom=0.2, right=0.88, left=0.12)
		cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.05])
		fig.colorbar(cs, cmap='plasma', orientation='horizontal',
			label='Kinetic Temperature Perturbation / K', cax=cbar_ax)		
		
		# Set figure title
		str_plt_title = 'SABER Kinetic Temperature Profiles'
		str_plt_title += '\n' + 'Date: ' + str_date + ' ' + str(complete_boxes+1)
		plt.title(str_plt_title, y=15)
		str_plt_savename = 'sabertest_' + str_date + '_' + str(complete_boxes+1)
		os.chdir('SABER')
		plt.savefig(str_plt_savename,dpi=300)
		print('Figure: ', str_plt_savename, " saved in:", os.getcwd())
		
		# Reset after box completed
		start_elmnt = 0
		end_elmnt = 0
		complete_boxes += 1
		
		plt.close()
		os.chdir('..')
