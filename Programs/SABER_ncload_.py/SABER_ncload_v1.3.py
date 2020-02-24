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
import os
os.putenv('CODA_DEFINITION',
'/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from itertools import groupby

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
strfile = 'SABER_Temp_O3_H2O_June2019_v2.0.nc'
infile = strdirectory + strfile
data = nc.Dataset(infile)

# Toggles
xaxis_type = 'time' # 'track'

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
		ax1 = plt.subplot2grid((1,1), (0,0))
		cs = plt.pcolormesh(x, y, z, cmap='plasma', vmin=160, vmax=300) # vmin = 120, vmax = 300
		fig.colorbar(cs, cmap='plasma', orientation='horizontal',
					label='Kinetic Temperature / K')
		ax1.set_ylabel('Altitude / km')
		if xaxis_type == 'time':
			fixmydateaxis(ax1, xaxis) # Uses my own function 'fixmydateaxis'
			ax1.set_xlabel('Time of day / HH:MM')
		elif xaxis_type == 'track':
			ax1.set_xlabel('Along track distance / km')
		str_plt_title = 'SABER Kinetic Temperature Profiles'
		str_plt_title += '\n' + 'Date: ' + str_date + ' ' + str(complete_boxes+1)
		plt.title(str_plt_title)
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
