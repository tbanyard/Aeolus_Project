#!/usr/bin/env python3
"""
Testbed
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Splitting by time, into sections------------------------------
----------[DEPRECATED]-There_is_a_newer_version_of_this_file-------
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
xaxis_type = 'track'

# Print full arrays without truncation
np.set_printoptions(threshold=sys.maxsize)

"""Find out which indices correspond to each date in the month"""
# Date
data_all_dates = data.variables['date'][:]
date = data_all_dates[0]
day_starts = np.empty((0,1), dtype=int)
for elmnt in range(1, len(data_all_dates)):
	curr_date = data_all_dates[elmnt]
	if curr_date != date:
		day_starts = np.append(day_starts, elmnt)
		date = curr_date

# Date loop	
end = 0
for date in range(len(day_starts)+1):
	start = end
	if date == len(day_starts):
		end = len(data_all_dates)
	else:
		end = day_starts[date]
	
	"""=========================Download Variables========================="""
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
	strdate = str(data_date[0])
	# Time
	data_time = data.variables['time'][start:end]
	# Converted time
	# ~ data_time = nc.num2date(data.variables['time'][:],\
	# ~ calendar = 'standard', units = 'milliseconds since 00:00:00')

	# Print full arrays without truncation
	# ~ np.set_printoptions(threshold=sys.maxsize)

	"""======================================================================"""
	"""=============================Plotting================================="""
	"""======================================================================"""
	os.chdir('..')
	os.chdir('Plots')
	
	"""=====================Track distance along x-axis?====================="""
	if xaxis_type == 'track':
		lats = np.deg2rad(data_lat[:,350])
		lons = np.deg2rad(data_lon[:,350])
		 # Levels are inverted
		levels = np.linspace(len(data_lat[0])-1, 0, len(data_lat[0]))
		alts = np.linspace(0, 154, 78)
		track = np.zeros(len(data_lat))
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
	
	# ~ elif xaxis_type == 'time':
	print(data_time)	
	
	"""Plotting"""

	# Initialise arrays
	z = np.zeros((len(alts), len(track)))
	z_itrn = np.zeros((len(alts), len(track)))

	# Placing temperature values into bins of height 2km
	for event in range(len(data_ktemp)):
		for level in range(len(data_ktemp[0])):
			val = find_nearest(alts, data_alt[event][level])
			alt_elmnt = np.where(alts == val)[0][0]
			z[alt_elmnt][event] += data_ktemp[event][level]
			z_itrn[alt_elmnt][event] += 1
			
	z /= z_itrn
	x, y = np.meshgrid(track, alts)
	x, y, z = x[10:51], y[10:51], z[10:51]
	fig = plt.figure()
	ax1 = plt.subplot2grid((1,1), (0,0))
	cs = plt.pcolormesh(x, y, z, cmap='plasma', vmin=120, vmax=300)
	fig.colorbar(cs, cmap='plasma', orientation='horizontal',
				label='Kinetic Temperature / K')
	ax1.set_ylabel('Altitude / km')
	ax1.set_xlabel('Along track distance / km')
	str_plt_title = 'SABER Kinetic Temperature Profiles'
	str_plt_title += '\n' + 'Date: ' + strdate
	plt.title(str_plt_title)
	str_plt_savename = 'sabertest_' + strdate
	os.chdir('SABER')
	plt.savefig(str_plt_savename,dpi=300)
	plt.close()
	os.chdir('..')
