#!/usr/bin/env python3
"""
Aeolus data load from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Looping through dataset and plotting all orbits that travel---
----------through the 40-80 box.----------------------------------------
---v1.2---3.2.20-Updates------------------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
produces plots of key parameters from the datasets
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
sys.path.append('/media/GWR/AEOLUS/')
from phdfunctions import timeseriesplot, find_nearest
from functions import ncload

# Here I need to iterate through all. nc files and plot all of them
# into jpgs to view one after another
"""Find directory and read netCDF data"""
strdirectory = '/media/GWR/AEOLUS/NC/'

directory = os.fsencode(strdirectory)
for file in os.listdir(directory):
	
	print("\n===========================================================")
	
	# Program Timing (Time taken to get through one file)	
	fstartTime = datetime.now()
	
	# Setting filename for dataload
	filename = os.fsdecode(file)
	print(str(filename), '\n')

	infile = strdirectory + str(filename) # Specifies data file
	print('netCDF file:')
	print(infile, '\n')
	data = nc.Dataset(infile)

	"""=============================================================="""
	"""======================Download Variables======================"""
	"""=============================================================="""
	# Longitude
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

	# Converted time
	# ~ data_time = nc.num2date(data.variables['time'][:],\
	# ~ calendar = 'standard', units = data.variables['time'].units)

	"""=============================================================="""
	"""=====Test to see if orbit is sufficiently within Andes box===="""
	"""=============================================================="""
	mnopib = 30 # minimum_number_of_profiles_in_box
	# Print full arrays without truncation
	np.set_printoptions(threshold=sys.maxsize)
	# ~ print(np.where(data_lat<-80, 0, (np.where(data_lat>-40, 0, 1))))
	# Find where the satellite is within the Andes box
	box = np.where(data_lat<-80, 0, (np.where(data_lat>-40, 0,
		(np.where(data_lon>320, 0, (np.where(data_lon<280, 0, 1)))))))
	diffs = np.diff(box) # Array of diffs for box
	# Grouping the differences between the elements in diffs
	# (2nd derivative)
	grouped_diffs = [(k, sum(1 for i in g)) for k,g in groupby(diffs)]
	# Returns:[(0, 6206),(1, 1),...,(0, 1748),(-1, 1),...,(0, 8617)]
	print(grouped_diffs)
	
	# Finding the start and end elements of the desired section
	itrn = 0
	start_elmnt = 0
	end_elmnt = 0
	# Iterate through grouped_diffs, adding up continually
	for u in grouped_diffs:
		if u[0] != 0: # Bypass 1s and -1s
			itrn += u[1]
		elif u[0] == 0:
			# Is this section the Andes box?
			if data_lat[itrn] < -40 and data_lat[itrn] > -80 and \
			box[itrn] == 1:
				# Are there enough profiles in the box?
				if u[1] > mnopib:
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
	data_lon = data_lon[start_elmnt:end_elmnt+1]
	data_lat = data_lat[start_elmnt:end_elmnt+1]
	data_alt = data_alt[start_elmnt:end_elmnt+1]
	data_HLOS = data_HLOS[start_elmnt:end_elmnt+1]
	rayleigh_times = rayleigh_times[start_elmnt:end_elmnt+1]

	"""=============================================================="""
	"""===================Creating arrays for plot==================="""
	"""=============================================================="""

	# Initialise meshgrids for x, y and z
	alts = np.linspace(0,20000, 21)
	# ~ # Lists
	# ~ z = [[0 for _ in range(len(RG))] for _ in range(len(alts))] 
	# ~ z_itrn = [[0 for _ in range(len(RG))] for _ in range(len(alts))]
	z = np.zeros((len(alts),len(RG))) # NumPy Arrays
	z_itrn = np.zeros((len(alts),len(RG)))
	# ~ print(np.shape(z))
	
	# Placing wind values into bins of height 1km and width
	# 1 rayleigh group
	lastgroupstarttime = 0
	RG_start = 0
	for RG_elmnt in range(len(RG)):
		for t in range(len(rayleigh_times)):
			# Find all elements inside this sandwich and add to z
			# and z_itrn:
			if rayleigh_times[t] < RG[RG_elmnt] and \
			rayleigh_times[t] >= lastgroupstarttime:
				if RG_start == 0:
					RG_start = RG_elmnt
				# Find the nearest altitude level
				val = find_nearest(alts, data_alt[t])
				alt_elmnt = np.where(alts == val)[0][0]
				# Cap wind speeds to 250 m/s
				if np.abs(data_HLOS[t]) < 25000:
					z[alt_elmnt][RG_elmnt] += data_HLOS[t]
					z_itrn[alt_elmnt][RG_elmnt] += 1
				RG_end = RG_elmnt
		lastgroupstarttime = RG[RG_elmnt]
		
	# Find the mean for each bin
	z /= 100 * z_itrn # Factor of 100 for conversion from cm/s to m/s
	print(z)

	# Amend RG array
	RG = RG[RG_start:RG_end+1]
	z = z[:, RG_start:RG_end+1]
	
	date_time = coda.time_to_utcstring(RG[:])
	date_time = np.array([datetime.strptime(date,
		'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
	rayleigh_times = coda.time_to_utcstring(rayleigh_times[:])
	rayleigh_times = np.array([datetime.strptime(date,
		'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times])
	x, y = np.meshgrid(date_time, alts)
	# ~ print(x)
	# ~ print(y)

	"""=============================================================="""
	"""=========================Plotting============================="""
	"""=============================================================="""
	os.chdir('..')
	os.chdir('Plots')

	YYYY = str(filename)[6:10]
	MM = str(filename)[11:13]
	print(os.getcwd())
	# Enter corresponding YYYY directory
	print('\n')
	try:
		os.mkdir(YYYY)
		print("Directory ", YYYY, " created")
	except OSError as e:
			if e.errno == errno.EEXIST:
				print("Directory ", YYYY, " already exists")
			else:
				raise
	os.chdir(YYYY)

	# Enter corresponding MM directory
	try:
		os.mkdir(MM)
		print("Directory ", MM, " created")
	except OSError as e:
			if e.errno == errno.EEXIST:
				print("Directory ", MM, " already exists")
			else:
				raise
	os.chdir(MM)
	print(infile, '\n')
	# Plotting data
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	cs = plt.contourf(x,y,z, cmap='RdBu',
		levels=np.linspace(-200, 200, 41))
	ax2 = ax1.twinx()
	ax2.plot(rayleigh_times, data_lat, c='black', marker='.',
		markersize='1', label='Latitude', linewidth=0.1)

	# Setting Date axis
	date_form = '%H:%M'
	minor_date_form = '%M'
	hours = dates.HourLocator()
	minutes = dates.MinuteLocator()
	date_form = dates.DateFormatter(date_form)
	minor_date_form = dates.DateFormatter(minor_date_form)
	ax1.xaxis.set_major_formatter(date_form)
	ax1.xaxis.set_minor_formatter(minor_date_form)
	ax1.set_xlim(date_time[0], date_time[-1])
	ax1.set_xlabel('Time')

	# Setting y axes
	ax1.set_ylabel('Altitude / m')
	ax2.set_ylabel('Latitude / $^\circ$')
	plt.title('Aeolus Orbit HLOS Rayleigh Wind Cross-section')
	fig.colorbar(cs, cmap='RdBu', ax=ax1, orientation='horizontal',
		label='HLOS Rayleigh Wind Speed / ms-1')
	# ~ plt.legend(loc=9)
	pngsavename = str(filename)[:-3] + '.png'
	plt.savefig(pngsavename,dpi=300)
	
	# Climb out of plot directory
	os.chdir('..')
	os.chdir('..')
	
	# Time taken for the file
	fduration = datetime.now() - fstartTime
	print('That file took ', fduration, ' seconds')
