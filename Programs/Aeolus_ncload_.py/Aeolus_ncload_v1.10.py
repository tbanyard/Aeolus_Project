#!/usr/bin/env python3
"""
Aeolus data load from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Looping through dataset and plotting all orbits that travel---
----------through the 40-80 box.----------------------------------------
---v1.2---3.2.20-Updates------------------------------------------------
---v1.3---Deals with multiple orbits------------------------------------
----------simple_contourf-----------------------------------------------
---v1.4---pcolormesh/imshow---------------------------------------------
---v1.5---Additional subplot with lat/lon grid showing satellite track--
---v1.6---23.03.20-Updates: Applying S-G Filter and 2D Boxcar on data---
---v1.7---Experimenting with 500m bins, with interpolation, different---
----------scheme for dealing with NaNs etc.-----------------------------
---v1.8---ERA5 interpolated onto Aeolus track---------------------------
---v1.9---Testing New NC_FullQC files, smoothing ERA5 first (in u and v)
---v1.10--Consolidating code and tidying up-----------------------------
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
from scipy.signal import savgol_filter
import scipy.ndimage as ndimage
from itertools import groupby
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Change current working directory to parent directory
os.chdir('..')

# Here I need to iterate through all. nc files and plot all of them
# into jpgs to view one after another
"""Find directory and read Aeolus netCDF data"""
strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC4_FullQC/'
"""ERA5 directory"""
ERA5_dir = '/home/tpb38/PhD/Bath/ERA5/'

# Choose pcolor or imshow
pc_or_im = 'im'
im_interp = 'none'

# Set oldinterp to True to use the old interpolation routine (nearest)
oldinterp = False

# Choose vertical resolution
vert_res = 500 # Given in metres

# Azores or Andes? (Comment in the desired region)
region = 'andes'
# ~ region = 'azores'

if region == 'andes':
	minlat, maxlat, minlon, maxlon = -80, -40, 280, 320
	bmlowlat, bmupperlat, bmleftlon, bmrightlon = -80, -40, -80, -40
elif region == 'azores':
	minlat, maxlat, minlon, maxlon = 25, 50, 300, 340
	bmlowlat, bmupperlat, bmleftlon, bmrightlon = 10, 70, -80, 0
	
# Iterate through files in directory
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
	data_HLOS = data.variables['Rayleigh_HLOS_wind_speed'][:]
	# ~ data_HLOS = data.variables['Zonal_wind_projection'][:]
	# ~ data_HLOS = data.variables['Meridional_wind_projection'][:]
	# ~ data_HLOS = data.variables['LOS_azimuth'][:]
	# Zonal Wind Projection
	data_u_wind = data.variables['Zonal_wind_projection'][:]
	# Meridional Wind Projection
	data_v_wind = data.variables['Meridional_wind_projection'][:]
	# Line-of-sight azimuth
	data_azimuth = data.variables['LOS_azimuth'][:]
		
	# Rayleigh_Grouping
	RG = data.variables['RG'][:]
	# Time
	rayleigh_times = data.variables['time'][:]

	# Converted time
	# ~ data_time = nc.num2date(data.variables['time'][:],\
	# ~ calendar = 'standard', units = data.variables['time'].units)
	data.close()

	"""=============================================================="""
	"""=====Test to see if orbit is sufficiently within Andes box===="""
	"""=============================================================="""
	mnopib = 150 # minimum_number_of_profiles_in_box
	# Print full arrays without truncation
	np.set_printoptions(threshold=sys.maxsize)
	# ~ print(np.where(data_lat<-80, 0, (np.where(data_lat>-40, 0, 1))))
	# Find where the satellite is within the Andes box
	box = np.where(data_lat<minlat, 0, (np.where(data_lat>maxlat, 0,
		(np.where(data_lon>maxlon, 0, (np.where(data_lon<minlon, 0, 1)))))))
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
	# Iterate through grouped_diffs, adding up continually
	for u in grouped_diffs:
		if u[0] != 0: # Bypass 1's and -1's
			itrn += u[1]
		elif u[0] == 0:
			# Is this section the Andes box?
			if data_lat[itrn] < maxlat and data_lat[itrn] > minlat and \
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
		data_lon_new = data_lon[start_elmnt:end_elmnt+1]
		data_lat_new = data_lat[start_elmnt:end_elmnt+1]
		data_alt_new = data_alt[start_elmnt:end_elmnt+1]
		data_HLOS_new = data_HLOS[start_elmnt:end_elmnt+1]
		rayleigh_times_new = rayleigh_times[start_elmnt:end_elmnt+1]
		data_u_wind = data_u_wind[start_elmnt:end_elmnt+1]
		data_v_wind = data_v_wind[start_elmnt:end_elmnt+1]
		data_azimuth = data_azimuth[start_elmnt:end_elmnt+1]

		"""=========================================================="""
		"""=====================ERA5 Interpolation==================="""
		"""=========================================================="""
		# =========================Load ERA5 Data=======================
		# Read data for the current date
		YYYY = str(filename)[6:10]
		MM = str(filename)[11:13]
		DD = str(filename)[14:16]
		doy = yyyymmdd_to_doy(YYYY, MM, DD)
		ncfile = ERA5_dir + str(YYYY) + '/era5_' + str(YYYY) + 'd' \
		+ doy + '.nc'
		print('ERA5 netCDF file:')
		print(ncfile, '\n')
		data = ERA5_dataload(ncfile)
		
		# Create special variables for ERA5 interpolation
		data_lon_new_ae = np.copy(data_lon_new)
		data_lat_new_ae = np.copy(data_lat_new)
		data_alt_new_ae = np.copy(data_alt_new)
		
		# Pass ERA5 data into their respective variables
		ERA5_data_lon, ERA5_data_lat, ERA5_data_lev, \
		ERA5_data_temp, ERA5_data_u, ERA5_data_v, \
		ERA5_data_time = data[0], data[1], data[2], data[3], \
		data[4], data[5], data[6]
		
		# Read data for the next date
		nextday_dt = yyyymmdd_nextday(YYYY, MM, DD)
		YYYY2, MM2, DD2 = nextday_dt.year, nextday_dt.month, \
		nextday_dt.day
		doy2 = yyyymmdd_to_doy(YYYY2, MM2, DD2)
		# Specify ERA5 data file
		ncfile = ERA5_dir + str(YYYY2) + '/era5_' + str(YYYY2) + \
		'd' + doy2 + '.nc'
		print('Following ERA5 netCDF file:')
		print(ncfile, '\n')
		data = ERA5_dataload(ncfile)
		
		# Pass ERA5 data for the next date into their respective
		# variables
		ERA5_data_temp, ERA5_data_u, ERA5_data_v, ERA5_data_time = \
		np.append(ERA5_data_temp, data[3], axis=0), \
		np.append(ERA5_data_u, data[4], axis=0), \
		np.append(ERA5_data_v, data[5], axis=0), \
		np.append(ERA5_data_time, data[6], axis=0)
		
		# Load any text files into arrays
		ERA5_p_levels = np.loadtxt("ERA5_pressure_levels.txt")*100
		ERA5_altitudes = np.loadtxt("ERA5_altitudes.txt")
		ERA5_levels = np.loadtxt("ERA5_pressure_levels_labels.txt")
		# ==============================================================
		
		# =====================Smooth u and v fields====================
		print(np.shape(ERA5_data_u))
		
		# Here, altdiff, latdiff and londiff give the fwhm window size in km.
		# From this, the fwhm window size in points is found and converted to sd
		
		# Calculate latitude difference
		lats = np.deg2rad(np.array([-50,-51.5]))
		lons = np.deg2rad(np.array([-70,-70]))
		latdiff = haversine(lats, lons)
		lat_sig = (3 / latdiff) / 2.355 # This is for 3km
		print("lat_sig: ", lat_sig)
		
		# Calculate longitude difference
		lats = np.deg2rad(np.array([-50,-50]))
		lons = np.deg2rad(np.array([-70,-71.5]))
		londiff = haversine(lats, lons)
		lon_sig = (3 / londiff) / 2.355 # This is for 3km
		print("lon_sig: ", lon_sig)
		
		ERA5_data_u = ndimage.gaussian_filter1d(ERA5_data_u, lat_sig, axis = 2)
		ERA5_data_u = ndimage.gaussian_filter1d(ERA5_data_u, lon_sig, axis = 3)
		ERA5_data_v = ndimage.gaussian_filter1d(ERA5_data_v, lat_sig, axis = 2)
		ERA5_data_v = ndimage.gaussian_filter1d(ERA5_data_v, lon_sig, axis = 3)
		
		new_ERA5_data_u = np.copy(ERA5_data_u)
		new_ERA5_data_v = np.copy(ERA5_data_v)
		
		# Iterate through model levels and find alt difference with next level
		for lev in range(111):
			if lev != 136:
				altdiff = ERA5_altitudes[lev] - ERA5_altitudes[lev+1]
			elif lev == 136:
				altdiff = ERA5_altitudes[lev]
			alt_sig = (10000 / altdiff) / 2.355
			print(lev + 1)
			print(alt_sig)
			
			if lev != 110:
				temp_ERA5_data_u = np.copy(ndimage.gaussian_filter1d(ERA5_data_u, alt_sig, axis = 1))
				new_ERA5_data_u[:, lev, :, :] = temp_ERA5_data_u[:, lev, :, :]
				temp_ERA5_data_v = np.copy(ndimage.gaussian_filter1d(ERA5_data_v, alt_sig, axis = 1))
				new_ERA5_data_v[:, lev, :, :] = temp_ERA5_data_v[:, lev, :, :]
			if lev == 110:
				temp_ERA5_data_u = np.copy(ndimage.gaussian_filter1d(ERA5_data_u, alt_sig, axis = 1))
				new_ERA5_data_u[:, 110:-1, :, :] = temp_ERA5_data_u[:, 110:-1, :, :]
				temp_ERA5_data_v = np.copy(ndimage.gaussian_filter1d(ERA5_data_v, alt_sig, axis = 1))
				new_ERA5_data_v[:, 110:-1, :, :] = temp_ERA5_data_v[:, 110:-1, :, :]
		
		ERA5_data_u = new_ERA5_data_u
		ERA5_data_v = new_ERA5_data_v
				
		# ==============================================================
				
		# =======================Run Interpolation======================
		# Create a datetime format version of rayleigh_times_new
		ERA5_u_interpolated = np.empty(0)
		ERA5_v_interpolated = np.empty(0)
		ERA5_interpolated = np.empty(0)
		no_result = 0
		rayleigh_times_new_dt = \
		coda.time_to_utcstring(rayleigh_times_new[:])
		rayleigh_times_new_dt = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times_new_dt])
		for t in range(len(rayleigh_times_new_dt)):
			# Finding the correct time array elements  
			rounded_down = round_down_3_hours(rayleigh_times_new_dt[t])
			rounded_up = round_up_3_hours(rayleigh_times_new_dt[t])
			down_elmt = np.where(ERA5_data_time == rounded_down)
			up_elmt = np.where(ERA5_data_time == rounded_up)
			
			for each in ['u', 'v', 'hlos']:
				# If 'hlos' then we have already got a u and v value so:
				if each == 'hlos':
					hlos = (-ERA5_u_interpolated[t] * np.sin(data_azimuth[t] * (np.pi/180))) - \
						(ERA5_v_interpolated[t] * np.cos(data_azimuth[t] * (np.pi/180)))
					ERA5_interpolated = np.append(ERA5_interpolated, hlos)
					# Now we skip the rest of the loop
					continue
					
				# Toggle between ERA5 sandwich times at 3hr intervals
				# and interpolate between them
				for t_E in [down_elmt, up_elmt]:
					t_E = t_E[0][0] # t_E = ERA5 time element
					if t_E == down_elmt[0][0]:
						time_toggle = 0
					elif t_E == up_elmt[0][0]:
						time_toggle = 1
					else:
						print("Error")
						
					# Call run_interpolation to execute 3D spatial
					# interpolation of ERA5 data onto Aeolus track
					if each == 'u':
						ERA5_var = ERA5_data_u
					elif each == 'v':
						ERA5_var = ERA5_data_v
					result = run_interpolation(data_lat_new_ae,
					data_lon_new_ae, data_alt_new_ae, 
					ERA5_data_lat, ERA5_data_lon, ERA5_altitudes, ERA5_var, t, t_E)
					# Raise flag if there is no result
					if result == None:
						no_result = 1
						continue
					# Set up and down times
					if time_toggle == 0:
						down = result
						# ~ print("Down = ", down[0])			
					elif time_toggle == 1:
						up = result
						# ~ print("Up = ", up[0])
					if no_result == 1:
						break
				
				# Break out of loop if there is no result		
				if no_result == 1:
					if each == 'u':
						ERA5_u_interpolated = np.append(ERA5_u_interpolated, 99999)
					elif each == 'v':
						ERA5_v_interpolated = np.append(ERA5_v_interpolated, 99999)
					no_result = 0
					continue
				
				# Linear interpolation in time using linearinterp2kp
				# function
				x_0 = toTimestamp(ERA5_data_time[down_elmt[0][0]])
				x_1 = toTimestamp(ERA5_data_time[up_elmt[0][0]])
				x = toTimestamp(rayleigh_times_new_dt[t])
				y_0 = down[0]
				y_1 = up[0]
				y = linearinterp2kp(x, x_0, x_1, y_0, y_1)
				
				if each == 'u':
					ERA5_u_interpolated = np.append(ERA5_u_interpolated, y)
				elif each == 'v':
					ERA5_v_interpolated = np.append(ERA5_v_interpolated, y)
			
		# ~ print(len(ERA5_interpolated))
		# ~ print(len(rayleigh_times_new))
		# ~ data_HLOS_new = np.copy(ERA5_interpolated)*100 
		# ~ print(data_HLOS_new)
		"""=========================================================="""
		"""=================Creating arrays for plot================="""
		"""=========================================================="""

		# Initialise meshgrids for x, y and z
		maxheight = 20000
		levnum = (maxheight / vert_res) + 1
		alts = np.linspace(0,maxheight, levnum)
		
		# Detrending to the smoothed model
		data_obs = np.copy(data_HLOS_new) # Store observations for later use
		data_mod = np.copy(ERA5_interpolated)*100 # Converts ERA5 data into cm/s
		for each in ['mod', 'obs']:
			if each == 'mod':
				data_HLOS_new = data_mod
			elif each == 'obs':
				data_HLOS_new = data_obs
		
			# ~ #Lists
			# ~ z = [[0 for _ in range(len(RG))] for _ in range(len(alts))]
			# ~ z_itrn = \
			# ~ [[0 for _ in range(len(RG))] for _ in range(len(alts))]
			z = np.zeros((len(alts),len(RG))) # NumPy Arrays
			z_itrn = np.zeros((len(alts),len(RG)))
			# ~ print(np.shape(z))
			
			# ============OLD INTERPOLATION ROUTINE==============
			if oldinterp == True:
				# Placing wind values into bins of height 1km and width
				# 1 rayleigh group
				lastgroupstarttime = 0
				RG_start = 0
				for RG_elmnt in range(len(RG)):
					for t in range(len(rayleigh_times_new)):
						# Find all elements inside this sandwich and add to z
						# and z_itrn:
						if rayleigh_times_new[t] < RG[RG_elmnt] and \
						rayleigh_times_new[t] >= lastgroupstarttime:
							if RG_start == 0:
								RG_start = RG_elmnt
							# Find the nearest altitude level
							val = find_nearest(alts, data_alt_new[t])
							alt_elmnt = np.where(alts == val)[0][0]
							# Cap wind speeds to 250 m/s
							if np.abs(data_HLOS_new[t]) < 25000:
								z[alt_elmnt][RG_elmnt] += data_HLOS_new[t]
								z_itrn[alt_elmnt][RG_elmnt] += 1
							RG_end = RG_elmnt
					lastgroupstarttime = RG[RG_elmnt]
			# ===================================================			
			else:
			# ============NEW INTERPOLATION ROUTINE==============
				z_new = np.zeros((len(alts),len(RG))) # NumPy Arrays
				points = np.empty(0)
				values = np.empty(0)
				lastgroupstarttime = 0
				RG_start = 0
				for RG_elmnt in range(len(RG)):
					for t in range(len(rayleigh_times_new)):
						# Find all elements inside this sandwich and add to z
						# and z_itrn:
						if rayleigh_times_new[t] < RG[RG_elmnt] and \
						rayleigh_times_new[t] >= lastgroupstarttime:
							if RG_start == 0:
								RG_start = RG_elmnt
							# Cap wind speeds to 250 m/s
							if np.abs(data_HLOS_new[t]) < 25000:
								if data_alt_new[t] in points:
									itis = np.where(points == data_alt_new[t])
								points = np.append(points, data_alt_new[t])
								values = np.append(values, data_HLOS_new[t])
							RG_end = RG_elmnt
					lastgroupstarttime = RG[RG_elmnt]
					z_new[:, RG_elmnt] = griddatainterpolation(points, values, alts)
					points = np.empty(0)
					values = np.empty(0)
			# ===================================================
								
			# Find the mean for each bin (Toggled on for old interpolation)
			if oldinterp == True:
				z /= 100 * z_itrn # Factor of 100 for conversion from cm/s - m/s
			else:
				# Using new interpolation routine
				z = np.copy(z_new) / 100
					
			# Amend RG array
			RG_new = RG[RG_start:RG_end+1]
			z = z[:, RG_start:RG_end+1]
			
			if each == 'mod':			
				data_mod = np.copy(z)
			elif each == 'obs':
				data_obs = np.copy(z)
	
		# Applying NaN mask
		grayhatchescmap = LinearSegmentedColormap('Grayhatchescmap',
			segmentdata=customcolormaps('grayhatches'), N=265)
		isnanarray = np.isnan(z)
		hatcharray = np.copy(z)
		for horidx in range(len(isnanarray)):
			for veridx in range(len(isnanarray[horidx])):
				if isnanarray[horidx][veridx] == True:
					hatcharray[horidx][veridx] = 1
				else:
					hatcharray[horidx][veridx] = 0
		
		# Fixing time dimension to match coda format
		date_time = coda.time_to_utcstring(RG_new[:])
		date_time = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
		rayleigh_times_new = \
		coda.time_to_utcstring(rayleigh_times_new[:])
		rayleigh_times_new = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times_new])
		x, y = np.meshgrid(date_time, alts)
					
		"""=========================================================="""
		"""=======================Plotting==========================="""
		"""=========================================================="""
		os.chdir('..')
		print(os.getcwd())
		os.chdir('Plots')
		# Save data (toggle on and off)
		# ~ savemat("data.mat", {'data':z})
		# ~ np.save("data.npy", [x, y, z], allow_pickle=True)
		# ~ np.save("rayleigh_times.npy", rayleigh_times_new, allow_pickle=True)
		# ~ np.save("date_time.npy", date_time, allow_pickle=True)
		# ~ np.save("data_track.npy", [data_lat_new, data_lon_new], allow_pickle=True)
		
		YYYY = str(filename)[6:10]
		MM = str(filename)[11:13]

		print(infile, '\n')
		# Plotting data
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
			y_lim_max = (maxheight + (vert_res/2)) / 1000
			y_lim_min = (0 - (vert_res/2)) / 1000
			y_lims = [y_lim_max, y_lim_min]
			
			# ======================= Pre ERA5 Scheme ======================== #
			fixnanswithmean(z) # Uses my own function 'fixnanswithmean'	
			
			# Test using Gaussian noise
			# ~ for i in range(len(z)):
				# ~ for j in range(len(z[0])):
					# ~ z[i][j] = random.gauss(0, 20)
			
			# ~ z = ndimage.uniform_filter(z, size=(3,7), mode = 'reflect') #N/A
			z1 = np.copy(z)
			# Calculate appropriate upper and lower bounds to band-pass
			sg_upper = int(np.ceil(15 * (1000/vert_res)))
			sg_upper += (1-(sg_upper % 2))
			sg_lower = int(np.floor(5 * (1000/vert_res)))
			sg_lower += (1-(sg_lower % 2))
			print(sg_upper * vert_res, "m upper bound")
			print(sg_lower * vert_res, "m lower bound")
			z2 = savgol_filter(z, sg_upper, 2, axis = 0) # + Vertical S-G filter
			z3 = savgol_filter(z, sg_lower, 2, axis = 0) # - Vertical S-G filter
			# ~ z2 = savgol_filter(z, 15, 2, axis = 1) # Horizontal S-G filter
			z = z3 - z2 # Return the perturbations
			# ~ z = z1 - z2
			# ======================= Pre ERA5 Scheme ======================== #
			# ~ fixnanswithmean(data_obs)
			# ~ z = data_obs - data_mod # Find perturbations relative to smoothed ERA5
			z = data_mod
			# ~ z = data_obs
						
			# Plots using imshow
			cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=-150, vmax=150,
				interpolation='none')
			mask = plt.imshow(hatcharray, aspect='auto', cmap=grayhatchescmap,
				extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], interpolation=im_interp)
			ax1.xaxis_date() # Initialises date axis
			date_form = dates.DateFormatter('%H:%M') # Sets date format
			ax1.xaxis.set_major_formatter(date_form)
			plt.gca().invert_yaxis() # Invert axis for imshow

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
		ax1.set_yticks(np.arange(21))
		ax1.yaxis.set_major_locator(plt.MaxNLocator(11))
		ax1.yaxis.set_minor_locator(plt.MaxNLocator(21))
		# ~ ax1.tick_params(axis='y', which='minor', left=True) # Minor ticks
		ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='y',
			which='both')

		# Satellite track plot
		ax2 = plt.subplot2grid((5,5), (0,3), colspan=2, rowspan=4)
		ax2.yaxis.set_label_position("right")
		ax2.yaxis.tick_right()
		map = Basemap(projection='cyl',llcrnrlat=bmlowlat,urcrnrlat=bmupperlat,\
					llcrnrlon=bmleftlon,urcrnrlon=bmrightlon,resolution='i', ax=ax2)
		# ~ map = Basemap(projection='cyl',llcrnrlat=-80,urcrnrlat=-40,\
					# ~ llcrnrlon=-80,urcrnrlon=-40,resolution='i', ax=ax2)
		map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
		map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
		map.drawcoastlines(linewidth=0.25, color='#666666')
		if region == 'andes':
			map.drawmeridians([-70, -60, -50], linewidth=0.3)
			map.drawparallels([-70, -60, -50], linewidth=0.3)
		elif region == 'azores':
			map.drawmeridians([-60, -40, -20], linewidth=0.3)
			map.drawparallels([20, 30, 40, 50, 60], linewidth=0.3)
		map.scatter(data_lon_new - 360, data_lat_new, marker = 'x', color = 'red',
			s=0.3, zorder=2)
		map.plot(data_lon_new - 360, data_lat_new, color = 'red', linewidth = '0.5')
		midpointlonindex = int(np.floor(len(data_lon_new)/2))
		midpointlatindex = int(np.floor(len(data_lat_new)/2))
		# ~ map.quiver(data_lon[midpointlonindex]-360, data_lat[midpointlatindex],
			# ~ data_lon[midpointlonindex+100] - data_lon[midpointlonindex-100],
			# ~ data_lat[midpointlatindex+100] - data_lat[midpointlatindex-100],
			# ~ color='black', zorder=3, headaxislength=0, headlength=0, pivot='middle')
		# ~ map.quiver(data_lon[-1]-360, data_lat[-1], data_lon[-2]-data_lon[-1],
			# ~ data_lat[-2]-data_lat[-1], color='black', zorder=4)

		# Add lines at time intervals corresponding to date_time ticks
		date_ticks = dates.num2date(ax1.xaxis.get_ticklocs())

		for i in range(len(date_ticks)):
			# ~ print(date_ticks[i].minute)
			# Co-locate between ticks and rayleigh_times array
			nearest_date = find_nearest(rayleigh_times_new, date_ticks[i])
			idx = (np.where(dates.date2num(rayleigh_times_new) == nearest_date)[0])[0]
			str_time = date_ticks[i].strftime('%H:%M') # Sets string form of each tick
			# Location of each line
			lon_loc = data_lon_new[idx]-360
			lat_loc = data_lat_new[idx]
			# ~ map.scatter([lon_loc], [lat_loc], s=50, marker='+', zorder=111,
				# ~ color='black')
			# Plots line on map and annotates with text
			try: 
				data_lat_new[idx+30] == False
				data_lat_new[idx-30] == False
				data_lon_new[idx+30] == False
				data_lon_new[idx-30] == False
			except:
				continue
			map.quiver(lon_loc, lat_loc, data_lat_new[idx-30] - data_lat_new[idx+30],
				data_lon_new[idx+30] - data_lon_new[idx-30], angles='xy', color='black',
				zorder=3, headaxislength=0, headlength=0, pivot='middle', units='xy')
			ax2.annotate(str_time, (lon_loc+1, lat_loc+1), fontsize=6, zorder=5)
			
			# ~ a = ax1.xaxis.get_ticklabels()
			# ~ print(text.Text(agg_filter=a))

		# ~ map.plot([-70, -65], [-70, -67], color='k')

		# Fix axes
		if region == 'andes':
			ax2.set_xticks([-80, -70, -60, -50, -40])
			ax2.set_yticks([-80, -70, -60, -50, -40])
		elif region == 'azores':
			ax2.set_xticks([-80, -60, -40, -20, 0])
			ax2.set_yticks([10, 20, 30, 40, 50, 60, 70])
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
		print("This file is:", str(filename)[6:-3]) # Here: 2019-08-11_224235
		strdate = str(filename)[6:16]
		print("strdate:", strdate)
		strtime = str(filename)[17:19] + ':' + str(filename)[19:21] + ':' + \
			str(filename)[21:23]
		print("strtime:", strtime)
		str_plt_title = 'Aeolus Orbit HLOS Rayleigh Wind Cross-section'
		# ~ str_plt_title = 'ERA5 interpolated onto Aeolus Orbit S-G 5-15km Band Pass'
		str_plt_title += '\n' + 'Orbit: ' + strdate + ' ' + strtime
		plt.title(str_plt_title, y=15)

		# ~ plt.legend(loc=9)
		pngsavename = str(filename)[:-3]
		if complete_boxes != 0:
			pngsavename += '_orb' + str(complete_boxes+1)
		if pc_or_im == 'pc':
			pngsavename += '_pc'
		elif pc_or_im == 'im':
			pngsavename += '_im_' + im_interp
		pngsavename += '.png'
		print(pngsavename)
		enterdirectory(pc_or_im)
		if pc_or_im == 'im':
			enterdirectory(im_interp)
			
			# Enter corresponding YYYY directory
			enterdirectory(YYYY)
		
			# Enter corresponding MM directory
			enterdirectory(MM)
			
			plt.savefig(pngsavename,dpi=300)
			os.chdir('..')
		else:
			# Enter corresponding YYYY directory
			enterdirectory(YYYY)
		
			# Enter corresponding MM directory
			enterdirectory(MM)
			
			# ~ plt.savefig(pngsavename,dpi=300)
		os.chdir('..')
		
		# Climb out of plot directory
		os.chdir('..')
		os.chdir('..')
		plt.savefig('testme7_10kmmod2.png',dpi=300)
		
		# Time taken for the file
		fduration = datetime.now() - fstartTime
		print('That file took ', fduration, ' seconds')
		
		# Reset after box completed
		start_elmnt = 0
		end_elmnt = 0
		complete_boxes += 1
		
		plt.close()
