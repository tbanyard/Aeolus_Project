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
"""Find directory and read netCDF data"""
strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC_FullQC_Jul2020/'

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
	# Rayleigh_Grouping
	RG = data.variables['RG'][:]
	# Time
	rayleigh_times = data.variables['time'][:]

	# Converted time
	# ~ data_time = nc.num2date(data.variables['time'][:],\
	# ~ calendar = 'standard', units = data.variables['time'].units)
	
	# Both QC Flags
	data_QC = data.variables['QC_Flag_Both'][:]

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
		data_QC_new = data_QC[start_elmnt:end_elmnt+1]

		"""=========================================================="""
		"""=================Creating arrays for plot================="""
		"""=========================================================="""

		# Initialise meshgrids for x, y and z
		maxheight = 20000
		levnum = (maxheight / vert_res) + 1
		alts = np.linspace(0,maxheight, levnum)
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
							if data_QC_new[t] == 1:
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
		print(RG_start)
		print(RG_end)
		
		# Fixing time dimension to match coda format
		date_time = coda.time_to_utcstring(RG_new[:])
		date_time = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
		rayleigh_times_new = \
		coda.time_to_utcstring(rayleigh_times_new[:])
		rayleigh_times_new = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times_new])
		x, y = np.meshgrid(date_time, alts)
		
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
			fixnanswithmean(z) # Uses my own function 'fixnanswithmean'	
			
			# Test using Gaussian noise
			# ~ for i in range(len(z)):
				# ~ for j in range(len(z[0])):
					# ~ z[i][j] = random.gauss(0, 20)
			
			# ~ z = ndimage.uniform_filter(z, size=(3,7), mode = 'reflect')
			z1 = np.copy(z)
			# Calculate appropriate upper and lower bounds to band-pass
			sg_upper = int(np.ceil(15 * (1000/vert_res)))
			sg_upper += (1-(sg_upper % 2))
			sg_lower = int(np.floor(5 * (1000/vert_res)))
			sg_lower += (1-(sg_lower % 2))
			print(sg_upper * vert_res, "m upper bound")
			print(sg_lower * vert_res, "m lower bound")
			# ~ z2 = savgol_filter(z, sg_upper, 2, axis = 0) # + Vertical S-G filter
			# ~ z3 = savgol_filter(z, sg_lower, 2, axis = 0) # - Vertical S-G filter
			try:
				z2a = savgol_filter(z, 25, 2, axis = 1) # Horizontal + S-G filter
				z2b = savgol_filter(z, 3, 2, axis = 1) # Horizontal - S-G filter
				z4 = savgol_filter(z, 21, 2, axis = 0) # Vertical S-G filter
				z6a = ndimage.gaussian_filter1d(z, 0.5, axis = 1)
				z6b = ndimage.gaussian_filter1d(z, 0.5, axis = 0)
				z6 = (z6a+z6b)/2
				z7a = ndimage.gaussian_filter1d(z, 1, axis = 1)
				z7b = ndimage.gaussian_filter1d(z, 1, axis = 0)
				z7 = (z7a+z7b)/2
				z5 = (z2b+z4)/2
			except:
				continue
			# ~ z = z3 - z2
			# ~ z = z1 - z2
			# ~ z = z2b - z2a
						
			# Plots using imshow
			# ~ cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
				# ~ x_lims[1], y_lims[0], y_lims[1]], vmin=-200, vmax=200,
				# ~ interpolation='none')
			cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=-20, vmax=20,
				interpolation='none')
				
			if np.mean(z)<0:
				cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
					x_lims[1], y_lims[0], y_lims[1]], vmin=-60, vmax=-10,
					interpolation='none')
				
			elif np.mean(z)>0:
				cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
					x_lims[1], y_lims[0], y_lims[1]], vmin=60, vmax=10,
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
		print("Here:", str(filename)[6:-3]) # Here: 2019-08-11_224235
		strdate = str(filename)[6:16]
		print("strdate", strdate)
		strtime = str(filename)[17:19] + ':' + str(filename)[19:21] + ':' + \
			str(filename)[21:23]
		print(strtime)
		str_plt_title = 'Aeolus Orbit HLOS Rayleigh Wind Cross-section'
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
		plt.savefig('testme2.png',dpi=300)
		
		# Time taken for the file
		fduration = datetime.now() - fstartTime
		print('That file took ', fduration, ' seconds')
		
		# Reset after box completed
		start_elmnt = 0
		end_elmnt = 0
		complete_boxes += 1
		
		plt.close()
