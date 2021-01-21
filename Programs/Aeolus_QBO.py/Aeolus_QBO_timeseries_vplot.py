#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Ascending/Descending Node Split-------------------------------
---v1.2---Daily_Means and creating .mat file for Neil-------------------
---v1.3---Plotting my own version of this plot extended out to March----
---v1.4---Using generated netCDF file from v1.3-------------------------
---v1.5---Editing and improving plot------------------------------------
---vplot--Plotting data from the .nc file created by vcalc--------------
----------[BRANCH]-This_is_a_working_branch_version_of_this_file--------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
generates zonal mean U plot at equator
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
from scipy.signal import savgol_filter
import time
from scipy.interpolate import griddata, interp2d
from scipy.io import savemat
import scipy.ndimage as ndimage
from itertools import groupby
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.contour import ClabelText
from matplotlib import colorbar

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Program Timing	
pstartTime = datetime.now()

# Change current working directory to parent directory
os.chdir('..')

# Find directory and read netCDF data
strdirectory = '/home/tpb38/PhD/Bath/Aeolus_Project/Programs/'
infile = strdirectory + 'qbo-jan13th2021_10deg.nc' # Specifies data file
print('netCDF file:')
print(infile, '\n')
data = nc.Dataset(infile)

qbocmap = LinearSegmentedColormap('QBOcustomcmap', segmentdata=customcolormaps('QBOcmap3'), N=265)
grayhatchescmap = LinearSegmentedColormap('Grayhatchescmap', segmentdata=customcolormaps('grayhatches'), N=265)
grayhatchescmap_r = grayhatchescmap.reversed()
whitehatchescmap_r = LinearSegmentedColormap('Whitehatchescmap', segmentdata=customcolormaps('whitehatches'), N=265)
whitehatchescmap = whitehatchescmap_r.reversed()

"""=============================================================="""
"""======================Download Variables======================"""
"""=============================================================="""
# Altitude
data_alt = data.variables['altitude'][:]
# Zonal Wind projection of the HLOS wind
data_u_proj = data.variables['Zonal_wind_projection'][:]
# Converted time
date_time = nc.num2date(data.variables['time'][:],\
calendar = 'standard', units = data.variables['time'].units)

alts = data_alt
# ~ print("Alts: ", alts)
# ~ print("Shape of Alts: ", np.shape(alts))
time = date_time[140:] # Change to [140::7] to see what one day per week is like
# ~ print("Time: ", time)
# ~ print("Shape of Time: ", np.shape(time))
z = data_u_proj[:,140:]
# ~ print("z: ", z)
# ~ print("Shape of z: ", np.shape(z))

# ~ date_time = np.array([dates.date2num(date) for date in date_time])

print(datetime.strftime(date_time[140], '%Y-%m-%d %H:%M:%S.%f'))

# These two commands take the real_datetime, convert to string, and then back to datetime.datetime
date_time = np.array([datetime.strftime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
date_time = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
# The above code is not needed, but is a nice trick to convert from
# real_datetime to datetime.datetime
			
date_time = np.array([dates.date2num(date) for date in date_time])
date_time = date_time[140:]

np.set_printoptions(threshold=sys.maxsize)

x, y = np.meshgrid(date_time, alts)

fig = plt.figure()
ax1 = fig.add_subplot(111)

#===========

imshow = False
interp = True

# Flatten x_meshgrid and convert from datetime format to mpl_compatible nums
x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
# ~ x_lims = dates.date2num(x_lims)
# Y limits
y_lim_max = 30.5
y_lim_min = -0.5
y_lims = [y_lim_max, y_lim_min]

# Set NaNs to mean and create binary array of NaNs.
isnanarray = np.isnan(z)
hatcharray = np.copy(z)
mean = np.nanmean(z)
for xidx in range(len(isnanarray)):
	for yidx in range(len(isnanarray[xidx])):
		if isnanarray[xidx][yidx] == True:
			# ~ z[xidx][yidx] = mean
			hatcharray[xidx][yidx] = 1
		else:
			hatcharray[xidx][yidx] = 0
hatcharray = ndimage.uniform_filter(hatcharray, size=(1,5), mode = 'reflect')
hatcharray = savgol_filter(hatcharray, 27, 2, axis = 1)
hatcharray = savgol_filter(hatcharray, 9, 2, axis = 1)
# ~ hatcharray = ndimage.gaussian_filter(hatcharray, sigma=0.5, order=0)
# ~ hatcharray = ndimage.gaussian_filter(hatcharray, sigma=0.5, order=0) 

fixnanswithmean(z)
z = ndimage.uniform_filter(z, size=(1,10), mode = 'reflect')
z = savgol_filter(z, 5, 2, axis = 1) # S-G filter

if interp == True:
	# 2D interpolation to round edges
	f = interp2d(x[0], y[:,0], z, kind='cubic') # Have tried linear and quintic
	xi2 = np.linspace(x[0][0], x[0][-1], 500)
	yi2 = np.linspace(y[:,0][0], y[:,0][-1], 100)
	zi2 = f(xi2, yi2)

# ~ z1 = np.copy(z)
# Calculate appropriate upper and lower bounds to band-pass
# ~ sg_upper = int(np.ceil(15 * (1000/vert_res)))
# ~ sg_upper += (1-(sg_upper % 2))
# ~ sg_lower = int(np.floor(5 * (1000/vert_res)))
# ~ sg_lower += (1-(sg_lower % 2))
# ~ print(sg_upper * vert_res, "m upper bound")
# ~ print(sg_lower * vert_res, "m lower bound")
# ~ z2 = savgol_filter(z, sg_upper, 2, axis = 0) # + Vertical S-G filter
# ~ z3 = savgol_filter(z, sg_lower, 2, axis = 0) # - Vertical S-G filter
# ~ z2 = savgol_filter(z, 15, 2, axis = 1) # Horizontal S-G filter
# ~ z = z3 - z2
# ~ z = z1 - z2

if imshow == True:
	# Plot using imshow
	cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
		x_lims[1], y_lims[0], y_lims[1]], vmin=-20, vmax=20,
		interpolation='none')
	plt.gca().invert_yaxis() # Invert axis for imshow

elif imshow == False:
	if interp == False:
		cs = plt.contourf(x,y/1000,z, cmap=qbocmap,
			levels=np.linspace(-25, 25, 51), vmin=-22, vmax=22)
	elif interp == True:
		cs = plt.contourf(xi2,yi2/1000,zi2, cmap=qbocmap,
			levels=np.linspace(-25, 25, 51), vmin=-22, vmax=22)
			
	try:
		cs4 = ax1.contour(xi2, yi2/1000, zi2, levels=[-20,-10,10,20], linewidths = 0.4, colors='k', vmin=-30, vmax=30, alpha=1)
	except:
		print("Unable to plot contours")
	else:
		cls = ax1.clabel(cs4, [-20,-10,10,20], fmt = '%i', fontsize=4, inline_spacing=1)
	try:		
		cs5 = ax1.contour(xi2, yi2/1000, zi2, levels=[0], linewidths = 0.8, colors = 'k', alpha=1)
	except:
		print("Unable to plot zero contour")
	else:
		# ~ ClabelText(x=xi2[140], y=5, text='0', color='r', zorder = 5)
		cls2 = ax1.clabel(cs5, [0], fmt = '%i', fontsize=4, inline_spacing=1)
		

mask = plt.contourf(x, y/1000, hatcharray, cmap = whitehatchescmap, zorder = 10, levels=1)
mask = plt.contour(x, y/1000, hatcharray, linestyles = 'dashed', linewidths = 1.0, cmap = grayhatchescmap, zorder = 10, levels=1)
# ~ mask = ax1.pcolor(x, y/1000, hatcharray, cmap = grayhatchescmap, facecolor = 'r', edgecolor = 'none')

# ~ mask = plt.imshow(hatcharray, aspect='auto', cmap=grayhatchescmap,
	# ~ extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], interpolation=im_interp)
ax1.xaxis_date() # Initialises date axis
date_form = dates.DateFormatter('%b') # Sets date format
ax1.xaxis.set_major_formatter(date_form)

print("here")
### X axis
# ~ ax1.set_xlabel('Time')
# Ensure the number of date ticks is sensible
# ~ timerange = date_time[-1] - date_time[0]
# ~ date_intvl = int(np.ceil(timerange.seconds/(5*60)))
# ~ ax1.xaxis.set_major_locator(plt.MaxNLocator(5)) # Maximum number of date ticks
# ~ ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=date_intvl))
# ~ ax1.grid(color='k', linestyle = 'dashed', linewidth = 0.25, axis='x')
print("here")
### Y axis
ax1.set_ylabel('Altitude / km')
ax1.set_yticks(np.arange(31))
ax1.yaxis.set_major_locator(plt.MaxNLocator(16))
ax1.yaxis.set_minor_locator(plt.MaxNLocator(31))
ax2 = ax1.twinx()
ax2.set_yticks(np.arange(31))
ax2.yaxis.set_major_locator(plt.MaxNLocator(16))
ax2.yaxis.set_minor_locator(plt.MaxNLocator(31))

# ~ ax1.tick_params(axis='y', which='minor', left=True) # Minor ticks
# ~ ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='y',
	# ~ which='both')
	
plt.title('Aeolus Zonal Mean U-component of HLOS Rayleigh Wind\n$\pm$5$^{{\circ}}$ Latitude (5-day mean) 2019-2020')

# Add colorbar to figure
fig.subplots_adjust(bottom=0.225, right=0.88, left=0.12)
cbar_ax = fig.add_axes([0.12, 0.125, 0.76, 0.03])
# ~ fig.colorbar(cs, cmap=qbocmap, orientation='horizontal',
	# ~ label='U-component of HLOS Rayleigh Wind Speed / ms$^{-1}$', cax=cbar_ax,
	# ~ boundaries = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], extend='both')

colorbar.ColorbarBase(cbar_ax, cmap = qbocmap, orientation='horizontal',
		label='U-component of HLOS Rayleigh Wind Speed / ms$^{-1}$', boundaries = np.linspace(-22,22,23), ticks=np.linspace(-30,30,7), extend='both')

ax1.grid(which='both', axis='y', color='k', linewidth=0.1, linestyle='dashed', zorder=2)

pngsavename = 'filejan13th.png'
plt.savefig(pngsavename,dpi=300)
print(os.getcwd())
print("here")

sys.exit(0) # Do not continue onto 2D Test Figure? (Toggle on/off)


#===========


# ~ plt.legend(loc=9)
pngsavename = 'file804.png'
plt.savefig(pngsavename,dpi=300)

# Climb out of plot directory
os.chdir('..')
os.chdir('..')



# Initialise meshgrids for x, y and z
alts = np.linspace(0,30000, 31)
z = np.zeros((len(alts), day_itrn+1))
# ~ date_time = np.empty(1, dtype='datetime64')
date_time = np.empty(1) # Use to create netcdf_file
daynum = 1 # Initial day number
monthnum = 1 # Initial month number
day_itrn = 0 # Day iterator
currnext = 0 # Start on the current day
daydt = 0 # Day in datetime format (Set to 0 initially for brevity)

# Initialise orbit arrays
y = np.zeros(len(alts))
y_itrn = np.zeros(len(alts))
y2 = np.zeros(len(alts))
y2_itrn = np.zeros(len(alts))
orbittrigger = 0

# Generate time_units string for .nc file using datetime.datetime array
time_units = "hours since " + "1900-01-01" + " 00:00:00"

for file in sorted(os.listdir(directory)):
	
	print("\n=========================================================")
	
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
	# Zonal Wind projection of the HLOS wind
	data_u_proj = data.variables['Zonal_wind_projection'][:]
	# Time
	rayleigh_times = data.variables['time'][:]
	
	"""=============================================================="""
	"""======Test to see where orbit is within Equatorial band======="""
	"""=============================================================="""
	# Print full arrays without truncation
	np.set_printoptions(threshold=sys.maxsize)
	# Find where the satellite is within the equatorial band
	lat_band = np.where(data_lat<-5, 0, (np.where(data_lat>5, 0, 1)))
	
	for t in range(len(rayleigh_times)):
		now = coda.time_to_utcstring(rayleigh_times[t])
		nowstrp = datetime.strptime(now, '%Y-%m-%d %H:%M:%S.%f')
		now_day = nowstrp.day
		now_month = nowstrp.month
			
		if daydt == 0:
			daydt = nowstrp
		
		# Compare current day (nowday) with the daynumber (daynum)
		if orbittrigger == 0:
			# Move forward in time:
			if now_day != daynum or now_month != monthnum:
				orbittrigger = 1 # Only do this once per orbit
				currnext = 1 # Use arrays corresponding to the next date
				date_stamp = daydt
				print(date_stamp)
				date_stamp = date_stamp.replace(hour=0)
				date_stamp = date_stamp.replace(minute=0)
				date_stamp = date_stamp.replace(second=0)
				date_stamp = date_stamp.replace(microsecond=0)
				daynum = now_day
				monthnum = now_month
				daydt = nowstrp
			
		# Is element within equatorial band?
		if lat_band[t] == 1:
			
			# Find the nearest altitude level
			val = find_nearest(alts, data_alt[t])
			alt_elmnt = np.where(alts == val)[0][0]
			
			# Cap wind speeds to 250 m/s
			if np.abs(data_u_proj[t]) < 250000:
				
				if currnext == 0:
					y[alt_elmnt] += data_u_proj[t]
					y_itrn[alt_elmnt] += 1
				
				elif currnext == 1:
					y2[alt_elmnt] += data_u_proj[t]
					y2_itrn[alt_elmnt] += 1
				
				"""# Ascending or Descending node?	
				try:
					data_lat[t+50]
				except:
					delta_lat = data_lat[t] - data_lat[t-50]
				else:
					delta_lat = data_lat[t+50] - data_lat[t]
									
				if delta_lat > 0:
					node = 1 # Ascending node
				elif delta_lat < 0:
					node = -1 # Descending node
				
				# Add data to orbit array
				if node == 1 or node == -1:
					y[alt_elmnt] += data_HLOS[t] * node
					y_itrn[alt_elmnt] += 1"""
		
	# Proceed with daily mean calculation?
	if currnext == 0:
		continue # Don't execute below code until the next day
	currnext = 0
	
	# Find the mean for each bin
	y /= 100 * y_itrn
	# ~ print(y)
	
	# Add data to plot array
	for h in range(len(y)):
		z[h][day_itrn] += y[h]
		print(day_itrn)
	
	# Sort out date array
	# ~ date_time[day_itrn] = date_stamp
	# ~ date_time = np.append(date_time, date_stamp)
	date_time = np.append(date_time, nc.date2num(date_stamp, time_units, calendar='standard'))
	day_itrn += 1
	
	# Reinitialise orbit arrays
	y1 = np.copy(y2)
	y1_itrn = np.copy(y2_itrn)
	y2 = np.zeros(len(alts))
	y2_itrn = np.zeros(len(alts))
	
	orbittrigger = 0
	# ~ print(z)

nc_dates = np.copy(date_time[:-1])
# ~ date_time = coda.time_to_utcstring(date_time)
# ~ date_time = np.array([datetime.strptime(date,
			# ~ '%Y-%m-%d %H:%M:%S.%f') for date in date_time])
# ~ print(date_time)
x, y = np.meshgrid(date_time, alts)

print("x: ", x)
print("y: ", y)
print("z: ", z)
print("Shape of z: ", np.shape(z))
# ~ np.savetxt("neildata.txt", z)
# ~ savemat("neildata.mat", {'data':z}) # <-- This is the one you need

# ~ file2write = open('filename.txt', 'w') # <-- This really didn't work
# ~ file2write.write(z)
# ~ file2write.close()

# Creating netCDF file


root = nc.Dataset('timdata.nc', 'w', format = "NETCDF4")
root.contact = "T. P. Banyard, tpb38@bath.ac.uk"
root.institution = \
"University of Bath, Claverton Down, Bath, BA2 7AY, United Kingdom"
root.title = "Daily mean of the Aeolus zonal wind projection over the Equator"
root.Aeolus_data_source = "https://aeolus-ds.eo.esa.int"
dim_time = root.createDimension("time", len(nc_dates))
dim_alt = root.createDimension("altitude", len(alts))

var_time = root.createVariable("time", "f8", ("time",))
var_time.standard_name = "time"
var_time.long_name = "time"
var_time.units = time_units

var_alt = root.createVariable("altitude", "f8", ("altitude",))
var_alt.standard_name = "altitude"
var_alt.long_name = "altitude"
var_alt.units = "m"

var_u_proj = root.createVariable("Zonal_wind_projection", "f8", ("altitude", "time",))
var_u_proj.standard_name = "zonal_wind_projection"
var_u_proj.long_name = "Zonal projection of the HLOS wind"
var_u_proj.units = "cm s-1"

print("Shape of var_u_proj: ", np.shape(var_u_proj))
print("Shape of nc_dates: ", np.shape(nc_dates))
print("Shape of var_time: ", np.shape(var_time))
print("Shape of alts: ", np.shape(alts))
print("Shape of var_alt: ", np.shape(var_alt))
var_time[:], var_alt[:], var_u_proj[:,:] = nc_dates, alts, z

print("Created file of type: ", root.data_model)
root.close()


"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')

fig = plt.figure()
ax1 = fig.add_subplot(111)
cs = plt.contourf(x,y,z)
plt.savefig("hi.png", dpi=300)
plt.close()

# Plotting data
fig = plt.figure()
ax1 = fig.add_subplot(111)
cs = plt.contourf(x,y,z, cmap='RdBu_r',
	levels=np.linspace(-20, 20, 41))
# ~ ax2 = ax1.twinx()
# ~ ax2.plot(rayleigh_times, data_lat, c='black',
	# ~ marker='.', markersize='1', label='Latitude', linewidth=0.1)

# Setting Date axis
date_form = '%b %Y'
# ~ date_form = '%d %b'
minor_date_form = '%H'
hours = dates.HourLocator()
minutes = dates.MinuteLocator()
date_form = dates.DateFormatter(date_form)
minor_date_form = dates.DateFormatter(minor_date_form)
ax1.xaxis.set_major_formatter(date_form)
ax1.xaxis.set_minor_formatter(minor_date_form)
ax1.set_xlim(date_time[0], date_time[-2])
ax1.set_xlabel('Time')

# Setting y axes
ax1.set_ylabel('Altitude / m')
# ~ ax2.set_ylabel('Latitude / $^\circ$')
plt.title('Aeolus Orbit HLOS Rayleigh Wind Cross-section')
fig.colorbar(cs, cmap='RdBu_r', ax=ax1, orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms-1')
# ~ plt.legend(loc=9)
pngsavename = 'file.png'
plt.savefig(pngsavename,dpi=300)

# Climb out of plot directory
os.chdir('..')
os.chdir('..')

# Time taken for the program
pduration = datetime.now() - pstartTime
print('That program took ', pduration, ' seconds')
