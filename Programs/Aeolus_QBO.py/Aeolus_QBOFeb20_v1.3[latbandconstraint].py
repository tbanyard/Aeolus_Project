#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Ascending/Descending Node Split-------------------------------
---v1.2---Daily_Means and creating .mat file for Neil-------------------
---v1.3---Plotting my own version of this plot extended out to March----
----------[CURRENT]-This_is_the_current_version_of_this_file------------
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
import time
from scipy.interpolate import griddata
from scipy.io import savemat
from itertools import groupby

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import timeseriesplot, find_nearest
from functions import ncload

# Program Timing	
pstartTime = datetime.now()

# Change current working directory to parent directory
os.chdir('..')

# Find directory and read netCDF data
strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC_FullQC/'
directory = os.fsencode(strdirectory)
ncflist = [name for name in os.listdir(strdirectory)] # NC File list
daynum = 1
day_itrn = 1

for file in sorted(os.listdir(directory)):
	filename = os.fsdecode(file)
	orbitdaynum = int(str(filename)[14:16])
	if daynum != orbitdaynum:
		day_itrn += 1
		daynum = orbitdaynum
		
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
y1 = np.zeros(len(alts))
y1_itrn = np.zeros(len(alts))
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
	lat_band = np.where(data_lat<-2, 0, (np.where(data_lat>2, 0, 1)))
	
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
				date_stamp = date_stamp.replace(hour=0)
				date_stamp = date_stamp.replace(minute=0)
				date_stamp = date_stamp.replace(second=0)
				date_stamp = date_stamp.replace(microsecond=0)
				daynum = now_day # Set the day number to be the current day
				monthnum = now_month # Set the month number to be the current month
				daydt = nowstrp
			
		# Is element within equatorial band?
		if lat_band[t] == 1:
			
			# Find the nearest altitude level
			val = find_nearest(alts, data_alt[t])
			alt_elmnt = np.where(alts == val)[0][0]
			
			# Cap wind speeds to 250 m/s
			if np.abs(data_u_proj[t]) < 250000:
				
				if currnext == 0:
					y1[alt_elmnt] += data_u_proj[t]
					y1_itrn[alt_elmnt] += 1
				
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
	y1 /= 100 * y1_itrn
	print(y1)
	print(y1_itrn)
	
	# Add data to plot array
	for h in range(len(y1)):
		z[h][day_itrn] += y1[h]
	
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


root = nc.Dataset('timlatband.nc', 'w', format = "NETCDF4")
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

x = x[:,:-1]

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
pngsavename = 'file02.png'
plt.savefig(pngsavename,dpi=300)

# Climb out of plot directory
os.chdir('..')
os.chdir('..')

# Time taken for the program
pduration = datetime.now() - pstartTime
print('That program took ', pduration, ' seconds')
