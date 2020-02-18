#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Ascending/Descending Node Split-------------------------------
---v1.2---Daily_Means---------------------------------------------------
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
strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC/'
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
z = np.zeros((len(alts), day_itrn))
date_time = np.zeros(day_itrn)
daynum = 1 # Initial day number
day_itrn = 0
# Initialise orbit array
y = np.zeros(len(alts))
y_itrn = np.zeros(len(alts))

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
	data_HLOS = data.variables['HLOS_wind_speed'][:]
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
		# Is element within equatorial band?
		if lat_band[t] == 1:
			
			# Find the nearest altitude level
			val = find_nearest(alts, data_alt[t])
			alt_elmnt = np.where(alts == val)[0][0]
			
			# Cap wind speeds to 250 m/s
			if np.abs(data_HLOS[t]) < 250000:
				
				# Ascending or Descending node?	
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
					y_itrn[alt_elmnt] += 1
		
	# Primitive Daily Mean Test
	orbitdaynum = int(str(filename)[14:16])
	if daynum == orbitdaynum:
		continue
	
	daynum = orbitdaynum
	
	# Find the mean for each bin
	y /= 100 * y_itrn
	# ~ print(y)
	
	# Add data to plot array
	for h in range(len(y)):
		z[h][day_itrn] += y[h]
	
	mid_index = int(np.floor(len(lat_band)/2))
	date_time[day_itrn] += rayleigh_times[mid_index]
	
	day_itrn += 1
	# Reinitialise orbit array
	y = np.zeros(len(alts))
	y_itrn = np.zeros(len(alts))
	
	print(z)

date_time = coda.time_to_utcstring(date_time)
date_time = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
print(date_time)
x, y = np.meshgrid(date_time, alts)

# ~ np.savetxt("neildata.txt", z)
savemat("neildata.mat", {'data':z})


file2write = open('filename.txt', 'w')
file2write.write(z)
file2write.close()



root = nc.Dataset('file.nc', 'w', format = "NETCDF4")
dim_time = root.createDimension("time", len(rayleigh_times))
dim_alt = root.createDimension("altitude", len(alts))
var_time = root.createVariable("time", "f8", ("time",))
var_time.standard_name = "time"
var_time.long_name = "time"
var_time.units = time_units

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('..')
print(os.getcwd())
os.chdir('Plots')

plt.figure()
ax1 = fig.add_subplot(111)
cs = plt.contourf(x,y,z)
plt.savefig("hi.png", dpi=300)

"""# Plotting data
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
print('That program took ', pduration, ' seconds')"""
