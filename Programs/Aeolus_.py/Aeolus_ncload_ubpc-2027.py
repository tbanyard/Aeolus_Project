#!/usr/bin/env python3
"""
Aeolus data load from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
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
os.putenv('CODA_DEFINITION', '/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata

# Import from functions file
import sys
# ~ sys.path.append('/home/tpb38/PhD/Bath/')
# ~ sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
sys.path.append('/media/GWR/AEOLUS/')
from phdfunctions import timeseriesplot, find_nearest
from functions import ncload

# Change current working directory to parent directory
# ~ os.chdir('..')

# Here I need to iterate through all. nc files and plot all of them into pngs to view one after another
"""Find directory and read netCDF data"""
# ~ infile = '/home/tpb38/PhD/Bath/Aeolus/NC/'
# ~ infile += 'AE_2019-07-19_221623.nc' # Specifies data file
infile = '/media/GWR/AEOLUS/NC/'
infile += 'AE_2B_2019-04-17_055717.nc'
print('netCDF file:')
print(infile, '\n')
data = nc.Dataset(infile)

"""Download variables"""
# Longitude
data_lon = data.variables['lon'][:]
# Latitude
data_lat = data.variables['lat'][:]
# Altitude
data_alt = data.variables['alt'][:]
# Horizontal Line of Sight Wind speed
data_HLOS = data.variables['HLOS_wind_speed'][:]
# Rayleigh_Grouping goes here
RG = data.variables['RG'][:]
# Time
rayleigh_times = data.variables['time'][:]

# Converted time
# ~ data_time = nc.num2date(data.variables['time'][:],\
# ~ calendar = 'standard', units = data.variables['time'].units)

"""===="""

# Initialise meshgrids for x, y and z
alts = np.linspace(0,20000, 21)
x, y = np.meshgrid(RG, alts)
print(x)
print(y)
# ~ z = [[0 for _ in range(len(RG))] for _ in range(len(alts))] # Lists
# ~ z_itrn = [[0 for _ in range(len(RG))] for _ in range(len(alts))]
z = np.zeros((len(alts),len(RG))) # NumPy Arrays
z_itrn = np.zeros((len(alts),len(RG)))
print(np.shape(z))

# Placing wind values into bins of height 1km and width 1 rayleigh group
lastgroupstarttime = 0
for RG_elmnt in range(len(RG)):
	for t in range(len(rayleigh_times)):
		# Find all elements inside this sandwich and add to z and z_itrn:
		if rayleigh_times[t] < RG[RG_elmnt] and rayleigh_times[t] >= lastgroupstarttime:
			val = find_nearest(alts, data_alt[t]) # Find the nearest altitude level
			alt_elmnt = np.where(alts == val)[0][0]
			if np.abs(data_HLOS[t]) < 25000: # Cap wind speeds to 250 m/s
				z[alt_elmnt][RG_elmnt] += data_HLOS[t]
				z_itrn[alt_elmnt][RG_elmnt] += 1
	lastgroupstarttime = RG[RG_elmnt]

# Find the mean for each bin
z /= 100 * z_itrn # Factor of 100 for conversion from cm/s to m/s
print(z)

date_time = coda.time_to_utcstring(RG[:])
date_time = np.array([datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f') for date in date_time])
rayleigh_times = coda.time_to_utcstring(rayleigh_times[:])
rayleigh_times = np.array([datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times])
x, y = np.meshgrid(date_time, alts)
# ~ print(x)
# ~ print(y)

"""=================================================================="""
"""===========================Plotting==============================="""
"""=================================================================="""
# ~ os.chdir('..')
os.chdir('Plots')

YYYY = '2019'
MM = '04'

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

# Plotting data
fig = plt.figure()
ax1 = fig.add_subplot(111)
cs = plt.contourf(x,y,z, cmap='RdBu', levels=np.linspace(-200, 200, 41))
ax2 = ax1.twinx()
ax2.plot(rayleigh_times, data_lat, c='black', marker='.', markersize='1', label='Latitude', linewidth=0.1)

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
fig.colorbar(cs, cmap='RdBu', ax=ax1, orientation='horizontal', label='HLOS Rayleigh Wind Speed / ms-1')
# ~ plt.legend(loc=9)
plt.savefig('test29.1.20.png',dpi=300)
