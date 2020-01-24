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
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import timeseriesplot, find_nearest
from functions import ncload


# Change current working directory to parent directory
os.chdir('..')

# Here I need to iterate through all. nc files and plot all of them into jpgs to view one after another
"""Find directory and read netCDF data"""
infile = '/home/tpb38/PhD/Bath/Aeolus/NC/'
infile += 'AE_2020-01-17_093235.nc' # Specifies ERA5 data file
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

# Converted time
data_time = nc.num2date(data.variables['time'][:],\
calendar = 'standard', units = data.variables['time'].units)

"""=================================================================="""
"""===========================Plotting==============================="""
"""=================================================================="""
os.chdir('..')
os.chdir('Plots')

YYYY = '2020'
MM = '01'

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
ax = fig.add_subplot(111)
#...continue code here to plot figures...
