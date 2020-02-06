#!/usr/bin/env python3
"""
Aeolus timeseries testing file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .DBL files downloaded from the Aeolus database and produces a
timeseries of the different variables available.

This program uses the timeseriesplot function from phdfunctions.py,
reads in data from a directory (e.g.'/home/tpb38/PhD/Bath/Aeolus/DATA/')
and plots two variables as a time series on different y axes (Y1 & Y2).
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
from datetime import datetime
# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
from phdfunctions import timeseriesplot

dbl = '/home/tpb38/PhD/Bath/Aeolus/DATA/'
hdr = dbl # Storing hdr and dbl files in same directory for now
dbl += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.DBL'
hdr += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.HDR'

# Opening Files
pf = coda.open(dbl)
pfhdr = coda.open(hdr)

# Fetching Data
mie_wind_velocity = coda.fetch(pf, 'mie_hloswind', -1,
	'windresult/mie_wind_velocity')
latitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/latitude_cog')
longitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/longitude_cog')
altitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/altitude_vcog')
date_time = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/datetime_cog')

# Retrieving Field Names
field_names = coda.get_field_names(pf)
mie_hloswind_field_names = coda.get_field_names(pf,
	'mie_hloswind', 0)
mie_geolocation_field_names = coda.get_field_names(pf,
	'mie_geolocation', 0)
windresult_geolocation_field_names = coda.get_field_names(pf,
	'mie_geolocation', 0, 'windresult_geolocation')
print(field_names)
print(mie_hloswind_field_names)
print(mie_geolocation_field_names)
print(windresult_geolocation_field_names)

# Arrays
print(mie_wind_velocity.shape)
print(mie_wind_velocity)

# Plotting data
os.chdir('..')
os.chdir('..')
os.chdir('Plots')

# Creating plot
X = date_time[11500:21500]	#[11500:21500] [11500:12500] [11540:11851]
Y = longitude[11500:21500]
Y2 = latitude[11500:21500]
variable = 'Longitude'
variable2 = 'Latitude'
timeseriesplot(X, Y, Y2,
plottitle = 'Variation in Aeolus\' longitude with latitude from 30N-30S',
	date_form = '%H:%M:%S', minor_date_form = '%M:%S', data_type = 'coda',
	size = 0.5, color = 'black', marker = '+', variable = variable,
	variable2 = variable2, legend = 0, l_adj = 0.15, r_adj=0.85)
# Saving plot
plt.savefig("timeseries.png", dpi=300)
print(coda.time_to_string(date_time[21000]))

