#!/usr/bin/env python3
"""
Aeolus scans testing file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
----------[DEPRECATED]-There_is_a_newer_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .DBL files downloaded from the Aeolus database and produces a
graphical plot of the satellite pass.
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
from datetime import datetime
# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs')
from phdfunctions import timeseriesplot
from functions import load_hdr_tags, load_dbl_tags
from scipy.interpolate import griddata

dbl = '/home/tpb38/PhD/Bath/Aeolus/DATA/'
hdr = dbl # Storing hdr and dbl files in same directory for now
dbl += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.DBL'
hdr += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.HDR'

# Opening Files
pf = coda.open(dbl)
pfhdr = coda.open(hdr)

# Load tags
load_hdr = False
load_dbl = True

# Fetching Data
mie_wind_velocity = coda.fetch(pf, 'mie_hloswind', -1, 'windresult/mie_wind_velocity')
latitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/latitude_cog')
longitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/longitude_cog')
altitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/altitude_vcog')
date_time = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/datetime_cog')

# Retrieving Field Names
field_names = coda.get_field_names(pf)
mie_hloswind_field_names = coda.get_field_names(pf, 'mie_hloswind', 0)
mie_geolocation_field_names = coda.get_field_names(pf, 'mie_geolocation', 0)
windresult_geolocation_field_names = coda.get_field_names(pf, 'mie_geolocation', 0, 'windresult_geolocation')
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

X = date_time[11500:12500]	#[11500:21500] [11500:12500]
Y = altitude[11500:12500]
Y2 = latitude[11500:12500]
Y3 = longitude[11500:12500]
Y4 = mie_wind_velocity[11500:12500]
variable = 'Altitude'
variable2 = 'Latitude'
timeseriesplot(X, Y3, Y2, plottitle = 'Variation in Aeolus\' longitude with latitude during a two minute period', date_form = '%H:%M', minor_date_form = '%M:%S', data_type = 'coda', size = 0.5, color = 'black', marker = '+', variable = variable, variable2 = variable2, legend = 0, l_adj = 0.15, r_adj=0.85)
plt.savefig("timeseries.png", dpi=300)
# ~ print(coda.time_to_string(date_time[21000]))

# Initialise arrays
alts = []
lats = []
lons = []
mwvs = []
alts2 = []
lats2 = []
lons2 = []
mwvs2 = []

# Separating into a single level
for t in range(len(X)):
	if Y[t]>9900 and Y[t]<10000:
		alts.append(Y[t])
		lats.append(Y2[t])
		lons.append(Y3[t])
		mwvs.append(Y4[t])

for t in range(len(X)):
	if Y[t]>8900 and Y[t]<9000:
		alts2.append(Y[t])
		lats2.append(Y2[t])
		lons2.append(Y3[t])
		mwvs2.append(Y4[t])

x = lats
y = lons
z = mwvs
x2 = lats2
y2 = lons2
z2 = mwvs2

xi = np.linspace(30, -30, 1000)
yi = np.linspace(342, 330, 1000)

# ~ xi = np.linspace(28.5, 26.5, 1000)
# ~ yi = np.linspace(341.6, 341.1, 1000)

zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
zi2 = griddata((x2, y2), z2, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()
ax = fig.add_subplot(111)
# ~ cs = ax.contour(yi, xi, zi, 5, linewidths = 0.1, colors = 'k')
# ~ cs = ax.contourf(yi, xi, zi, 5, cmap='jet')
# ~ ax.contourf(x, y, z, cmap='RdBu', shading='interp')
# ~ plt.savefig("test23.png", dpi=300)

ax.scatter(y, x, s=5, c='blue', marker='.')
ax.scatter(y2, x2, s=5, c='red', marker='.')
ax.set_xlim(340, 341)
ax.set_ylim(20, 26)
plt.savefig("latlon2.png", dpi=300)


