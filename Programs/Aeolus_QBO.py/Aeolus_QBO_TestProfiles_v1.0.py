#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---[DEPRECATED]-There_is_a_newer_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
generates zonal mean U plot at the equator
========================================================================
"""
########################################################################
# Please note, for brevity, most imports are not shown here. Please read
# the README.md file and the .pythonstartup file for more information.
# Imports required for this specific program can be found in imports.py

# Imports
from pathlib import Path
home = str(Path.home()) # Writes the home directory as a string
pythonstartup = home + '/.pythonstartup' # Or any alternative directory
exec(open(pythonstartup).read()) # Reads .pythonstartup file
imports = os.getcwd() + '/imports.py'
exec(open(imports).read()) # Reads imports file

# Program Timing
pstartTime = datetime.now()

# Change current working directory to parent directory
os.chdir('..') # Partly so that plots can be saved in their own folder
########################################################################

"""=================================================================="""
"""================USER PREFERENCES (Toggles/Variables)=============="""
"""=================================================================="""
# Please uncomment the required program mode. DEFAULT = Compute.
mode = 'compute'
# mode = 'plot'

# Replace the below with the directory containing the AE netCDF files.
ncdir = '/home/tim/Documents/Bath/NC_AEOLUS_DATA_JUN'

# Find and read netCDF data
ds = xr.open_mfdataset(ncdir + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
# Create a function which places data variables into a dictionary
### data_dict::phdfunctions
data = {
    'data_time': ds.coords['time'][:],
    'data_lat': ds.data_vars['lat'][:],
    'data_lon': ds.data_vars['lon'][:],
    'data_alt': ds.data_vars['alt'][:],
    'data_u_proj': ds.data_vars['Zonal_wind_projection'][:],
    }

# Create a binary array which restricts to an equatorial band between 15N-15S    
lat_band = xr.where(data['data_lat']<-10, 0, (xr.where(data['data_lat']>10, 0, 1)))

# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, None, (data['data_u_proj']))

# Now restrict all arrays to be None outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, None, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, None, data['data_lat'])
data['data_alt_band'] = xr.where(lat_band==0, None, data['data_alt'])

# Labelling of some important arrays for testing
u_proj_2 = data['data_u_proj_capped'][lat_band.values]
lat_1 = data['data_lat'][:]
lat_2 = data['data_lat'][lat_band.values]
print("lat_1: ", lat_1.values)
print("lat_2: ", lat_2.values)
print("==/==")
print("data_time: ", data['data_time'])

# Slicing to leave only final few orbits
# print("len: ", len(data['data_u_proj_capped_band']))
# data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][810000:]
# print("len: ", len(data['data_u_proj_capped_band']))
# data['data_time'] = data['data_time'][810000:]
# data['data_alt_band'] = data['data_alt_band'][810000:]

# pippin = data['data_time'].sel(time=slice(0,1))
# print("pippin: ", pippin)

# Plotting
plt.figure()
ax1 = plt.subplot(111)
ax1.scatter(data['data_time'].values[:], data['data_alt_band'].values[:], marker='x', s=0.25, color='black')
# plt.scatter(lat_1.values[:], lat_band.values[:], marker='x', s=0.5, color='black')
ax1.set_ylabel('Altitude / m')
ax1.set_xlabel('Date')
# date_form = dates.DateFormatter('%d %b') # Sets date format
# ax1.xaxis.set_major_formatter(date_form)
# ax1.xaxis.set_major_locator(plt.MaxNLocator(4)) # Maximum number of date ticks
plt.title("Aeolus Data from 1-3 Jun $\pm$10$^\circ$ about equator")
plt.savefig("testtesttest.png", dpi=300)
print(os.getcwd())
