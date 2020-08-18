#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Trimmed to look at the four equator passes of Jun 3rd---------
----------[DEPRECATED]-There_is_a_newer_version_of_this_file------------
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
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, -9999, (data['data_u_proj']))

# Now restrict all arrays to be None outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, -9999, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, -9999, data['data_lat'])
data['data_lon_band'] = xr.where(lat_band==0, -9999, data['data_lon'])
data['data_alt_band'] = xr.where(lat_band==0, -9999, data['data_alt'])

# Labelling of some important arrays for testing
u_proj_2 = data['data_u_proj_capped'][lat_band.values]
lat_1 = data['data_lat'][:]
lat_2 = data['data_lat'][lat_band.values]
print("lat_1: ", lat_1.values)
print("lat_2: ", lat_2.values)
print("==/==")
print("data_time: ", data['data_time'])

# Slicing to leave only final few orbits
data['data_time'] = data['data_time'][810000:]
data['data_lat_band'] = data['data_lat_band'][810000:]
data['data_lon_band'] = data['data_lon_band'][810000:]
data['data_alt_band'] = data['data_alt_band'][810000:]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][810000:]

print("lennn: ", len(data['data_time']))

# Single pass
data['data_time'] = data['data_time'][5000:15000]
data['data_lat_band'] = data['data_lat_band'][5000:15000]
data['data_lon_band'] = data['data_lon_band'][5000:15000]
data['data_alt_band'] = data['data_alt_band'][5000:15000]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][5000:15000]

# pippin = data['data_time'].sel(time=slice(0,1))
# print("pippin: ", pippin)

# Test using xr.where to split into altitude segments
p = np.zeros(26)
j_itrn = 0

# better1 = data['data_alt_band'].dropna('time')
# better2 = data['data_u_proj_capped_band'].dropna('time')
# better1 = data['data_alt_band'].fillna(-9999)
# better2 = data['data_u_proj_capped_band'].fillna(-9999)

for j in np.linspace(500, 25500, 26):
    new_winds = xr.where((data['data_alt_band'] < j) & (data['data_alt_band'] > j-1000), data['data_u_proj_capped_band'], 0)
    p[j_itrn] = np.nanmean(new_winds)
    j_itrn += 1
    
    
# new_winds = xr.where((data['data_alt_band'].values >= j-1000) & (data['data_alt_band'].values < j), data['data_u_proj_capped_band'], 0)

# Create numpy copies of orbit segment
np_time = np.copy(data['data_time'].values[:])
np_lat_band = np.copy(data['data_lat_band'].values[:])
np_lon_band = np.copy(data['data_lon_band'].values[:])
np_alt_band = np.copy(data['data_alt_band'].values[:])
np_u_proj = np.copy(data['data_u_proj_capped_band'].values[:])

alts = np.linspace(0,25000, 26)
x = np.zeros(len(alts))
x_itrn = np.zeros(len(alts))
tot=0

# Attempt at iterating through...
# for i in range(len(np_time)):
    # if np_alt_band[i] == None:
        # tot +=1
        # continue
    # val = find_nearest(alts, np_alt_band[i])
    # alt_elmnt = np.where(alts == val)[0][0]
    
    # x[alt_elmnt] += np_u_proj[i]
    # x_itrn[alt_elmnt] += 1
    
# x /= 100 * x_itrn
# print("tot: ", tot)

data['data_alt_band'] = xr.where(data['data_alt_band'] == -9999, None, data['data_alt_band']) 

import geopandas
import geoplot
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Plotting
fig = plt.figure()
gridspec.GridSpec(5,5)
ax1 = plt.subplot2grid((5,5), (0,0), colspan=1, rowspan=5, projection=ccrs.PlateCarree(100.0))
# ax1 = plt.subplot(111, projection=ccrs.PlateCarree(100.0))
ax1.coastlines(resolution='10m', linewidth=1)
ax1.set_extent([100, 105, -10, 10], ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN)
# plt.axes(projection=ccrs.PlateCarree(100.0))

# ax1.set_aspect('equal')
# world = geopandas.read_file(
    # geopandas.datasets.get_path('naturalearth_lowres')
# )
# geoplot.polyplot(world, extent=(100, -10, 105, 10), edgecolor='black')

ax1.scatter(data['data_lon_band'].values[:], data['data_lat_band'].values[:], marker='x', s=0.15, color='red', zorder=1, transform=ccrs.PlateCarree())
# plt.scatter(lat_1.values[:], lat_band.values[:], marker='x', s=0.5, color='black')

# Latitude axis
ax1.set_ylabel('Latitude / $^\circ$')
ax1.set_yticks([-10, -5, 0, 5, 10], crs=ccrs.PlateCarree())
lat_formatter = LatitudeFormatter()
ax1.yaxis.set_major_formatter(lat_formatter)

# Longitude axis
ax1.set_xlabel('Longitude / $^\circ$')
ax1.set_xticks([100, 105], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
ax1.xaxis.set_major_formatter(lon_formatter)


ax1.tick_params(bottom=True,top=True,left=True,right=True,
            labelbottom=True,labeltop=False,labelleft=True,labelright=False)
ax1.xaxis.set_visible(True)
ax1.yaxis.set_visible(True)

ax1.set_aspect('auto') # Stretch map to fill subplot

# date_form = dates.DateFormatter('%d %b') # Sets date format
# ax1.xaxis.set_major_formatter(date_form)
# ax1.xaxis.set_major_locator(plt.MaxNLocator(4)) # Maximum number of date ticks


ax2 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
ax2.scatter(data['data_time'].values[:], data['data_alt_band'].values[:], marker='x', s=0.15, color='blue', zorder=1)
ax2.set_ylabel('Altitude / m')
ax2.set_xlabel('Time / HH:MM')
date_form = dates.DateFormatter('%H:%M') # Sets date format
ax2.xaxis.set_major_formatter(date_form)

ax3 = ax2.twiny()
# ax3 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
ax3.plot(p[1:], alts[1:], color='red')
ax3.set_xticks([-25,-20,-15,-10,-5,0,5])

fig.suptitle("Aeolus Data for Singapore pass on 3rd June 2020 $\pm$10$^\circ$ about equator")
plt.savefig("testtesttest.png", dpi=300)
print(os.getcwd())
