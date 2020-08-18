#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Trimmed to look at the four equator passes of Jun 3rd---------
---v1.2---Completed subplots with profile, fixing code------------------
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
    'data_HLOS_wind': ds.data_vars['Rayleigh_HLOS_wind_speed'][:],
    'data_QC_Flag_Both': ds.data_vars['QC_Flag_Both'][:]
    }

# Create a binary array which restricts to an equatorial band between 15N-15S    
lat_band = xr.where(data['data_lat']<-10, 0, (xr.where(data['data_lat']>10, 0, 1)))

# Applying QCs
data['data_u_proj'] = xr.where(data['data_QC_Flag_Both'] == 0, -99999999, (data['data_u_proj']))

# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, -99999999, (data['data_u_proj']))
# np.set_printoptions(threshold=sys.maxsize)
# print(data['data_u_proj_capped'].values)

# Now restrict all arrays to be None (or -99999999) outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, -99999999, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, -99999999, data['data_lat'])
data['data_lon_band'] = xr.where(lat_band==0, -99999999, data['data_lon'])
data['data_alt_band'] = xr.where(lat_band==0, -99999999, data['data_alt'])

# Labelling of some important arrays for testing
u_proj_2 = data['data_u_proj_capped'][lat_band.values]
lat_1 = data['data_lat'][:]
lat_2 = data['data_lat'][lat_band.values]

# Slicing to leave only final few orbits
data['data_time'] = data['data_time'][810000:]
data['data_lat_band'] = data['data_lat_band'][810000:]
data['data_lon_band'] = data['data_lon_band'][810000:]
data['data_alt_band'] = data['data_alt_band'][810000:]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][810000:]

# Single pass
data['data_time'] = data['data_time'][5000:15000]
data['data_lat_band'] = data['data_lat_band'][5000:15000]
data['data_lon_band'] = data['data_lon_band'][5000:15000]
data['data_alt_band'] = data['data_alt_band'][5000:15000]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][5000:15000]

# Test using xr.where to split into altitude segments
    
# Binning into actual range bins
bin_edges = [-104, 400, 1408, 2416, 3424, 4432, 5440, 6448, 7456, 8464, 9472, 10480, 11488, 12496, 13504, 14764, 16024, 17284, 18544, 19804, 21064, 22576, 24088, 25600]
bin_mids = np.zeros(len(bin_edges)-1)
for num in range(len(bin_mids)):
    bin_mids[num] = np.floor((bin_edges[num+1] + bin_edges[num]) / 2)

# Creating fill variables
p = np.zeros(len(bin_mids))
j_itrn = 0

# Binning into altitude bins
for k in bin_mids:
    btm = k-500
    top = k+500
    new_winds = xr.where((data['data_alt_band'] < top) & (data['data_alt_band'] > btm), data['data_u_proj_capped_band'], None)
    new_winds = xr.where(new_winds == -99999999, None, new_winds)
    p[j_itrn] = new_winds.mean()
    j_itrn += 1
    
# Create numpy copies of orbit segment
np_time = np.copy(data['data_time'].values[:])
np_lat_band = np.copy(data['data_lat_band'].values[:])
np_lon_band = np.copy(data['data_lon_band'].values[:])
np_alt_band = np.copy(data['data_alt_band'].values[:])
np_u_proj = np.copy(data['data_u_proj_capped_band'].values[:])

# Resetting NaNs for plotting
data['data_alt_band'] = xr.where(data['data_alt_band'] == -99999999, None, data['data_alt_band']) 

# Plotting Imports
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
# ax1.add_feature(cfeature.OCEAN)
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

# Sorting geographical plot ticks
ax1.tick_params(bottom=True,top=True,left=True,right=True,
            labelbottom=True,labeltop=False,labelleft=True,labelright=False)
ax1.xaxis.set_visible(True)
ax1.yaxis.set_visible(True)

ax1.set_aspect('auto') # Stretch map to fill subplot

# Profile Plot
ax2 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
ax2.scatter(data['data_time'].values[:], data['data_alt_band'].values[:], marker='x', s=0.15, color='blue', zorder=1)
ax2.set_ylabel('Altitude / m')
ax2.set_xlabel('Time / HH:MM')
date_form = dates.DateFormatter('%H:%M') # Sets date format
ax2.xaxis.set_major_formatter(date_form)
# ax1.xaxis.set_major_locator(plt.MaxNLocator(4)) # Maximum number of date ticks

ax3 = ax2.twiny()
# ax3 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
ax3.plot(p[1:], bin_mids[1:], color='red')
# ax3.set_xticks([-25,-20,-15,-10,-5,0,5])

# Finishing figure
fig.suptitle("Aeolus Data for Singapore pass on 3rd June 2020 $\pm$10$^\circ$ about equator")
plt.savefig("testtesttest.png", dpi=300)
print(os.getcwd())
