#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Trimmed to look at the four equator passes of Jun 3rd---------
---v1.2---Completed subplots with profile, fixing code------------------
---v1.3---Focusing on one or two Singapore profiles---------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
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
data['data_u_proj'] = xr.where(data['data_QC_Flag_Both'] == 0, -9999, (data['data_u_proj']))

# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, -9999, (data['data_u_proj']))

# Now restrict all arrays to be None (or -99999999) outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, -9999, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, -9999, data['data_lat'])
data['data_lon_band'] = xr.where(lat_band==0, -9999, data['data_lon'])
data['data_alt_band'] = xr.where(lat_band==0, -9999, data['data_alt'])

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
data['data_time'] = data['data_time'][12657:12738]
data['data_lat_band'] = data['data_lat_band'][12657:12738]
data['data_lon_band'] = data['data_lon_band'][12657:12738]
data['data_alt_band'] = data['data_alt_band'][12657:12738]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][12657:12738]

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
    new_winds = xr.where(new_winds == -9999, None, new_winds)
    p[j_itrn] = new_winds.mean()
    j_itrn += 1
    
# Create numpy copies of orbit segment
np_time = np.copy(data['data_time'].values[:])
np_lat_band = np.copy(data['data_lat_band'].values[:])
np_lon_band = np.copy(data['data_lon_band'].values[:])
np_alt_band = np.copy(data['data_alt_band'].values[:])
np_u_proj = np.copy(data['data_u_proj_capped_band'].values[:])

# Resetting NaNs for plotting
data['data_alt_band'] = xr.where(data['data_alt_band'] == -9999, None, data['data_alt_band']) 
# data['data_u_proj_capped_band'] = xr.where(data['data_u_proj_capped_band'] == -99999999, None, data['data_u_proj_capped_band'])

# Read in Singapore Radiosonde Data
SR_pressure = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 0]
SR_wind_dir = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 1]
SR_wind_spd = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 2]
SR_u_wind = SR_wind_spd * np.sin(SR_wind_dir*np.pi/180)

# Correct Radiosonde Data
# SR_pressure = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-0zwind-data.txt")[:, 2]
# SR_wind_dir = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-0zwind-data.txt")[:, 7]
# SR_wind_spd = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-0zwind-data.txt")[:, 8]
# SR_u_wind = SR_wind_spd * np.sin(SR_wind_dir*np.pi/180)


print("SR_pressure: ", SR_pressure)
print("SR_wind_dir: ", SR_wind_dir)
print("SR_wind_spd: ", SR_wind_spd)

# Plotting Imports
import geopandas
import geoplot
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib

# Plotting
fig = plt.figure()
gridspec.GridSpec(5,5)
ax1 = plt.subplot2grid((5,5), (0,0), colspan=1, rowspan=5, projection=ccrs.PlateCarree(103.5))
# ax1 = plt.subplot(111, projection=ccrs.PlateCarree(100.0))
ax1.coastlines(resolution='10m', linewidth=1)
ax1.set_extent([103, 104, 0, 2], ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN)
ax1.add_feature(cfeature.COASTLINE)
# plt.axes(projection=ccrs.PlateCarree(100.0))

# ax1.set_aspect('equal')
# world = geopandas.read_file(
    # geopandas.datasets.get_path('naturalearth_lowres')
# )
# geoplot.polyplot(world, extent=(100, -10, 105, 10), edgecolor='black')

ax1.scatter(data['data_lon_band'].values[:], data['data_lat_band'].values[:], marker='+', s=1, color='blue', zorder=1, transform=ccrs.PlateCarree())
# plt.scatter(lat_1.values[:], lat_band.values[:], marker='x', s=0.5, color='black')

# Latitude axis
ax1.set_ylabel('Latitude / $^\circ$')
ax1.set_yticks([0, 0.5, 1, 1.5, 2], crs=ccrs.PlateCarree())
lat_formatter = LatitudeFormatter()
ax1.yaxis.set_major_formatter(lat_formatter)

# Longitude axis
ax1.set_xlabel('Longitude / $^\circ$')
ax1.set_xticks([103, 104], crs=ccrs.PlateCarree())
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

# Convert m to km for altitude
data['data_alt_band'] = xr.where(data['data_alt_band']==None, None, data['data_alt_band'].values/1000)
bin_mids /= 1000

# Convert cm/s to m/s for wind speed
p /= 100

ax2.plot(p[1:], bin_mids[1:], color='red')
ax2.set_ylabel('Altitude / km')
ax2.set_yticks([0, 5, 10, 15, 20, 25])
ax2.set_xlabel('Zonal Wind Speed / ms$^{-1}$')
ax2.set_xticks([-40,-30,-20,-10,0,10,20])
ax2.set_xticklabels([-40,-30,-20,-10,0,10,20])
ax2.set_xlim([-40,20])

ax3 = ax2.twiny()
# ax3 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
ax3.scatter(data['data_time'].values[:], data['data_alt_band'].values[:], marker='+', s=1, color='blue', zorder=1)
date_form = dates.DateFormatter('%H:%M:%S') # Sets date format
ax3.xaxis.set_major_formatter(date_form)
ax3.xaxis.set_major_locator(plt.MaxNLocator(4)) # Maximum number of date ticks
# ax3.set_xticks([-25,-20,-15,-10,-5,0,5])

SR_u_wind = np.where(SR_u_wind == -9875.895717610794, None, SR_u_wind/10)
SR_u_wind = np.asarray(SR_u_wind, dtype = 'float64')
print(SR_u_wind.dtype)

# Interpolate through nans
"""def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

nans, x = nan_helper(SR_u_wind)
SR_u_wind[nans] = np.interp(x(nans), x(~nans), SR_u_wind[~nans])"""

winds = xr.DataArray(-SR_u_wind, dims = 'pressure')
print(winds)
winds = winds.interpolate_na(dim = ('pressure'), method = 'linear')
print("winds: ", winds.values)

# ax4 = ax2.twinx()
# ax4.plot(winds.values[:], SR_pressure/100, color='green')
# ax4.invert_yaxis()
# ax4.set_yscale('log')
# ax4.set_ylabel('Pressure / hPa', rotation=270, labelpad=10)
# ax4.set_yticks([1000, 500, 250, 100, 50, 30, 20, 10])
# ax4.set_yticklabels([1000, 500, 250, 100, 50, 30, 20, 10])
# ax4.set_ylim([1500,21])
# ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# Finishing figure
fig.suptitle("Aeolus Data for Singapore pass on 3rd June 2020 $\pm$10$^\circ$ about equator")
plt.savefig("testtesttest2.png", dpi=300)
print(os.getcwd())
