#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Trimmed to look at the four equator passes of Jun 3rd---------
---v1.2---Completed subplots with profile, fixing code------------------
---v1.3---Focusing on one or two Singapore profiles---------------------
---v1.4---Aug2020-Reusing code for different radiosonde launches--------
---v_5km+---Only plotting from 5km upwards, corrected altitude/pressure-
----------[BRANCH]-This_is_a_working_branch_version_of_this_file--------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
colocates with one or two Singapore radiosonde profiles for testing
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
ncdir = '/home/tim/Documents/Bath/NC_AEOLUS_DATA_AUG'

"""============================Aeolus Data==========================="""
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

# Create a binary array which restricts to a latitude band between 0N-2.5N
lat_band = xr.where(data['data_lat']<0.5, 0, (xr.where(data['data_lat']>2, 0, 1)))
# Do the same for a longitude range between 100E-105E (Mainly to separate orbits)
lon_band = xr.where(data['data_lon']<103, 0, (xr.where(data['data_lon']>104, 0, 1)))

# Applying QCs
data['data_u_proj'] = xr.where(data['data_QC_Flag_Both'] == 0, -9999, (data['data_u_proj']))

# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, -9999, (data['data_u_proj']))

# Toggle additional check to limit data to within 5 m/s of Radiosonde extremities
data['data_u_proj_capped'] = xr.where(data['data_u_proj_capped'] > 2000, -9999, (xr.where(data['data_u_proj_capped'] < -3000, -9999, data['data_u_proj_capped'])))

# Now restrict all arrays to be None (or -99999999) outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, -9999, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, -9999, data['data_lat'])
data['data_lon_band'] = xr.where(lat_band==0, -9999, data['data_lon'])
data['data_alt_band'] = xr.where(lat_band==0, -9999, data['data_alt'])

# Do the same outside of the longitude band
data['data_u_proj_capped_band'] = xr.where(lon_band==0, -9999, data['data_u_proj_capped_band'])
data['data_lat_band'] = xr.where(lon_band==0, -9999, data['data_lat_band'])
data['data_lon_band'] = xr.where(lon_band==0, -9999, data['data_lon_band'])
data['data_alt_band'] = xr.where(lon_band==0, -9999, data['data_alt_band'])

printfullarrays()
print(data['data_lon_band'].values)

# Labelling of some important arrays for testing
u_proj_2 = data['data_u_proj_capped'][lat_band.values]
lat_1 = data['data_lat'][:]
lat_2 = data['data_lat'][lat_band.values]

# Slicing to leave only final few orbits
data['data_time'] = data['data_time'][:]
data['data_lat_band'] = data['data_lat_band'][:]
data['data_lon_band'] = data['data_lon_band'][:]
data['data_alt_band'] = data['data_alt_band'][:]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][:] # [810000:]

# Single pass
data['data_time'] = data['data_time'][:] # [12657:12738]
data['data_lat_band'] = data['data_lat_band'][:]
data['data_lon_band'] = data['data_lon_band'][:]
data['data_alt_band'] = data['data_alt_band'][:]
data['data_u_proj_capped_band'] = data['data_u_proj_capped_band'][:]

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

"""======================Singapore Radiosonde Data==================="""
# Read in Singapore Radiosonde Data
# SR_pressure = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 0]
# SR_wind_dir = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 1]
# SR_wind_spd = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data.txt")[:, 2]
# SR_u_wind = SR_wind_spd * np.sin(SR_wind_dir*np.pi/180)

# Correct Radiosonde Data
SR_pressure = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data_2020-08-13_0z.txt")[:, 2]
SR_wind_dir = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data_2020-08-13_0z.txt")[:, 7]
SR_wind_spd = np.loadtxt("/home/tim/Documents/Bath/SNM00048698/SNM00048698-wind-data_2020-08-13_0z.txt")[:, 8]
SR_u_wind = SR_wind_spd * np.sin(SR_wind_dir*np.pi/180) # Calculate u wind

print("SR_pressure: ", SR_pressure)
print("SR_wind_dir: ", SR_wind_dir)
print("SR_wind_spd: ", SR_wind_spd)

"""============================Plotting=============================="""
# Plotting Imports
import geopandas
import geoplot
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib

"""+++++++++++++++++++++++Geographical Map Plot++++++++++++++++++++++"""
fig = plt.figure()
gridspec.GridSpec(5,5)
ax1 = plt.subplot2grid((5,5), (0,0), colspan=1, rowspan=5, projection=ccrs.PlateCarree(103.5))
# ax1 = plt.subplot(111, projection=ccrs.PlateCarree(100.0))
# ax1.coastlines(resolution='50m', linewidth=0.5)
ax1.set_extent([103, 104, 0, 2], ccrs.PlateCarree())
ax1.add_feature(cfeature.OCEAN)
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.LAND)
# plt.axes(projection=ccrs.PlateCarree(100.0))

# ax1.set_aspect('equal')
# world = geopandas.read_file(
    # geopandas.datasets.get_path('naturalearth_lowres')
# )
# geoplot.polyplot(world, extent=(100, -10, 105, 10), edgecolor='black')

# Aeolus data positions on map
ax1.scatter(data['data_lon_band'].values[:], data['data_lat_band'].values[:], marker='.', s=15, linewidths = 0.5, edgecolor = 'k', color='red', zorder=1, transform=ccrs.PlateCarree())
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
ax1.text(103.8, 1.35, 'x', color = 'green', transform=ccrs.PlateCarree())

ax1.set_aspect('auto') # Stretch map to fill subplot

"""+++++++++++++++Radiosonde vs Aeolus profile plot++++++++++++++++++"""
ax2 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)

# Convert m to km for altitude
printfullarrays()
# data['data_alt_band'] = np.where(data['data_alt_band']==None, None, data['data_alt_band'].values/1000.0)
values = data['data_alt_band'].values
for i in range(len(data['data_alt_band'])):
    if str(type(data['data_alt_band'].values[i]))[8:16] == 'NoneType':
        values[i] = None
    else:
        values[i] = data['data_alt_band'].values[i]/1000.0
bin_mids /= 1000

# Convert cm/s to m/s for wind speed for Aeolus data
p /= 100

# Plotting mean Aeolus profile
l_AE = ax2.plot(p[1:], bin_mids[1:], color='red', label='Aeolus')
ax2.set_ylabel('Altitude / km')
ax2.set_yticks([5, 10, 15, 20, 25])
ax2.set_xlabel('Zonal Wind Speed / ms$^{-1}$')
ax2.set_xticks([-40,-30,-20,-10,0,10,20])
ax2.set_xticklabels([-40,-30,-20,-10,0,10,20])
ax2.set_xlim([-40,20])
ax2.set_ylim([4.5, 25.5])
# ax2.grid(which='major', color = 'red')

# ax3 = ax2.twiny() # Toggle
# ax3 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=5)
# ax3.scatter(data['data_time'].values[:], values, marker='+', s=1, color='blue', zorder=1) # Toggle
date_form = dates.DateFormatter('%H:%M:%S') # Sets date format
# ax3.xaxis.set_major_formatter(date_form) # Toggle
# ax3.xaxis.set_major_locator(plt.MaxNLocator(4)) # Maximum number of date ticks # Toggle
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

ax4 = ax2.twinx()
l_RS = ax4.plot(winds.values[:], SR_pressure/100, color='green', label='0z Singapore Radiosonde')
ax4.invert_yaxis()
ax4.set_yscale('log')
ax4.set_ylabel('Pressure / hPa', rotation=270, labelpad=10)
ax4.set_yticks([500, 250, 100, 50, 30, 20, 10])
ax4.set_yticklabels([500, 250, 100, 50, 30, 20, 10])
ax4.set_ylim([599,24.5])
ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax4.grid(which='major', color = 'green')

l_all = l_AE + l_RS
labs = [l.get_label() for l in l_all]
ax2.legend(l_all, labs, loc = 'upper left', fontsize = 'small', framealpha=0, frameon=False)

# Finishing figure
fig.suptitle("Singapore Radiosonde validation for Aeolus pass: 12th August 2020 22:58")
plt.savefig("testtesttest3.png", dpi=300)
print(os.getcwd())
