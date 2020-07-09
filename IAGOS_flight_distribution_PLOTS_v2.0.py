#!/usr/bin/env python3
"""
IAGOS Flight distribution plotting file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Plotting code for a single .npy file of data------------------
---v1.2---Combine each month into one plot------------------------------
---v1.3---Neater plot with logarithmic scale and the option to----------
----------produce an orthographic projection of the data----------------
---v2.0---Tidied_up_code: Complete finalised version of code------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .npy files produced from IAGOS_flight_distribution_vx.x.py
and produces plots of IAGOS flight distribution.
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.colors as colors
import os
import errno
from mpl_toolkits.basemap import Basemap

# Change current working directory to parent directory
os.chdir('..')

# Set desired global variables
res = 0.625 # Resolution of density plot in degrees

# Set file directory for IAGOS data
infile = '/home/tpb38/PhD/Bath/IAGOS/Timeseries/' # IAGOS directory
itrn = 0 # Flight number iterator
resrcp = 1/res # Resolution Reciprocal
# ~ print('resrcp = ', resrcp)

# Check that the input resolution is acceptable
if 180 % res != 0:
	raise ValueError('Invalid resolution. Must be both a factor of 180 \
and non-residual in binary.')

# Calculate the number of latitude and longitude points
num_lat = int(180 * resrcp)+1
# ~ print('num_lat: ', num_lat)
num_lon = int(360 * resrcp)+1
# ~ print('num_lon: ', num_lon)

# Build lat/lon grid for data to be placed in
data = np.asarray([[0 for _ in np.linspace(-180,180,num_lon)] \
for _ in np.linspace(-90,90,num_lat)]) # itrn_array[latitude][longitude]

# Change directory to point to the correct data files
os.chdir('DATA_0.625deg.npy')

for YYYY in range(1994,2021):
	npsavename = str(YYYY) + '_flight_distributions.npy'
	itrn_array = np.load(npsavename, allow_pickle=True)
	# ~ print(itrn_array)
	data += itrn_array

os.chdir('..')

# Load lat/lon grids from file
data_lat = np.load("lat_grid.npy", allow_pickle=True)
data_lon = np.load("lon_grid.npy", allow_pickle=True)

"""=================================================================="""
"""===========================Plotting==============================="""
"""=================================================================="""
os.chdir('..')
os.chdir('Plots')

# Plotting data
fig = plt.figure()
ax = fig.add_subplot(111)
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='l', ax=ax)
# ~ map = Basemap(projection='ortho',lat_0=45,lon_0=0,\
            # ~ llcrnrlon=-180,urcrnrlon=180,resolution='i', ax=ax)
map.drawcoastlines(linewidth=0.5)
# draw parallels and meridians.
map.drawparallels(np.arange(-90.,91.,30.),linewidth=0.5)
map.drawmeridians(np.arange(-180.,180.,60.),linewidth=0.5)
lons, lats = np.meshgrid(data_lon, data_lat)
x,y = map(lons,lats)
cs = map.contourf(x,y,np.log(data), cmap='Reds',
levels=[0,1,2,3,4,5,6,7,8,9,10], vmin=0, vmax=10, extend="both")
# ~ cs = map.contourf(x,y,data==0, cmap='Greys', levels=[0,0.5,1],
# ~ vmin=0, vmax=1)
fig.colorbar(cs, cmap='Reds', ax=ax, orientation='horizontal')
plt.title("Density Distribution of all flights in IAGOS dataset")
plt.xticks([-180,-120,-60,0,60,120,180])
plt.yticks([-90,-60,-30,0,30,60,90])
plt.xlabel('Longitude (Degrees)')
plt.ylabel('Latitude (Degrees)')

plt.savefig("IAGOS_Flight_Density_Distribution_0.625deg.png", dpi=300)
print(os.getcwd())
