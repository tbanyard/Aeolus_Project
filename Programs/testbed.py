#!/usr/bin/env python3
"""
General testbed for plotting data from saved .npy files
========================================================================
------------------------------------------------------------------------
---v1.0---[1st PhD CODE REVAMP]-----------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .npy files created from programs and plots them to suit the
immediate requirement
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

# .npy file to load for plotting
os.chdir('Programs')
print(os.getcwd())
npyfile = 'u_proj_timeseries4.npy'

# Load data for plotting
z = np.load(npyfile, allow_pickle=True)
x = np.linspace(1, 31, 31)
y = np.linspace(0, 30000, 31)

# Change to plot directory to save plots
os.chdir('..')
os.chdir('Plots')

# Plotting data
fig = plt.figure()
ax = fig.add_subplot(111)
cs = plt.contourf(x, y/1000, z, levels=[-1000,-750,-500,-250,0,250,500,750,1000], vmin= -1000, vmax = 1000, cmap='RdBu_r')
fig.colorbar(cs, cmap='RdBu_r', orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms-1')
plt.savefig("xarraytest6.png", dpi=300)

"""
# Load data for plotting
data = np.load('2000_flight_distributions.npy', allow_pickle=True)
data_lat = np.load("lat_grid.npy", allow_pickle=True)
data_lon = np.load("lon_grid.npy", allow_pickle=True)

np.save('data.npy', data)

printfullarrays()

# Change to plot directory to save plots
os.chdir('Plots')

# Plotting data
fig = plt.figure()
ax = fig.add_subplot(111)
# map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            # llcrnrlon=-180,urcrnrlon=180,resolution='l', ax=ax)
# ~ map = Basemap(projection='ortho',lat_0=45,lon_0=0,\
            # ~ llcrnrlon=-180,urcrnrlon=180,resolution='i', ax=ax)
# map.drawcoastlines(linewidth=0.5)
# draw parallels and meridians.
# map.drawparallels(np.arange(-90.,91.,30.),linewidth=0.5)
# map.drawmeridians(np.arange(-180.,180.,60.),linewidth=0.5)
lons, lats = np.meshgrid(data_lon, data_lat)
x,y = lons, lats
cs = plt.contourf(x,y,np.log(data), cmap='Reds',
levels=[0,1,2,3,4,5,6,7,8,9,10], vmin=0, vmax=10, extend="both")
# ~ cs = map.contourf(x,y,data==0, cmap='Greys', levels=[0,0.5,1],
# ~ vmin=0, vmax=1)
fig.colorbar(cs, cmap='Reds', ax=ax, orientation='horizontal')
plt.title("Density Distribution of all flights in IAGOS dataset")
plt.xticks([-180,-120,-60,0,60,120,180])
plt.yticks([-90,-60,-30,0,30,60,90])
plt.xlabel('Longitude (Degrees)')
plt.ylabel('Latitude (Degrees)')

plt.savefig("IAGOS_Flight_Density_Distribution_0.625deg.png", dpi=300)"""

