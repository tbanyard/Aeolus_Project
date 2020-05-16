#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
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
ncdir = '/home/tpb38/PhD/Bath/Aeolus/NC_FullQC'
ncdir2 = '/home/tim/Documents/Bath/NC_AEOLUS_DATA_APR'

"""=================================================================="""
"""===============Reading Data and Array Initialisation=============="""
"""=================================================================="""
# Generate dictionary to store directories
dirs = dict([('Programs', os.getcwd())])
dirs['Plots'] = dirs['Programs'][:-8] + 'Plots'
dirs['ncdir'], dirs['ncdir2'] = ncdir, ncdir2

print(dirs['Programs']) # Generalise as much as possible!! Multiple .nc files???
print(dirs['Plots']) # Note, when Plots is created, I can use the exceptions.
print(dirs['ncdir'])

@delayed(pure=True)
def add(a, b):
    for t in range(20000000):
        y = t + 4
    return a + b
    
@delayed(pure=True)
def inc(a):
    for t in range(20000000):
        y = t + 4
    return a + 1
    
x1 = inc(2)
y1 = inc(1)
z1 = add(x1, y1)
x2 = inc(3)
y2 = inc(4)
z2 = add(x2, y2)
a = add(z1, z2)
a.visualize(filename='daskgraph.png')
one = a + a
two = a * a
three = a + a
four = a * a
five = one + two + three + four
five.compute()
with ProgressBar():
    a.compute()

with ProgressBar():
    print(x1.compute())

# time.sleep(1000)

# Find and read netCDF data
ds = xr.open_mfdataset(dirs['Programs'] + '/timdatamarQCdDesc2.nc', 
    combine = 'nested', concat_dim = 'time')
ds2 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')

print(ds2.nbytes/1e6)
print(pstartTime)
print(ds.data_vars)

ds.data_vars['Zonal_wind_projection'][:,1:].plot()
plt.savefig('test0123.png', dpi=600)
print(os.getcwd())

ds3 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds4 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds5 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds6 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds7 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds8 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds9 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
ds10 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    


# tpbcmaps
