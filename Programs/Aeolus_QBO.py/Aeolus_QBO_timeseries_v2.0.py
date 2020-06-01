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

# Find and read netCDF data
ds = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
# Create plot array (z) with a size matching start and end dates of data
### deltatime::phdfunctions
day_pop = deltatime(ds.time.values[0], ds.time.values[-1], 'days')
alts = np.linspace(0,30000, 31)
z = np.zeros((len(alts), day_pop+1))

# Create a function which places data variables into a dictionary
### data_dict::phdfunctions
data = {
    'data_lat': ds.data_vars['lat'][:],
    'data_lon': ds.data_vars['lon'][:],
    'data_alt': ds.data_vars['alt'][:],
    'data_u_proj': ds.data_vars['Zonal_wind_projection'][:],
    }

# data_lat = ds.data_vars['lat'][:]
# data_lon = ds.data_vars['lon'][:]
lat_band = xr.where(data['data_lat']<-5, 0, (xr.where(data['data_lat']>5, 0, 1)))

####### READ THIS TOMORROW #######
# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, 0, (data['data_u_proj']))

# Now I use groupby to group by the day, test2 simply gives the mean for all alts.
test = data['data_u_proj_capped'].groupby("time.day")
test2 = data['data_u_proj_capped'].groupby("time.day").mean()
print(test2.values)
print("==/1/==")

# returns an array pointing to the bin that each value is in
def altgroup(x):
    binplace = np.digitize(x, alts)
    return binplace

# groupby which applies the altgroup
test3 = data['data_alt'].groupby("time.day").apply(altgroup)
print(test3)
print("==/2/==")

# Where I need to start next!
test3b = data['data_alt'].groupby("time.day")
u_proj = data['data_u_proj_capped'].groupby("time.day")
herea = xr.apply_ufunc(da.digitize, test3b, alts, dask='allowed', output_dtypes = [float])
print("herea: ", np.mean(herea.values))
print("==/3/==")

@delayed(pure=True, nout=31)
def digitize(i):
    return xr.apply_ufunc(da.digitize, i, alts, dask='allowed', output_dtypes = [float])

@delayed(pure=True, nout=31)
def loop(altplaces, u_proj):
    for i in altplaces:
        hereb = digitize(i)
        print("len(i): ", len(i))
        for j in range(31):
            single_alt = xr.where(hereb == j, u_proj, 0)
            

answeris = loop(test3b, u_proj)
answeris.compute()
    

print("==/4/==")


time.sleep(100000)

def wrap_digitize(data):
    return np.digitize(data, alts)
    
here = xr.apply_ufunc(wrap_digitize, test3b, dask='parallelized', output_dtypes = [float])
print(here)
print("==//==")


test4 = data['data_u_proj_capped'].groupby("time.day").apply(altgroup)

@delayed
def eachday(i):
    x = np.asarray(i)
    print(x)

for i in test:
    intermediate = eachday(i)
    ans = intermediate.compute

# alts = data['data_u_proj_capped'].isel(alt=i).groupby(

print(data['data_u_proj_capped'])

time.sleep(100000)


# Find and read netCDF data
ds = xr.open_mfdataset(dirs['Programs'] + '/timdatamarQCdDesc2.nc', 
    combine = 'nested', concat_dim = 'time')
ds2 = xr.open_mfdataset(dirs['ncdir2'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')

print(ds2.nbytes/1e6)
print(ds2.data_vars['lat'].nbytes)
print(pstartTime)
print(ds2.data_vars)

print('------')

data_lat = ds2.data_vars['lat'][:]
print(data_lat)
lat_band = xr.where(data_lat<-5, 0, (xr.where(data_lat>5, 0, 1)))
print(lat_band)



time.sleep(100000)

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

# Time taken for the program
pduration = datetime.now() - pstartTime
print('That program took ', pduration, ' seconds')
