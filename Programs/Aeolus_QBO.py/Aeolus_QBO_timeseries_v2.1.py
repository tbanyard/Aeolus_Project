#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
---v2.1---Rewriting code to deal with multiprocessing-------------------
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
### deltatime::phdfunctions - Finds difference between two times
day_pop = deltatime(ds.time.values[0], ds.time.values[-1], 'days')
alts = np.linspace(0,30000, 31)
z = np.zeros((len(alts), day_pop))
print(day_pop)

# Create a function which places data variables into a dictionary
### data_dict::phdfunctions
data = {
    'data_lat': ds.data_vars['lat'][:],
    'data_lon': ds.data_vars['lon'][:],
    'data_alt': ds.data_vars['alt'][:],
    'data_u_proj': ds.data_vars['Zonal_wind_projection'][:],
    'data_HLOS_wind': ds.data_vars['Rayleigh_HLOS_wind_speed'][:],
    'data_QC_Flag_Both': ds.data_vars['QC_Flag_Both'][:]
    }

# Applying QCs
data['data_u_proj'] = xr.where(data['data_QC_Flag_Both'] == 0, -99999999, (data['data_u_proj']))

# Below I cap u_proj using xr.where
data['data_u_proj_capped'] = xr.where(data['data_u_proj']>250000, -99999999, (data['data_u_proj']))

# Create a binary array which restricts to an equatorial band between 10N-10S 
lat_band = xr.where(data['data_lat']<-10, 0, (xr.where(data['data_lat']>10, 0, 1)))

# Now restrict all arrays to be None (or -99999999) outside of the latitude band
data['data_u_proj_capped_band'] = xr.where(lat_band==0, -99999999, data['data_u_proj_capped'])
data['data_lat_band'] = xr.where(lat_band==0, -99999999, data['data_lat'])
data['data_lon_band'] = xr.where(lat_band==0, -99999999, data['data_lon'])
data['data_alt_band'] = xr.where(lat_band==0, -99999999, data['data_alt'])

"""# Now I use groupby to group by the day, test2 simply gives the mean for all alts.
test = data['data_u_proj_capped'].groupby("time.day")
test2 = data['data_u_proj_capped'].groupby("time.day").mean()
print(test2.values)
print("==/1/==")"""

"""# returns an array pointing to the bin that each value is in
def altgroup(x):
    binplace = np.digitize(x, alts)
    return binplace

# groupby which applies the altgroup
test3 = data['data_alt'].groupby("time.day").apply(altgroup)
print(test3)
print("==/2/==")"""

# Grouping and Binning data
test3b = data['data_alt_band'].groupby("time.day")
# .group_bins ***LIKELY DEPRECATED***
# test3c = data['data_alt'].groupby_bins("time.day", bins=alts/1000)
# print("test3c: ", test3c.groups)

# Binning the altitudes according to the array alts
binned_alts = xr.apply_ufunc(da.digitize, data['data_alt_band'], alts, dask='allowed', output_dtypes = [float])

u_proj = data['data_u_proj_capped_band'].groupby("time.day")
print(u_proj._unique_coord)

hered = binned_alts.groupby("time.day")

print("u_proj: ", u_proj)
print("binned_alts: ", binned_alts)
print("hered: ", hered)
print("==/3/==")

# @delayed(pure=True)
def mean_each_alt(j):
    new_winds = xr.where(binned_alts == j, data['data_u_proj_capped_band'], None)
    new_winds = xr.where(new_winds == -99999999, None, new_winds)
    new_winds2 = new_winds.groupby("time.day")
    return new_winds2.mean()
    
# @delayed(pure=True)
# def loop_it(z):
    # for j in range(31):
        # print("z: ", z[0])
        # print("mean: ", mean_each_alt(j))
        # z[j][:] = mean_each_alt(j)

# @delayed(pure=True)
def buildz(z, j, num):
    z[j][:] = num
    return z

for j in range(31):
    p = multiprocessing.Process(target = mean_each_alt, args = (j,))
    processes.append(p)
    p.start()
    first = last
    last += segment_size

# Testing
def runbuildz(alt_segment, procfile):
    z = np.zeros(len(alt_segment))
    for j in alt_segment:
        num = mean_each_alt(j)
        z = buildz(z, j, num)
    
# Run as multiple processes for speed
processes = [] # Initialise list for processes
ncores = 4 # Number of CPU cores / workers used, 2GB required per core!
segment_size = np.ceil(31 / ncores) # Alts in worker segment, except last
alts_km = alts / 1000 # List of altitudes in km
first, last = 0, segment_size # First and last elements of first segment

# Multiple processes loop
for n in range(ncores):
    procfile = "proc" + str(n) + ".npy"
    alt_segment = alts_km[first:last]
    p = multiprocessing.Process(target = runbuildz, args = (alt_segment, procfile))
    processes.append(p)
    p.start()
    first = last
    last += segment_size
    if last > 31:
        last == 31
    
for j in range(31):
    num = mean_each_alt(j)
    z = buildz(z, j, num)
    print(j)

time.sleep(3)

print(z)

with ProgressBar():
    z.compute()

time.sleep(30)

# df = z.persist()

# time.sleep(10)

# df = client.persist(z)

z = visualize(filename='daskgraph.png')
# z.compute()

print(z)

# Saving timeseries to npy file
np.save('u_proj_timeseries2.npy', z)

# Time since start of program
timestamp(pstartTime)

for j in range(31):
    print("a")
    new_winds = xr.where(binned_alts == j, data['data_u_proj_capped_band'], None)
    print("b")
    new_winds = xr.where(new_winds == -99999999, None, new_winds)
    print("c")
    new_winds2 = new_winds.groupby("time.day")
    print("d")
    z[j][:] = new_winds2.mean()
    print("j = ", j, " Done.")

# Saving timeseries to npy file
np.save('u_proj_timeseries.npy', z)

print("==/4/==")
timestamp(pstartTime)
days = np.linspace(1, 30, 30)
x, y = days, alts
fig = plt.figure()
cs = plt.contourf(x, y/1000, z, cmap='RdBu_r')
fig.colorbar(cs, cmap='RdBu_r', orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms-1')
plt.savefig("xarraytest3.png", dpi=300)


# @delayed(pure=True, nout=31)
def digitize(i):
    return xr.apply_ufunc(da.digitize, i, alts, dask='allowed', output_dtypes = [float])

# @delayed(pure=True, nout=31)
def loop(altplaces, u_proj):
    print(len(altplaces))
    for i in altplaces:
        hereb = digitize(i)
        print("len(i): ", len(i))
        print(i)
        for j in range(31):
            single_alt = xr.where(hereb == j, u_proj, 0)
            

answeris = loop(test3b, u_proj)
# answeris.compute()
    

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
