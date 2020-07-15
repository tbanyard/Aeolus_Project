#!/usr/bin/env python3
"""
Aeolus data load from netCDF format for QBO test
========================================================================
------------------------------------------------------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
---v2.1---Rewriting code to deal with multiprocessing-------------------
---v2.2---Tidying code, fixing groupby day issue------------------------
---v2.3---Removing print statements...----------------------------------
---v2.4---Tidied version of timeseries code using xarray and dask-------
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

print(dirs['Programs']) # Generalise as much as possible. Multiple .nc files?
print(dirs['Plots']) # Note, when Plots is created, I can use the exceptions.
print(dirs['ncdir'])

# Find and read netCDF data
print("Opening Data from directory ", ncdir, " lazily...")
ds = xr.open_mfdataset(dirs['ncdir'] + '/*.nc',
    combine = 'nested', concat_dim = 'time')
    
# Create plot array (z) with a size matching start and end dates of data
### deltatime::phdfunctions - Finds difference between two times
day_pop = deltatime(ds.time.values[0], ds.time.values[-1], 'days')
alts = np.linspace(0,30000, 31)
z = np.zeros((len(alts), day_pop+1))
print("Data loaded and plot array generated.")

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

# Binning the altitudes according to the array alts
binned_alts = xr.apply_ufunc(da.digitize, data['data_alt_band'], alts, dask='allowed', output_dtypes = [float])

print("==/1/==")

# Function to take the mean for each altitude
def mean_each_alt(j):
    new_winds = xr.where(binned_alts == j, data['data_u_proj_capped_band'], None)
    new_winds = xr.where(new_winds == -99999999, None, new_winds)
    new_winds2 = new_winds.sortby('time').resample(time='1D')
    answer = new_winds2.mean().values
    return answer

# Function to build z using the answer from mean_each_alt
def runbuildz(alt_segment, day_pop, procfile):
    z = np.zeros((len(alt_segment), day_pop+1))
    itrn = 0
    for j in alt_segment:
        j = int(j)
        z[itrn][:] = mean_each_alt(j)
        itrn += 1
    np.save(procfile, z, allow_pickle=True)
    
# Run as multiple processes for speed
processes = [] # Initialise list for processes
ncores = 1 # Number of CPU cores / workers used, 2GB required per core!
segment_size = np.ceil(31 / ncores) # Alts in worker segment, except last
alts_km = alts / 1000 # List of altitudes in km
first, last = 0, segment_size # First and last elements of first segment

# Execute runbuildz, and hence mean_each_alt on multiple cores
with ProgressBar(minimum=3, dt=0.2):
    if __name__ == '__main__':
        # Enables function to be executed using multiple processes/cores
        for n in range(ncores):
            first, last = int(first), int(last)
            procfile = "proc" + str(n) + ".npy"
            alt_segment = alts_km[first:last]
            p = multiprocessing.Process(target = runbuildz,
                args = (alt_segment, day_pop, procfile))
            processes.append(p)
            p.start()
            first = last
            last += segment_size
            if last > 31.0:
                last = 31.0
        
        # Join Processes                        
        for process in processes:
            process.join()
        
        # Concatenate data from temporary .npy files
        j0, jn = 0, int(segment_size)
        for n in range(ncores):
            j0, jn = int(j0), int(jn)
            procfile = "proc" + str(n) + ".npy"
            z[j0:jn][:] = np.load(procfile)
            removefile(procfile)
            j0 = jn
            jn += segment_size
            if jn > 31:
                jn == 31

print(timestamp(pstartTime))

# Saving timeseries to npy file
np.save('u_proj_timeseries6.npy', z)

# Time taken for the program
pduration = datetime.now() - pstartTime
print('That program took ', pduration, ' seconds')
