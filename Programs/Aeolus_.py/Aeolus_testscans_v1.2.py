#!/usr/bin/env python3
"""
Aeolus scans testing file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Testing-------------------------------------------------------
---v1.2---N.B. v1.0 and v1.1 array indices will no longer work----------
----------Plots a cross-section of HLOS winds for a single orbit--------
----------[DEPRECATED]-There_is_a_newer_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .DBL files downloaded from the Aeolus database and produces a
graphical plot of the satellite pass.
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import os
os.putenv('CODA_DEFINITION',
'/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import datetime
# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs')
from phdfunctions import timeseriesplot, find_nearest
from functions import load_hdr_tags, load_dbl_tags
from scipy.interpolate import griddata

dbl = '/home/tpb38/PhD/Bath/Aeolus/DATA/'
hdr = dbl # Storing hdr and dbl files in same directory for now
dbl += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.DBL'
hdr += 'AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.HDR'

# Opening Files
pf = coda.open(dbl)
pfhdr = coda.open(hdr)

# Load tags if required
load_hdr = False
load_dbl = False

if load_hdr == True:
	Earth_Explorer_Header, Fixed_Header, File_Name, \
	File_Description, Notes, Mission, File_Class, File_Type, \
	File_Version, Validity_Period, Validity_Start, Validity_Stop, \
	Source, System, Creator, Creator_Version, Creation_Date, \
	Variable_Header, Main_Product_Header, Product, Proc_Stage, \
	Ref_Doc, Acquisition_Station, Proc_Center, Proc_Time, \
	Software_Ver, Baseline, Sensing_Start, Sensing_Stop, Phase, \
	Cycle, Rel_Orbit, Abs_Orbit, State_Vector_Time, Delta_UT1, \
	X_Position, Y_Position, Z_Position, X_Velocity, Y_Velocity, \
	Z_Velocity, Vector_Source, Utc_Sbt_Time, Sat_Binary_Time, \
	Clock_Step, Leap_Utc, Gps_Utc_Time_Difference, Leap_Sign, \
	Leap_Err, Product_Err, Tot_Size, Sph_Size, Num_Dsd, Dsd_Size, \
	Num_Data_Sets, Specific_Product_Header, Sph_Descriptor, \
	NumMeasurements, NumMieGroups, NumRayleighGroups, \
	NumMieWindResults, NumRayleighWindResults, NumMieProfiles, \
	NumRayleighProfiles, NumAMDprofiles, Intersect_Start_Lat, \
	Intersect_Start_Long, Intersect_Stop_Lat, Intersect_Stop_Long, \
	Sat_Track, List_of_Valid_Mie_Profile_Counts, \
	Valid_Mie_Profile_Count, List_of_Valid_Rayleigh_Profile_Counts, \
	Valid_Rayleigh_Profile_Count, List_of_Invalid_Mie_Profile_Counts, \
	Invalid_Mie_Profile_Count, List_of_Invalid_Rayleigh_Profile_Counts,\
	Invalid_Rayleigh_Profile_Count, Num_Profiles_Surface_Mie, \
	Num_Profiles_Surface_Ray, List_of_Valid_L2B_Mie_Wind_Counts, \
	Valid_L2B_Mie_Wind_Count, List_of_Valid_L2B_Rayleigh_Wind_Counts, \
	Valid_L2B_Rayleigh_Wind_Count, List_of_Invalid_L2B_Mie_Wind_Counts,\
	Invalid_L2B_Mie_Wind_Count, \
	List_of_Invalid_L2B_Rayleigh_Wind_Counts, \
	Invalid_L2B_Rayleigh_Wind_Count, List_of_Dsds, Dsd, Meas_Map_ADS, \
	Mie_Grouping_ADS, Rayleigh_Grouping_Map, Mie_Geolocation_ADS, \
	Rayleigh_Geolocation_ADS, AMD_Product_Confid_Data_ADS, \
	Meas_Product_Confid_Data_ADS, Mie_Wind_Product_Conf_Data_ADS, \
	Rayl_Wind_Prod_Conf_Data_ADS, Mie_Wind_MDS, Rayleigh_Wind_MDS, \
	Mie_Profile_MDS, Rayleigh_Profile_MDS, Aeolus_Level_1B_Product, \
	Aux_Met_Product, Aeolus_RBC, Clim_Product, Cal_Product, \
	Level_2B_Proc_Params = load_hdr_tags(hdr)
	
if load_dbl == True:
	Product, Proc_Stage, Ref_Doc, Acquisition_Station, \
	Proc_Center, Proc_Time, Software_Ver, Baseline, Sensing_Start, \
	Sensing_Stop, Phase, Cycle, Rel_Orbit, Abs_Orbit, \
	State_Vector_Time, Delta_UT1, X_Position, Y_Position, Z_Position, \
	X_Velocity, Y_Velocity, Z_Velocity, Vector_Source, Utc_Sbt_Time, \
	Sat_Binary_Time, Clock_Step, Leap_Utc, Gps_Utc_Time_Difference, \
	Leap_Sign, Leap_Err, Product_Err, Tot_Size, Sph_Size, Num_Dsd, \
	Dsd_Size, Num_Data_Sets, Sph_Descriptor, NumMeasurements, \
	NumMieGroups, NumRayleighGroups, NumMieWindResults, \
	NumRayleighWindResults, NumMieProfiles,	NumRayleighProfiles, \
	NumAMDprofiles, Intersect_Start_Lat, Intersect_Start_Long, \
	Intersect_Stop_Lat, Intersect_Stop_Long, Sat_Track, \
	Valid_Mie_Profile_Count, Valid_Rayleigh_Profile_Count, \
	Invalid_Mie_Profile_Count, Invalid_Rayleigh_Profile_Count, \
	Num_Profiles_Surface_Mie, Num_Profiles_Surface_Ray, \
	Valid_L2B_Mie_Wind_Count, Valid_L2B_Rayleigh_Wind_Count, \
	Invalid_L2B_Mie_Wind_Count, Invalid_L2B_Rayleigh_Wind_Count, \
	Dsd, Meas_Map, Mie_Map_of_L1B_Meas_Used, \
	Rayleigh_Map_of_L1B_Meas_Used, Mie_Grouping, Rayleigh_Grouping, \
	Mie_Geolocation, Rayleigh_Geolocation, AMD_Product_Confid_Data, \
	Meas_Product_Confid_Data, Mie_Wind_Prod_Conf_Data, \
	Rayleigh_Wind_Prod_Conf_Data, Mie_HLOS_Wind, Rayleigh_HLOS_Wind, \
	Mie_Profile, Rayleigh_Profile = load_dbl_tags(dbl)

# Fetching Data
rayleigh_wind_velocity = coda.fetch(pf, 'rayleigh_hloswind', -1,
	'windresult/rayleigh_wind_velocity')
rayleigh_latitude = coda.fetch(pf, 'rayleigh_geolocation', -1,
	'windresult_geolocation/latitude_cog')
rayleigh_longitude = coda.fetch(pf, 'rayleigh_geolocation', -1,
	'windresult_geolocation/longitude_cog')
rayleigh_altitude = coda.fetch(pf, 'rayleigh_geolocation', -1,
	'windresult_geolocation/altitude_vcog')
rayleigh_date_time = coda.fetch(pf, 'rayleigh_geolocation', -1,
	'windresult_geolocation/datetime_cog')
mie_wind_velocity = coda.fetch(pf, 'mie_hloswind', -1,
	'windresult/mie_wind_velocity')
mie_latitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/latitude_cog')
mie_longitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/longitude_cog')
mie_altitude = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/altitude_vcog')
mie_date_time = coda.fetch(pf, 'mie_geolocation', -1,
	'windresult_geolocation/datetime_cog')

# Mie and Rayleigh Grouping
Mie_Grouping = coda.fetch(pf,
	'mie_grouping')
Rayleigh_Grouping = coda.fetch(pf,
	'rayleigh_grouping')

# Validity Flags
Mie_Wind_Prod_Conf_Data = coda.fetch(pf,
	'mie_wind_prod_conf_data')
Rayleigh_Wind_Prod_Conf_Data = coda.fetch(pf,
	'rayleigh_wind_prod_conf_data')
Mie_HLOS_Wind = coda.fetch(pf,
	'mie_hloswind')
Rayleigh_HLOS_Wind = coda.fetch(pf,
	'rayleigh_hloswind')

# Mie Wind Product Confidence Data
itrn1 = 0
itrn2 = 0
for i in range(len(mie_wind_velocity)):
	flags1 = Mie_Wind_Prod_Conf_Data[i][2][1]
	flags2 = Mie_Wind_Prod_Conf_Data[i][2][2]
	flags3 = Mie_Wind_Prod_Conf_Data[i][2][3]
	flags4 = Mie_Wind_Prod_Conf_Data[i][2][4]
	tot = flags1 + flags2 + flags3 + flags4
	if tot != 0:
		itrn1 += 1
	else:
		itrn2 += 1

print("Mie Wind Product Conf (Invalid) = ", itrn1)
print("Mie Wind Product Conf (Valid) = ", itrn2)

# Rayleigh Wind Product Confidence Data
itrn3 = 0
itrn4 = 0
for i in range(len(rayleigh_wind_velocity)):
	flags1 = Rayleigh_Wind_Prod_Conf_Data[i][2][1]
	flags2 = Rayleigh_Wind_Prod_Conf_Data[i][2][2]
	flags3 = Rayleigh_Wind_Prod_Conf_Data[i][2][3]
	flags4 = Rayleigh_Wind_Prod_Conf_Data[i][2][4]
	tot = flags1 + flags2 + flags3 + flags4
	if tot != 0:
		itrn3 += 1
	else:
		itrn4 += 1

print("Rayleigh Wind Product Conf (Invalid) = ", itrn3)
print("Rayleigh Wind Product Conf (Valid) = ", itrn4)

# Mie HLOS Wind Validity
itrn5 = 0
itrn6 = 0
for i in range(len(Mie_HLOS_Wind)):
	validity_flag = Mie_HLOS_Wind[i][2][2]
	if validity_flag == 0:
		itrn5 += 1
	elif validity_flag == 1:
		itrn6 += 1

print("Mie HLOS Wind Validity (Invalid) = ", itrn5)
print("Mie HLOS Wind Validity (Valid) = ", itrn6)

# Rayleigh HLOS Wind Validity
itrn7 = 0
itrn8 = 0
for i in range(len(Rayleigh_HLOS_Wind)):
	validity_flag = Rayleigh_HLOS_Wind[i][2][2]
	if validity_flag == 0:
		itrn7 += 1
	elif validity_flag == 1:
		itrn8 += 1

print("Rayleigh HLOS Wind Validity (Invalid) = ", itrn7)
print("Rayleigh HLOS Wind Validity (Valid) = ", itrn8)

# Mie HLOS Observation Type
itrn9 = 0
itrn10 = 0
for i in range(len(Mie_HLOS_Wind)):
	validity_flag = Mie_HLOS_Wind[i][2][1]
	if validity_flag == 1:
		itrn9 += 1
	elif validity_flag == 2:
		itrn10 += 1

print("Mie HLOS Observation Type (Cloudy) = ", itrn9)
print("Mie HLOS Observation Type (Clear) = ", itrn10)

# Rayleigh HLOS Observation Type
itrn11 = 0
itrn12 = 0
for i in range(len(Rayleigh_HLOS_Wind)):
	validity_flag = Rayleigh_HLOS_Wind[i][2][1]
	if validity_flag == 1:
		itrn11 += 1
	elif validity_flag == 2:
		itrn12 += 1

print("Rayleigh HLOS Observation Type (Cloudy) = ", itrn11)
print("Rayleigh HLOS Observation Type (Clear) = ", itrn12)

# Initialise arrays
mie_times = []
mie_alts = []
mie_lats = []
mie_lons = []
mie_wvs = []
rayleigh_times = []
rayleigh_alts = []
rayleigh_lats = []
rayleigh_lons = []
rayleigh_wvs = []

# Build arrays
for j in range(len(mie_date_time)):
	observation_type = Mie_HLOS_Wind[j][2][1]
	validity_flag = Mie_HLOS_Wind[j][2][2]
	if observation_type == 1: # i.e. cloudy
		if validity_flag == 1:
			mie_times.append(mie_date_time[j])
			mie_alts.append(mie_altitude[j])
			mie_lats.append(mie_latitude[j])
			mie_lons.append(mie_longitude[j])
			mie_wvs.append(mie_wind_velocity[j])
			
for k in range(len(rayleigh_date_time)):
	observation_type = Rayleigh_HLOS_Wind[k][2][1]
	validity_flag = Rayleigh_HLOS_Wind[k][2][2]
	if observation_type == 2: # i.e. clear
		if validity_flag == 1:
			rayleigh_times.append(rayleigh_date_time[k])
			rayleigh_alts.append(rayleigh_altitude[k])
			rayleigh_lats.append(rayleigh_latitude[k])
			rayleigh_lons.append(rayleigh_longitude[k])
			rayleigh_wvs.append(rayleigh_wind_velocity[k])
			
# Plotting data
os.chdir('..')
os.chdir('..')
os.chdir('Plots')

print("Mie Times has length: ", len(mie_times))
print("Rayleigh Times has length: ", len(rayleigh_times))
print("Rayleigh Grouping has length: ", len(Rayleigh_Grouping))

for t in range(len(rayleigh_times)): # testing
	print(rayleigh_times[t])

# ~ RG_elmnt = 0
# ~ for t in range(len(rayleigh_times)):
	# ~ if rayleigh_times[t] < Rayleigh_Grouping[RG_elmnt][1]:
		# ~ print(rayleigh_times[t])
	# ~ else:
		# ~ RG_elmnt += 1

# Convert list of Rayleigh Group start times into a sensible format
RG = np.zeros(len(Rayleigh_Grouping))
for g in range(len(RG)):
	RG[g] = (Rayleigh_Grouping[g][1])

# Initialise meshgrids for x, y and z
alts = np.linspace(0,20000, 21)
x, y = np.meshgrid(RG, alts)
print(x)
print(y)
# ~ z = [[0 for _ in range(len(RG))] for _ in range(len(alts))] # Lists
# ~ z_itrn = [[0 for _ in range(len(RG))] for _ in range(len(alts))]
z = np.zeros((len(alts),len(RG))) # NumPy Arrays
z_itrn = np.zeros((len(alts),len(RG)))
print(np.shape(z))

# Placing wind values into bins of height 1km and width 1 rayleigh group
lastgroupstarttime = 0
for RG_elmnt in range(len(RG)):
	for t in range(len(rayleigh_times)):
		# Find all elements inside this sandwich and add to z & z_itrn:
		if rayleigh_times[t] < Rayleigh_Grouping[RG_elmnt][1] and \
		rayleigh_times[t] >= lastgroupstarttime:
			# Find the nearest altitude level
			val = find_nearest(alts, rayleigh_alts[t])
			alt_elmnt = np.where(alts == val)[0][0]
			# Cap wind speeds to 100 m/s
			if np.abs(rayleigh_wvs[t]) < 10000:
				z[alt_elmnt][RG_elmnt] += rayleigh_wvs[t]
				z_itrn[alt_elmnt][RG_elmnt] += 1
	lastgroupstarttime = Rayleigh_Grouping[RG_elmnt][1]

# Find the mean for each bin
z /= 100 * z_itrn # Factor of 100 for conversion from cm/s to m/s
print(z)

# Plotting
fig = plt.figure()
ax1 = fig.add_subplot(111)
cs = plt.contourf(x,y,z, cmap='RdBu')
ax2 = ax1.twinx()
ax2.plot(rayleigh_times, rayleigh_lats, c='black', marker='.',
	label='latitude', linewidth=0.1)
ax1.set_xlabel('Time')
ax1.set_ylabel('Altitude / m')
ax2.set_ylabel('Latitude / $^\circ$')
plt.title('Aeolus Orbit HLOS Wind Cross-section')
fig.colorbar(cs, cmap='RdBu', ax=ax1, orientation='horizontal',
	label='HLOS Wind Speed / ms-1')
plt.savefig('test2.png',dpi=300)

"""
X1 = mie_times[860:910]
X2 = rayleigh_times
Y1 = mie_alts[860:910]
Y2 = mie_lats[860:910]
Y3 = mie_lons[860:910]
Y4 = mie_wvs[860:910]
Y5 = rayleigh_alts
Y6 = rayleigh_lats
Y7 = rayleigh_lons
Y8 = rayleigh_wvs
variable = 'Longitude'
variable2 = 'Latitude'
# ~ timeseriesplot(X1, Y3,
	plottitle = \
	'Variation in Aeolus\' longitude and latitude \n during one full orbit for Mie (QC)',
	date_form = '%H:%M', minor_date_form = '%M:%S', data_type = 'coda',
	size = 0.5, color = 'black', marker = '+', variable = variable,
	variable2 = variable2, legend = 0, l_adj = 0.15, r_adj=0.85)
# ~ plt.savefig("timeseries.png", dpi=300)
# ~ print(coda.time_to_string(date_time[21000]))

xi = np.linspace(rayleigh_times[0], rayleigh_times[-1], 100)
yi = np.linspace(rayleigh_alts[0], rayleigh_alts[-1], 10)
xi, yi = np.meshgrid(xi, yi)
z = rayleigh_wvs
zi = griddata((rayleigh_times, rayleigh_alts),z,(xi,yi),method='linear')
# ~ zi[mask] = np.nan
date_time = coda.time_to_utcstring(mie_times[:])
date_time = np.array([datetime.strptime(date,
	'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
fig = plt.figure()
ax = fig.add_subplot(111)
cs = plt.contourf(date_time,yi,zi, cmap='RdBu')
# ~ plt.plot(rayleigh_times,rayleigh_alts,'k.')
plt.xlabel('xi',fontsize=16)
plt.ylabel('yi',fontsize=16)
fig.colorbar(cs, cmap='RdBu', ax=ax, orientation='vertical')
plt.savefig('interpolated.png',dpi=300)
"""

"""
# Initialise arrays
alts = []
lats = []
lons = []
mwvs = []

# Separating into a single level
for t in range(len(X)):
	if Y1[t]>9900 and Y1[t]<10000:
		alts.append(Y1[t])
		lats.append(Y2[t])
		lons.append(Y3[t])
		mwvs.append(Y4[t])

date_time = coda.time_to_utcstring(X[:])
date_time = np.array([datetime.strptime(date,
	'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
X = date_time

xx,yy = np.meshgrid(X, Y1)
# ~ print(xx)

"""	
"""
xi = np.linspace(date_time[0], date_time[-1], 1000)
yi = np.linspace(0, 20000, 1000)
Z = Y4

X=X/100000

print(X)

# ~ date_time = coda.time_to_utcstring(X[:])
# ~ date_time = np.array([datetime.strptime(date,
	'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
# ~ X = date_time

# ~ xi = np.linspace(28.5, 26.5, 1000)
# ~ yi = np.linspace(341.6, 341.1, 1000)

zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.contour(xi, yi, zi, linewidths = 0.1, colors = 'k')
cs = ax.contourf(xi, yi, zi, cmap='jet')
# ~ ax.contourf(x, y, z, cmap='RdBu', shading='interp')
plt.savefig("test1.png", dpi=300)

# ~ ax.scatter(Y3, Y2, s=5, c='blue', marker='.')
# ~ ax.set_xlim(340, 341)
# ~ ax.set_ylim(20, 26)
# ~ plt.savefig("latlon.png", dpi=300)"""


