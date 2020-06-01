#!/usr/bin/env python3
"""
PhD Functions file
================================================================================
--------------------------------------------------------------------------------
----------------------------No version control----------------------------------
--------------------------------------------------------------------------------
================================================================================
Provides functions for all files for my PhD
================================================================================
"""

"""
Note: This file has an 80 column width to cater for the size of some of the 
programs.
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import os
import errno
from datetime import timedelta, datetime
from scipy.interpolate import griddata
import calendar
import time
import coda
import warnings

"""=========================================================================="""
"""------------------Functions for dealing with datasets---------------------"""
"""=========================================================================="""

def round_down_3_hours(t):
	"""Rounds down to nearest 3 hours"""
	return t.replace(second=0, microsecond=0, minute=0, \
	hour=(t.hour//3)*3)
    
def round_up_3_hours(t):
	"""Rounds up to nearest 3 hours"""
	if t.hour >= 21:
		time = t.replace(second=0, microsecond=0, minute=0, \
		hour=((t.hour+3)//3)*3-24)
		YYYY, MM, DD = str(time.year), str(time.month), str(time.day)
		time = yyyymmdd_nextday(YYYY, MM, DD)
	elif t.hour < 21:
		time = t.replace(second=0, microsecond=0, minute=0, \
		hour=((t.hour+3)//3)*3)
	return time
	
def find_nearest(array, value):
	"""Finds the nearest element in an array to a specified value."""
	if isinstance(value, datetime) == True:
		array = np.asarray(array)
		array = dates.date2num(array)
		value = dates.date2num(value)
		idx = (np.abs(array - value)).argmin()
	else:
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
	return array[idx]
	
def toTimestamp(d):
	"""Converts datetime object to UNIX timestamp value since 1970"""
	return calendar.timegm(d.timetuple())

def linearinterp2kp(x, x_0, x_1, y_0, y_1):
	"""Simple linear interpolation between two known points.
	x is the position of the unknown value, y, to be found between the
	known positions x_0 and x_1, with values y_0 and y_1
	correspondingly."""
	y = (y_0*(x_1-x)+y_1*(x-x_0))/(x_1-x_0)
	return y
	
def yyyymmdd_to_doy(YYYY, MM, DD):
	"""Converts date in yyyymmdd format to doy format."""
	doy = (datetime(int(YYYY), int(MM),
	int(DD))-datetime(int(YYYY), 1, 1)).days+1
	if len(str(doy)) == 1:
		doy = '00' + str(doy)
	elif len(str(doy)) == 2:
		doy = '0' + str(doy)
	elif len(str(doy)) == 3:
		doy = str(doy)
	return doy
	
def doy_to_yyyymmdd(YYYY, doy):
	"""Converts date in yyyydoy format to yyyymmdd format"""
	dtdate = datetime(YYYY, 1, 1) + timedelta(doy - 1)
	return dtdate
	
def yyyymmdd_nextday(YYYY, MM, DD):
	"""Gives the next day for any date in YYYYMMDD format."""
	date = datetime(int(YYYY), int(MM), int(DD))
	date += timedelta(days=1)
	return date
	
def haversine(lats, lons):
	"""Returns the distance between two lat/lon points"""
	R = 6371 # Earth's radius in km
	# Haversine formula
	a = (np.sin((lats[1]-lats[0])/2) ** 2) + (np.cos(lats[0]) * \
		np.cos(lats[1]) * (np.sin((lons[1]-lons[0])/2) ** 2))
	c = 2 * np.arcsin(np.sqrt(a))
	d = R * c
	return d

class time_type_switcher(object):
    
    def __init__(self):
        self.help = """The class time_type_switcher can be used to run
                    either dtdtfunc or npdtfunc depending on the type of
                    the times input into deltatime. Simply use the
                    function: time_type_switcher.switch"""
    
    def switch(self, first, last, time_type, freq):
        func = time_type.get('func_code') + 'func'
        args = [first, last, freq]
        method = getattr(self, func, lambda: "Invalid switcher input.")(*args)
        return method

    def dtdtfunc(self, first, last, freq):
        """datetime.timedelta function"""
        # N.B. This function is preliminarily deprecated.
        switcher = {
            'weeks': int((last-first).days/7),
            'days': int((last-first).days),
            'hours': int((last-first).seconds/3600),
            'minutes': int((last-first).seconds/60),
            'seconds': int((last-first).seconds),
            'microseconds': int((last-first).microseconds)
            }
        return switcher.get(freq, lambda: "Invalid deltatime frequency.")
         
    def npdtfunc(self, first, last, freq):
        """numpy.timedelta64 function"""
        switcher = {
            # 'years': (last-first).astype('timedelta64[Y]').astype(int)
            # 'months': (last-first).astype('timedelta64[M]').astype(int)
            'weeks': (last-first).astype('timedelta64[W]').astype(int),
            'days': (last-first).astype('timedelta64[D]').astype(int),
            'hours': (last-first).astype('timedelta64[h]').astype(int),
            'minutes': (last-first).astype('timedelta64[m]').astype(int),
            'seconds': (last-first).astype('timedelta64[s]').astype(int),
            'microseconds': (last-first).astype('timedelta64[us]').astype(int)
            }
        
        return switcher.get(freq, lambda: "Invalid deltatime frequency.")
    
def deltatime(first, last, freq = 'days'):
    """Returns the difference between two times"""
   
    # Conditions for condition statement
    # Checking the first time
    conditionsA = {
        "cond1": type(first) == 'numpy.datetime64',
        "cond2": isinstance(first, datetime),
        "cond3": isinstance(first, np.datetime64)
        }
    
    # Checking the last time
    conditionsB = {
        "cond4": type(last) == 'numpy.datetime64',
        "cond5": isinstance(last, datetime),
        "cond6": isinstance(last, np.datetime64)
        }
    
    # Combining as OR statement
    conditions = {
        "condA": not False in conditionsB.values(),
        "condB": not False in conditionsA.values()
        }
    
    # Raise error for invalid type
    if not True in conditions.values():
        raise TypeError('One or both of the input times is of an invalid type')
    
    # Testing that the input times are of the same type
    if type(first) != type(last):
        warnings.warn('Input times are not of the same type')
        
    # Converting time type to np.datetime64 to give accurate timedelta
    first, last = np.datetime64(first), np.datetime64(last)
    
    # Setting the time_type
    try:
        mod = type(last-first).__module__
        name = type(last-first).__name__
        ttype = mod + '.' + name
    except:
        ttype = type(last-first)
    
    # Generating func_code
    if ttype == 'datetime.timedelta':
        fcode = 'dtdt'
    elif ttype == 'numpy.timedelta64':
        fcode = 'npdt'
    
    # Creating time_type dict
    time_type = {
        'type': ttype,
        'func_code': fcode
        }
    
    # Running time_type_switcher to calculate difference with given freq
    p = time_type_switcher()
    answer = p.switch(first, last, time_type, freq)
    
    return answer

"""=========================================================================="""
"""--------------Functions for dealing with specific datasets----------------"""
"""=========================================================================="""
def ERA5_dataload(infile):
	"""This function loads the required data from a netCDF file."""

	# Load data from file
	data = nc.Dataset(infile)

	"""Download variables"""
	# Longitude
	data_lon = data.variables['longitude'][:]
	# Latitude
	data_lat = data.variables['latitude'][:]
	# Level
	data_lev = data.variables['level'][:]
	# Air temperature
	data_temp = data.variables['t'][:]
	# Zonal wind
	data_u = data.variables['u'][:]
	# Meridional wind
	data_v = data.variables['v'][:]

	# Converted time
	data_time = nc.num2date(data.variables['time'][:],\
	calendar = 'standard', units = data.variables['time'].units)

	# Close dataset and save data
	data.close()
	return data_lon, data_lat, data_lev, data_temp, data_u, data_v, \
	data_time
	
"""=========================================================================="""
"""------------------Functions for navigating directories--------------------"""
"""=========================================================================="""

def enterdirectory(strdirectory):
# Enter corresponding directory
	print('\n')
	try:
		os.mkdir(strdirectory)
		print("Directory ", strdirectory, " created")
	except OSError as e:
			if e.errno == errno.EEXIST:
				print("Directory ", strdirectory, " already exists")
			else:
				raise
	os.chdir(strdirectory)
	print(os.getcwd())
	return
	
"""=========================================================================="""
"""-------------------------Functions for Plotting---------------------------"""
"""=========================================================================="""

def timeseriesplot(X, Y, Y2, plottype = 'scatter', plottitle = 'Sample', date_form = '%H:%M',
minor_date_form = '%M', data_type = 'coda', size = 0.5, size2 = 0.5,
color = 'black', color2 = 'blue', marker = '+', marker2 = 'o',
variable = 'variable', variable2 = 'variable2', legend = 0, l_adj = 0.1,
r_adj = 0.9, b_adj = 0.1, t_adj = 0.9, linewidth1 = 0.5, linewidth2 = 0.5):
	# Add in if statements here depending on format of time (X) array.
	if data_type == 'coda':
		print("Data type: coda")
		date_time = coda.time_to_utcstring(X[:])
		date_time = np.array([datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f') for date in date_time])		
		print("date_time length: ", len(date_time))
	else:
		date_time = X
	X = date_time
	print("X: ", X)
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	if plottype == 'scatter':
		ax1.scatter(X, Y, s=size, c=color, marker=marker, label=variable)
		if Y2 != []:
			ax2 = ax1.twinx()
			ax2.scatter(X, Y2, s=size2, c=color2, marker=marker2, label=variable2)
	elif plottype == 'line':
		ax1.plot(X, Y, s=size, c=color, marker=marker, label=variable, linewidth=linewidth1)
		if Y2 != []:
			ax2 = ax1.twinx()
			ax2.plot(X, Y2, s=size2, c=color2, marker=marker2, label=variable2, linewidth=linewidth2)
	if legend == 1:
		plt.legend()
	# ~ ax1.set_title('Sample')
	plt.title(plottitle)
	
	# Setting Date axis
	hours = dates.HourLocator()
	minutes = dates.MinuteLocator()
	date_form = dates.DateFormatter(date_form)
	minor_date_form = dates.DateFormatter(minor_date_form)
	ax1.xaxis.set_major_formatter(date_form)
	ax1.xaxis.set_minor_formatter(minor_date_form)
	ax1.set_xlim(date_time[0], date_time[-1])
	ax1.set_xlabel('Time')
	ax1.set_ylabel(variable)
	if Y2 != []:
		ax2.set_ylabel(variable2)
	# Adjust y axis for big y values
	plt.gcf().subplots_adjust(left=l_adj)
	plt.gcf().subplots_adjust(right=r_adj)
	# Adjust x axis for big x values
	plt.gcf().subplots_adjust(bottom=b_adj)
	plt.gcf().subplots_adjust(top=t_adj)
	# ~ ax1.set_ylim(8920, 8950)
	
	return
	
def fixmydateaxis(ax1, date_time, date_form = '%H:%M', minor_date_form = '%M'):
	# Setting Date axis
	hours = dates.HourLocator()
	minutes = dates.MinuteLocator()
	date_form = dates.DateFormatter(date_form)
	minor_date_form = dates.DateFormatter(minor_date_form)
	ax1.xaxis.set_major_formatter(date_form)
	ax1.xaxis.set_minor_formatter(minor_date_form)
	ax1.set_xlim(date_time[0], date_time[-1])
	ax1.set_xlabel('Time')
	
	return

def fixnanswithmean(z):
	where_are_NaNs = np.isnan(z)
	mean = np.nanmean(z)
	z[where_are_NaNs] = mean
	return z
	
def customcolormaps(cmap):
	# Original
	"""cdict = {'red': [[0.0,  0.0, 0.0],
					[0.5,  1.0, 1.0],
					[1.0,  1.0, 1.0]],
			'green':[[0.0,  0.0, 0.0],
					[0.25, 0.0, 0.0],
					[0.75, 1.0, 1.0],
					[1.0,  1.0, 1.0]],
			'blue': [[0.0,  0.0, 0.0],
					[0.5,  0.0, 0.0],
					[1.0,  1.0, 1.0]]}"""
	
	# Red-Green-Blue		
	if cmap == 'RGB':				
		cdict = {'red': [[0.0,  1.0, 1.0],
						[1/3,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'green':[[0.0,  0.0, 0.0],
						[1/3, 0.0, 1.0],
						[2/3, 1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'blue': [[0.0,  0.0, 0.0],
						[2/3,  0.0, 1.0],
						[1.0,  1.0, 1.0]]}
	
	# White-Green
	if cmap == 'wg':
		cdict = {'red': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'green':[[0.0,  1.0, 1.0],
						[0.5, 1.0, 0.4],
						[1.0,  0.4, 0.4]],
				'blue': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.0],
						[1.0,  0.0, 0.0]]}
	
	# Colormap for 2dplottest				
	if cmap == '2dplottest':
		cdict = {'red': [[0.0, 1.0, 1.0],
						[0.1,  1.0, 0.9],
						[0.25, 0.9, 0.75],
						[0.5,  0.75, 0.5],
						[0.75, 0.5, 0.25],
						[0.9,  0.25, 0.1],
						[1.0,  0.1, 0.0]],
				'green':[[0.0, 1.0, 1.0],
						[0.1,  1.0, 1.0],
						[0.25, 1.0, 0.9],
						[0.5,  0.9, 0.7],
						[0.75, 0.7, 0.5],
						[0.9,  0.5, 0.4],
						[1.0,  0.4, 0.4]],
				'blue': [[0.0, 1.0, 1.0],
						[0.1,  1.0, 0.9],
						[0.25, 0.9, 0.75],
						[0.5,  0.75, 0.5],
						[0.75, 0.5, 0.25],
						[0.9,  0.25, 0.1],
						[1.0,  0.1, 0.0]]}
						
	# Testcmap		
	if cmap == 'QBOcmap':
		cdict = {'red': [[0.0,  0.0, 0.1],
						[0.125, 0.1, 0.2],
						[0.25, 0.2, 0.4],
						[0.375, 0.4, 0.7],
						[0.5, 0.7, 1.0],
						[0.625, 1.0, 1.0],
						[0.75, 1.0, 0.8],
						[0.875, 0.8, 0.4],
						[1.0,  0.4, 0.0]],
				'green':[[0.0,  0.0, 0.2],
						[0.125, 0.2, 0.4],
						[0.25, 0.4, 0.6],
						[0.375, 0.6, 0.8],
						[0.5, 0.8, 0.6],
						[0.625, 0.6, 0.2],
						[0.75, 0.2, 0.0],
						[0.875, 0.0, 0.0],
						[1.0,  0.0, 0.0]],
				'blue': [[0.0,  0.0, 0.3],
						[0.125, 0.3, 0.6],
						[0.25, 0.6, 0.8],
						[0.375, 0.8, 0.9],
						[0.5, 0.9, 0.6],
						[0.625, 0.6, 0.2],
						[0.75, 0.2, 0.0],
						[0.875, 0.0, 0.0],
						[1.0,  0.0, 0.0]]}
	
	# Testcmap		
	if cmap == 'QBOcmap2':
		cdict = {'red': [[0.0,  0.0, 0.1],
						[1/12, 0.1, 0.2],
						[2/12, 0.2, 0.4],
						[3/12, 0.4, 0.7],
						[4/12, 0.7, 0.85],
						[5/12, 0.85, 0.92],
						[6/12, 0.92, 1.0],
						[7/12, 1.0, 1.0],
						[8/12, 1.0, 1.0],
						[9/12, 1.0, 1.0],
						[10/12, 1.0, 0.8],
						[11/12, 0.8, 0.4],
						[1.0,  0.4, 0.0]],
				'green':[[0.0,  0.0, 0.2],
						[1/12, 0.2, 0.4],
						[2/12, 0.4, 0.6],
						[3/12, 0.6, 0.8],
						[4/12, 0.8, 0.9],
						[5/12, 0.9, 0.95],
						[6/12, 0.95, 0.9],
						[7/12, 0.9, 0.8],
						[8/12, 0.8, 0.6],
						[9/12, 0.6, 0.2],
						[10/12, 0.2, 0.0],
						[11/12, 0.0, 0.0],
						[1.0,  0.0, 0.0]],
				'blue': [[0.0,  0.0, 0.3],
						[1/12, 0.3, 0.6],
						[2/12, 0.6, 0.8],
						[3/12, 0.8, 0.9],
						[4/12, 0.9, 0.95],
						[5/12, 0.95, 0.975],
						[6/12, 0.975, 0.9],
						[7/12, 0.9, 0.8],
						[8/12, 0.8, 0.6],
						[9/12, 0.6, 0.2],
						[10/12, 0.2, 0.0],
						[11/12, 0.0, 0.0],
						[1.0,  0.0, 0.0]]}
						
	# Testcmap		
	if cmap == 'QBOcmap3':
		cdict = {'red': [[0.0,  0.0, 0.05],
						[1/22, 0.05, 0.1],
						[2/22, 0.1, 0.15],
						[3/22, 0.15, 0.2],
						[4/22, 0.2, 0.3],
						[5/22, 0.3, 0.4],
						[6/22, 0.4, 0.55],
						[7/22, 0.55, 0.7],
						[8/22, 0.7, 0.85],
						[9/22, 0.85, 0.92],
						[10/22, 0.92, 0.96],
						[11/22, 0.96, 1.0],
                        [12/22, 1.0, 1.0],
						[13/22, 1.0, 1.0],
						[14/22, 1.0, 1.0],
						[15/22, 1.0, 1.0],
						[16/22, 1.0, 1.0],
						[17/22, 1.0, 0.9],
						[18/22, 0.9, 0.8],
						[19/22, 0.8, 0.6],
						[20/22, 0.6, 0.4],
						[21/22, 0.4, 0.2],
						[1.0, 0.2, 0.0]],
				'green':[[0.0,  0.0, 0.1],
						[1/22, 0.1, 0.2],
						[2/22, 0.2, 0.3],
						[3/22, 0.3, 0.4],
						[4/22, 0.4, 0.5],
						[5/22, 0.5, 0.6],
						[6/22, 0.6, 0.7],
						[7/22, 0.7, 0.8],
						[8/22, 0.8, 0.9],
						[9/22, 0.9, 0.95],
						[10/22, 0.95, 0.975],
						[11/22, 0.975, 0.95],
                        [12/22, 0.95, 0.9],
						[13/22, 0.9, 0.8],
						[14/22, 0.8, 0.6],
						[15/22, 0.6, 0.4],
						[16/22, 0.4, 0.2],
						[17/22, 0.2, 0.1],
						[18/22, 0.1, 0.0],
						[19/22, 0.0, 0.0],
						[20/22, 0.0, 0.0],
						[21/22, 0.0, 0.0],
						[1.0, 0.0, 0.0]],
				'blue': [[0.0,  0.0, 0.15],
						[1/22, 0.15, 0.3],
						[2/22, 0.3, 0.45],
						[3/22, 0.45, 0.6],
						[4/22, 0.6, 0.7],
						[5/22, 0.7, 0.8],
						[6/22, 0.8, 0.85],
						[7/22, 0.85, 0.9],
						[8/22, 0.9, 0.925],
						[9/22, 0.95, 0.975],
						[10/22, 0.975, 0.988],
						[11/22, 0.988, 0.95],
						[12/22,  0.95, 0.9],
						[13/22, 0.9, 0.8],
						[14/22, 0.8, 0.6],
						[15/22, 0.6, 0.4],
						[16/22, 0.4, 0.2],
						[17/22, 0.2, 0.1],
						[18/22, 0.1, 0.0],
						[19/22, 0.0, 0.0],
						[20/22, 0.0, 0.0],
						[21/22, 0.0, 0.0],
						[1.0, 0.0, 0.0]]}
						
	# Whitehatches
	if cmap == "whitehatches":
		cdict = {'red': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'green':[[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'blue': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'alpha':[[0.0, 1.0, 1.0],
						[0.5, 1.0, 0.0],
						[1.0, 0.0, 0.0]]}
	
	# Grayhatches		
	if cmap == 'grayhatches':
		cdict = {'red': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'green':[[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'blue': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.6],
						[1.0,  0.6, 0.6]],
				'alpha':[[0.0, 0.0, 0.0],
						[0.5, 0.0, 1.0],
						[1.0, 1.0, 1.0]]}
						
	# Blackhatches		
	if cmap == 'blackhatches':
		cdict1 = {'red': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'green':[[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'blue': [[0.0,  1.0, 1.0],
						[0.5,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'alpha':[[0.0, 0.0, 0.0],
						[0.5, 0.0, 1.0],
						[1.0, 1.0, 1.0]]}
								
		cdict = cdict1.copy()
		# ~ cdict['alpha'] = ((0.0, 1.0, 1.0),
						 # ~ (0.5, 1.0, 0.0),
						 # ~ (1.0, 0.0, 0.0))
		
	
	# Testcmap		
	if cmap == 'testcmap':
		cdict = {'red': [[0.0,  1.0, 1.0],
						[1/3,  1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'green':[[0.0,  0.0, 0.0],
						[1/3, 0.0, 1.0],
						[2/3, 1.0, 0.0],
						[1.0,  0.0, 0.0]],
				'blue': [[0.0,  0.0, 0.0],
						[2/3,  0.0, 1.0],
						[1.0,  1.0, 1.0]]}
	return cdict

def lathemi(x):
	if x > 0:
		hemisphere = "N"
	elif x < 0:
		hemisphere = "S"
	else:
		hemisphere = ""
	return hemisphere
	
def lonhemi(x):
	if x > 0:
		hemisphere = "E"
	elif x < 0:
		hemisphere = "W"
	else:
		hemisphere = ""
	return hemisphere
