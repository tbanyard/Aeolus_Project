#!/usr/bin/env python3
"""
Aeolus data conversion to netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .DBL files downloaded from the Aeolus database and produces netCDF
files of key parameters from the datasets
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import os 
import errno
import multiprocessing
from datetime import timedelta, datetime
import time

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from functions import createAeolusnc, load_rayleigh_data

# Change current working directory to parent directory
os.chdir('..')

def eachmonth(MM):
	"""
	This function...
	"""
	# Aeolus directory
	parent = '/home/tpb38/PhD/Bath/Aeolus/'
	# NetCDF file save directory
	nc_dir = 'NC/'

	# Set format for month
	if MM < 10:
		MM = '0' + str(MM)
	else:
		MM = str(MM)
		
	# Year
	YYYY = 2020
	
	# Set up directory format
	datetag = str(YYYY) + '/' + str(MM) + '/'
	strdirectory = parent + 'DATA2/' + datetag
	directory = os.fsencode(strdirectory)
	
	# Loop through directory
	for file in os.listdir(directory):
		
		print("\n===========================================================")
		
		# Program Timing (Time taken to get through one file)	
		startTime = datetime.now()
		
		# Setting filename for dataload
		filename = os.fsdecode(file)
		print(str(filename), '\n')

		# Set dbl link
		dbl = strdirectory + str(filename)

		# Scan for Day of Month
		dd = str(str(filename)[25:27])
		
		# Scan for Hour of Day
		HH = str(str(filename)[28:30])
		
		# Scan for Minute of Hour
		mm = str(str(filename)[30:32])
		
		# Scan for Second of Minute
		ss = str(str(filename)[32:34])
		
		print('YYYY-MM-dd HH:mm:ss = ', YYYY, '-', MM, '-', dd, ' ', HH, ':', mm, ':', ss, '\n')
		
		ncfilename = 'AE_' + str(YYYY) + '-' + str(MM) + '-' + str(dd) + '_' + str(HH) + str(mm) + str(ss) + '.nc'
		sub = nc_dir + ncfilename
		outfile = parent + sub
		print(outfile)
		createAeolusnc(dbl, outfile)

if __name__ == '__main__':
	"""Enables program to be executed using multiple processes/cores"""
	startTime = datetime.now()
	processes = []
	for MM in range(1,13):
		p = multiprocessing.Process(target=eachmonth, args=(MM,))
		processes.append(p)
		p.start()
			
	for process in processes:
		process.join()
	
	# Time taken for the entire program
	duration = datetime.now() - startTime
	print('That took ', duration, ' seconds')
