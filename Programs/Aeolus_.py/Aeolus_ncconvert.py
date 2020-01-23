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
from functions import createAeolusnc

# Change current working directory to parent directory
os.chdir('..')

def eachset(SS):
	"""
	This function...
	"""
	# Aeolus directory
	parent = '/home/tpb38/PhD/Bath/Aeolus/'
	# NetCDF file save directory
	nc_dir = 'NC/'
	
	directory = parent + 'DBL/'
	directory = os.fsencode(directory)		
	for file in os.listdir(directory):
		# Program Timing (Time taken to get through one file)	
		startTime = datetime.now()
		
		# Setting filename for dataload
		filename = os.fsdecode(file)
		print(str(filename), '\n')
		
	
	sub = nc_dir + 'example'
	outfile = parent + sub
	creataeolusnc(outfile, rayleigh_time, rayleigh_alt, rayleigh_lat,
	rayleigh_lon, rayleigh_wind)

if __name__ == '__main__':
	"""Enables program to be executed using multiple processes/cores"""
	startTime = datetime.now()
	processes = []
	for YYYY in range(1994,2021):
		p = multiprocessing.Process(target=eachyear, args=(YYYY,))
		processes.append(p)
		p.start()
			
	for process in processes:
		process.join()
	
	# Time taken for the entire program
	duration = datetime.now() - startTime
	print('That took ', duration, ' seconds')
