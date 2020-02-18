#!/usr/bin/env python3
"""
Plot Aeolus data segment as pcolormesh or imagesc/imshow from npy file
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
------------------------------------------------------------------------
========================================================================
Plots the Aeolus data from an npy file for an orbital segment as either
a pcolormesh plot or an imagesc/imshow plot.
========================================================================
"""

"""-----------------------------------------------------------------------------
----------------------------------Imports---------------------------------------
-----------------------------------------------------------------------------"""

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
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from itertools import groupby

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import ncload

# ~ enterdirectory("test")
os.chdir('..')
os.chdir('..')
os.chdir('Plots')
print(os.getcwd())

"""-----------------------------------------------------------------------------
----------------------------------Program---------------------------------------
-----------------------------------------------------------------------------"""
# Choose pcolor or imshow
pc_or_im = 'im'

# Load data
data = np.load('data.npy', allow_pickle=True)
date_time = np.load('date_time.npy', allow_pickle=True)
x, y, z = data[0], data[1], data[2]
z = np.array(z, dtype='float64')

# Initialise figure
fig = plt.figure()
ax1 = fig.add_subplot(111)

if pc_or_im == 'pc':
	cs = plt.pcolormesh(x, y/1000, z, cmap='RdBu_r', vmin=-200, vmax=200)
	fixmydateaxis(ax1, date_time) # Uses my own function 'fixmydateaxis'

elif pc_or_im == 'im':
	# Plot imshow plot
	x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
	x_lims = dates.date2num(x_lims)
	y_lims = [20.5, -0.5]
	fixnanswithmean(z) # Uses my own function 'fixnanswithmean'
	cs = plt.imshow(z, aspect='auto', cmap='RdBu_r', extent=[x_lims[0],
		x_lims[1], y_lims[0], y_lims[1]], vmin=-200, vmax=200,
		interpolation='sinc')
	ax1.xaxis_date()
	date_form = dates.DateFormatter('%H:%M')
	ax1.xaxis.set_major_formatter(date_form)
	plt.gca().invert_yaxis()

# z array check
print(z)

# Fix labels and title
ax1.set_xlabel('Time')
ax1.set_ylabel('Altitude / km')
plt.title('Aeolus Orbit HLOS Rayleigh Wind Cross-section')

# Fix colorbar and axes
fig.colorbar(cs, cmap='RdBu_r', ax=ax1, orientation='horizontal',
	label='HLOS Rayleigh Wind Speed / ms-1')
ax1.set_yticks(np.arange(len(y)))

# Save figure
os.chdir('/home/tpb38/PhD/Bath/Aeolus_Project/Plots/dump')
plt.savefig('testbed.png',dpi=300)
