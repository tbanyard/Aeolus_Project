#!/usr/bin/env python3
"""
Interpolate ERA5 onto Aeolus profiles from .nc files for each dataset
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
.nc files from MARS at the ECMWF for ERA5
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import matplotlib.text as text
import os
os.putenv('CODA_DEFINITION',
'/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.io import savemat
from scipy.signal import savgol_filter
import scipy.ndimage as ndimage
from itertools import groupby
from mpl_toolkits.basemap import Basemap
import random

# Import from functions file
import sys
sys.path.append('/home/tpb38/PhD/Bath/')
sys.path.append('/home/tpb38/PhD/Bath/Aeolus_Project/Programs/')
from phdfunctions import *
from functions import *

# Change current working directory to parent directory
os.chdir('..')

# Iterate through all. Aeolus nc files
"""Find directory and read netCDF data"""
strdirectory = '/home/tpb38/PhD/Bath/Aeolus/NC4/'
