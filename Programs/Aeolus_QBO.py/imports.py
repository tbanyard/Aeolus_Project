"""This file contains all of the modules imported into these programs
in addition to those found in the .pythonstartup file."""

# Matplotlib
import matplotlib.dates as dates
from matplotlib import colorbar
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.contour import ClabelText

# Scipy
import scipy.ndimage as ndimage
from scipy.signal import savgol_filter
from scipy.interpolate import griddata, interp2d
from scipy.io import savemat

# Itertools
from itertools import groupby

# Xarray
import xarray as xr

# Dask
from dask import array as da
from dask import delayed, visualize
from dask.diagnostics import ProgressBar
pbar = ProgressBar()
from dask.distributed import Client, LocalCluster
# cluster = LocalCluster(n_workers=2, processes=False, threads_per_worker=1, memory_limit='7GB', scheduler_port=0)
# client = Client(cluster)
# client = Client(processes=False)

from dask.distributed import performance_report

# Multiprocessing
import multiprocessing

# Warnings
import warnings
warnings.filterwarnings('ignore')
