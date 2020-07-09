========================================================================
These programs investigate the Aeolus representation of the 2019-2020
disruption of the QBO (Quasi-Biennial Oscillation). The following files
can be found in this directory:
========================================================================
Aeolus_QBO_timeseries_vx.x.py| Produces timeseries plots of disruption
Aeolus_QBO_Feb20_vx.x.py| Precursor programs to the above (made Feb2020)
Aeolus_QBO_CW_EQ_newbins| Produces cartographical plots of disruption
Aeolus_QBO_TestProfiles_vx.x.py| Profiles of early June 2020 test orbits
========================================================================

Details for each can be found below:
________________________________________________________________________
For any extraneous modules not explicitly imported in these programs,
please use the .pythonstartup file provided. Either export this on your
own system, or copy the imported modules manually. This way, a
comprehensive list of commonly imported modules can be kept separate
outside each individual program.
________________________________________________________________________

========================================================================
Aeolus_QBO_timeseries_vx.x.py and Aeolus_QBO_Feb20_vx.x.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Ascending/Descending Node Split-------------------------------
---v1.2---Daily_Means and creating .mat file for Neil-------------------
---v1.3---Plotting my own version of this plot extended out to March----
---v1.4---Using generated netCDF file from v1.3-------------------------
---v1.5---Editing and improving plot------------------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
generates zonal mean U plot at equator
========================================================================


========================================================================
Aeolus_QBO_CW_EQ_newbins
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Testing, still with pcolormesh only---------------------------
---v1.2---Improving plot quality----------------------------------------
---v1.3---Creating 2D Test figure to see data gaps----------------------
---v1.4---Improving quality of maps-------------------------------------
---v1.5---Looping over height levels and days---------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
------------------------------------------------------------------------
========================================================================
Reads .mat files converted from .nc files from the Aeolus database and 
produces cartographical plots
========================================================================


========================================================================
Aeolus_QBO_TestProfiles_vx.x.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Trimmed to look at the four equator passes of Jun 3rd---------
---v1.2---Completed subplots with profile, fixing code------------------
---v1.3---Focusing on one or two Singapore profiles---------------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
colocates with one or two Singapore radiosonde profiles for testing
========================================================================


