========================================================================
These programs loop through the entire Aeolus dataset and pick out those
orbital sections that pass through any specified lat/lon box. Various
filters are applied and ERA5 data is interpolated onto the Aeolus track
to calculate wind perturbations. The following files can be found in
this directory:
========================================================================
Aeolus_ncload_vx.x.py| Produces perturbation plots from Aeolus' database
Aeolus_ncload_ubpc-2027_vx.x.py| As above but for computer UBPC-2027
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
Aeolus_ncload_vx.x.py and Aeolus_QBO_Feb20_vx.x.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Looping through dataset and plotting all orbits that travel---
----------through the 40-80 box.----------------------------------------
---v1.2---3.2.20-Updates------------------------------------------------
---v1.3---Deals with multiple orbits------------------------------------
----------simple_contourf-----------------------------------------------
---v1.4---pcolormesh/imshow---------------------------------------------
---v1.5---Additional subplot with lat/lon grid showing satellite track--
---v1.6---23.03.20-Updates: Applying S-G Filter and 2D Boxcar on data---
---v1.7---Experimenting with 500m bins, with interpolation, different---
----------scheme for dealing with NaNs etc.-----------------------------
---v1.8---ERA5 interpolated onto Aeolus track---------------------------
---v1.9---Testing New NC_FullQC files, smoothing ERA5 first (in u and v)
---v1.10--Consolidating code and tidying up-----------------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
plots orbital cross-sections of the wind perturbations
========================================================================


