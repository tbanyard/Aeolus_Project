========================================================================
These programs pick out SABER orbital segments that pass through a
specified orbital box, apply S-G filtering and produce cross-sectional
profiles to colocate with Aeolus data. The following files
can be found in this directory:
========================================================================
SABER_ncload_vx.x.py| Produces timeseries plots of disruption
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
SABER_ncload_vx.x.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Splitting by time, into sections------------------------------
---v1.2---40-80 Andes box: Preparing code for this experiment.----------
----------Solved issue with the time formatting and overlapping times.--
---v1.3---Complete 40-80 Andes box.-------------------------------------
---v1.4---Satellite track subplot.--------------------------------------
---v2.0---23.03.20-Updates: Applying S-G Filtering----------------------
------------------------------------------------------------------------
========================================================================
Reads .nc files from the SABER database and generates profiles of the
kinetic temperature from teh SABER instrument, with a map for context.
========================================================================

