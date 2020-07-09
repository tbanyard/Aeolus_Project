========================================================================
These programs convert the Aeolus data in DBL files into the netCDF4 
format, primarily for use on the PC UBPC-2027. A header ncdump for the
files created by this program is given at the bottom of this README.md
file. The following files can be found in this directory:
========================================================================
Aeolus_ncconvert_ubpc_vx.x.py| Converts DBL files into netCDF34 format
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
Aeolus_wind_variance_vx.x.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Fixing array operations and indexing issues-------------------
---v2.0---[1st PhD CODE REVAMP]-----------------------------------------
------------------------------------------------------------------------
========================================================================
Reads .DBL files unzipped from TGZ files from the Aeolus database and 
converts them into .nc files for further use. See below for file format.
========================================================================


========================================================================
$ ncdump -h AE_2B_2020-06-02_015135.nc
========================================================================
dimensions:
	time = 27606 ;
	RG = 710 ;
variables:
	double time(time) ;
		time:standard_name = "time" ;
		time:long_name = "time" ;
		time:units = "seconds since 2000-01-01 00:00:00" ;
	double lon(time) ;
		lon:standard_name = "longitude" ;
		lon:long_name = "Longitude" ;
		lon:units = "degree_east" ;
	double lat(time) ;
		lat:standard_name = "latitude" ;
		lat:long_name = "Latitude" ;
		lat:units = "degree_north" ;
	double alt(time) ;
		alt:standard_name = "altitude" ;
		alt:long_name = "Altitude" ;
		alt:units = "m" ;
	double Rayleigh_HLOS_wind_speed(time) ;
		Rayleigh_HLOS_wind_speed:standard_name = "wind_speed" ;
		Rayleigh_HLOS_wind_speed:long_name = "Rayleigh_Horizontal_Line_of_Sight_Wind_speed" ;
		Rayleigh_HLOS_wind_speed:units = "cm s-1" ;
	double LOS_azimuth(time) ;
		LOS_azimuth:standard_name = "los_azimuth" ;
		LOS_azimuth:long_name = "Line of sight azimuth angle of the target-to-satellite pointing vector" ;
		LOS_azimuth:notes = "Measured in degrees from north" ;
		LOS_azimuth:units = "deg" ;
	double Zonal_wind_projection(time) ;
		Zonal_wind_projection:standard_name = "zonal_wind_projection" ;
		Zonal_wind_projection:long_name = "Zonal projection of the HLOS wind" ;
		Zonal_wind_projection:units = "cm s-1" ;
	double Meridional_wind_projection(time) ;
		Meridional_wind_projection:standard_name = "meridional_wind_projection" ;
		Meridional_wind_projection:long_name = "Meridional projection of the HLOS wind" ;
		Meridional_wind_projection:units = "cm s-1" ;
	double Satellite_Velocity(time) ;
		Satellite_Velocity:standard_name = "satellite_velocity" ;
		Satellite_Velocity:long_name = "Line of sight velocity of the satellite" ;
		Satellite_Velocity:units = "m s-1" ;
	double QC_Flag_Both(time) ;
		QC_Flag_Both:standard_name = "qc_flag_both" ;
		QC_Flag_Both:long_name = "Binary flag corresponding to both QC filters" ;
		QC_Flag_Both:notes = "The following quality controls have been applied to this variable: \t1) Rayleigh Cloudy data removed, leaving only Rayleigh Clear. 2) L2B_Rayleigh_Hlos_Error_Estimate must not equal 1.7e+38." ;
		QC_Flag_Both:key = "1 = GOOD, 0 = BAD" ;
		QC_Flag_Both:units = "unitless" ;
	double QC_Flag_ObsType(time) ;
		QC_Flag_ObsType:standard_name = "qc_flag_obstype" ;
		QC_Flag_ObsType:long_name = "Binary flag corresponding to only the observation type QC filter" ;
		QC_Flag_ObsType:notes = "The following quality controls have been applied to this variable: \t1) Rayleigh Cloudy data removed, leaving only Rayleigh Clear." ;
		QC_Flag_ObsType:key = "1 = GOOD (Clear Sky), 0 = BAD (Cloudy)" ;
		QC_Flag_ObsType:units = "unitless" ;
	double QC_Flag_HLOSErr(time) ;
		QC_Flag_HLOSErr:standard_name = "qc_flag_hloserr" ;
		QC_Flag_HLOSErr:long_name = "Binary flag corresponding to only the L2B Rayleigh HLOS Error QC filter" ;
		QC_Flag_HLOSErr:notes = "The following quality controls have been applied to this variable: \t1) L2B_Rayleigh_Hlos_Error_Estimate must not equal 1.7e+38." ;
		QC_Flag_HLOSErr:key = "1 = GOOD, 0 = BAD" ;
		QC_Flag_HLOSErr:units = "unitless" ;
	double RG(RG) ;
		RG:standard_name = "rayleigh_grouping" ;
		RG:long_name = "Rayleigh_Grouping" ;
		RG:units = "unitless" ;

// global attributes:
		:title = "Aeolus HLOS Rayleigh Wind Data" ;
		:contact = "T. P. Banyard, tpb38@bath.ac.uk" ;
		:institution = "University of Bath, Claverton Down, Bath, BA2 7AY, United Kingdom" ;
		:Aeolus_data_source = "https://aeolus-ds.eo.esa.int" ;
		:Date_of_creation = "02 Apr 2020" ;
}

