#!/usr/bin/env python3
"""
Aeolus Functions file
================================================================================
--------------------------------------------------------------------------------
----------------------------No version control----------------------------------
--------------------------------------------------------------------------------
================================================================================
Provides functions for all files for the Aeolus project
================================================================================
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


def load_hdr_tags(hdr):
	# Open HDR file
	pfhdr = coda.open(hdr)
	
	# Get field names and fetch data according to the ICD document:
	# L-2B/2C_I/O_Data_Definitions
	Earth_Explorer_Header = coda.fetch(pfhdr, 'Earth_Explorer_Header')
	
	# Fixed_Header
	Fixed_Header = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Fixed_Header')
	File_Name = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/File_Name')
	File_Description = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/File_Description')
	Notes = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Notes')
	Mission = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Mission')
	File_Class = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/File_Class')
	File_Type = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/File_Type')
	File_Version = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/File_Version')
	
	#	 Validity_Period:
	Validity_Period = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Validity_Period')
	Validity_Start = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Validity_Period/Validity_Start')
	Validity_Stop = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Validity_Period/Validity_Stop')
	
	#	 Source:
	Source = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Source')
	System = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Source/System')
	Creator = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Source/Creator')
	Creator_Version = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Source/Creator_Version')
	Creation_Date = coda.fetch(pfhdr, 'Earth_Explorer_Header/Fixed_Header/Source/Creation_Date')
	
	# Variable_Header
	Variable_Header = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header')
	
	#	 Main_Product_Header
	Main_Product_Header = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header')
	Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Product')
	Proc_Stage = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Proc_Stage')
	Ref_Doc = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Ref_Doc')
	Acquisition_Station = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Acquisition_Station')
	Proc_Center = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Proc_Center')
	Proc_Time = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Proc_Time')
	Software_Ver = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Software_Ver')
	Baseline = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Baseline')
	Sensing_Start = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Sensing_Start')
	Sensing_Stop = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Sensing_Stop')
	Phase = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Phase')
	Cycle = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Cycle')
	Rel_Orbit = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Rel_Orbit')
	Abs_Orbit = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Abs_Orbit')
	State_Vector_Time = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/State_Vector_Time')
	Delta_UT1 = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Delta_UT1')
	X_Position = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/X_Position')
	Y_Position = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Y_Position')
	Z_Position = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Z_Position')
	X_Velocity = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/X_Velocity')
	Y_Velocity = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Y_Velocity')
	Z_Velocity = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Z_Velocity')
	Vector_Source = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Vector_Source')
	Utc_Sbt_Time = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Utc_Sbt_Time')
	Sat_Binary_Time = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Sat_Binary_Time')
	Clock_Step = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Clock_Step')
	Leap_Utc = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Leap_Utc')
	Gps_Utc_Time_Difference = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Gps_Utc_Time_Difference')
	Leap_Sign = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Leap_Sign')
	Leap_Err = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Leap_Err')
	Product_Err = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Product_Err')
	Tot_Size = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Tot_Size')
	Sph_Size = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Sph_Size')
	Num_Dsd = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Num_Dsd')
	Dsd_Size = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Dsd_Size')
	Num_Data_Sets = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Main_Product_Header/Num_Data_Sets')
	
	#	 Specific_Product_Header
	Specific_Product_Header = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header')
	Sph_Descriptor = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Sph_Descriptor')
	NumMeasurements = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumMeasurements')
	NumMieGroups = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumMieGroups')
	NumRayleighGroups = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumRayleighGroups')
	NumMieWindResults = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumMieWindResults')
	NumRayleighWindResults = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumRayleighWindResults')
	NumMieProfiles = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumMieProfiles')
	NumRayleighProfiles = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumRayleighProfiles')
	NumAMDprofiles = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/NumAMDprofiles')
	Intersect_Start_Lat = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Intersect_Start_Lat')
	Intersect_Start_Long = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Intersect_Start_Long')
	Intersect_Stop_Lat = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Intersect_Stop_Lat')
	Intersect_Stop_Long = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Intersect_Stop_Long')
	Sat_Track = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Sat_Track')
	
	#		 Counts
	List_of_Valid_Mie_Profile_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_Mie_Profile_Counts')
	Valid_Mie_Profile_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_Mie_Profile_Counts/Valid_Mie_Profile_Count')
	List_of_Valid_Rayleigh_Profile_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_Rayleigh_Profiles_Counts')
	Valid_Rayleigh_Profile_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_Rayleigh_Profiles_Counts/Valid_Rayleigh_Profile_Count')
	List_of_Invalid_Mie_Profile_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_Mie_Profile_Counts')
	Invalid_Mie_Profile_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_Mie_Profile_Counts/Invalid_Mie_Profile_Count')
	List_of_Invalid_Rayleigh_Profile_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_Rayleigh_Profile_Counts')
	Invalid_Rayleigh_Profile_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_Rayleigh_Profile_Counts/Invalid_Rayleigh_Profile_Count')
	Num_Profiles_Surface_Mie = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Num_Profiles_Surface_Mie')
	Num_Profiles_Surface_Ray = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/Num_Profiles_Surface_Ray')
	List_of_Valid_L2B_Mie_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2B_Mie_Wind_Counts')
	Valid_L2B_Mie_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2B_Mie_Wind_Counts/Valid_L2B_Mie_Wind_Count')
	List_of_Valid_L2B_Rayleigh_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2B_Rayleigh_Wind_Counts')
	Valid_L2B_Rayleigh_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2B_Rayleigh_Wind_Counts/Valid_L2B_Rayleigh_Wind_Count')
	List_of_Invalid_L2B_Mie_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2B_Mie_Wind_Counts')
	Invalid_L2B_Mie_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2B_Mie_Wind_Counts/Invalid_L2B_Mie_Wind_Count')
	List_of_Invalid_L2B_Rayleigh_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2B_Rayleigh_Wind_Counts')
	Invalid_L2B_Rayleigh_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2B_Rayleigh_Wind_Counts/Invalid_L2B_Rayleigh_Wind_Count')
	
	
	#		 L2C Products Only: [N.B. These will need adding to the return statement if used]
	# ~ List_of_Valid_L2C_Mie_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2C_Mie_Wind_Counts')
	# ~ Valid_L2C_Mie_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2C_Mie_Wind_Counts/Valid_L2C_Mie_Wind_Count')
	# ~ List_of_Valid_L2C_Rayleigh_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2C_Rayleigh_Wind_Counts')
	# ~ Valid_L2C_Rayleigh_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Valid_L2C_Rayleigh_Wind_Counts/Valid_L2C_Rayleigh_Wind_Count')
	# ~ List_of_Invalid_L2C_Mie_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2C_Mie_Wind_Counts')
	# ~ Invalid_L2C_Mie_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2C_Mie_Wind_Counts/Invalid_L2C_Mie_Wind_Count')
	# ~ List_of_Invalid_L2C_Rayleigh_Wind_Counts = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2C_Rayleigh_Wind_Counts')
	# ~ Invalid_L2C_Rayleigh_Wind_Count = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Invalid_L2C_Rayleigh_Wind_Counts/Invalid_L2C_Rayleigh_Wind_Count')
	
	#		 Dsds N.B. [For DSD Number add 1 to the python indexes used]
	List_of_Dsds = coda.get_field_names(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds')
	Dsd = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd')
	Meas_Map_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 0)
	Mie_Grouping_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 1)
	Rayleigh_Grouping_Map = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 2)
	Mie_Geolocation_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 3)
	Rayleigh_Geolocation_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 4)
	AMD_Product_Confid_Data_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 5)
	Meas_Product_Confid_Data_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 6)
	Mie_Wind_Product_Conf_Data_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 7)
	Rayl_Wind_Prod_Conf_Data_ADS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 8)
	Mie_Wind_MDS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 9)
	Rayleigh_Wind_MDS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 10)
	Mie_Profile_MDS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 11)
	Rayleigh_Profile_MDS = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 12)
	Aeolus_Level_1B_Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 13)
	Aux_Met_Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 14)
	Aeolus_RBC = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 15)
	Clim_Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 16)
	Cal_Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 17)
	Level_2B_Proc_Params = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 18)
	
	#		 L2C Products Only: [N.B. These will need adding to the return statement if used]
	# ~ Aeolus_Level_2B_Product = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 19)
	# ~ Level_2C_Proc_Params = coda.fetch(pfhdr, 'Earth_Explorer_Header/Variable_Header/Specific_Product_Header/List_of_Dsds/Dsd', 20)
	
	return Earth_Explorer_Header, Fixed_Header, File_Name, \
	File_Description, Notes, Mission, File_Class, File_Type, \
	File_Version, Validity_Period, Validity_Start, Validity_Stop, \
	Source, System, Creator, Creator_Version, Creation_Date, \
	Variable_Header, Main_Product_Header, Product, Proc_Stage, \
	Ref_Doc, Acquisition_Station, Proc_Center, Proc_Time, \
	Software_Ver, Baseline, Sensing_Start, Sensing_Stop, Phase, \
	Cycle, Rel_Orbit, Abs_Orbit, State_Vector_Time, Delta_UT1, \
	X_Position, Y_Position, Z_Position, X_Velocity, Y_Velocity, \
	Z_Velocity, Vector_Source, Utc_Sbt_Time, Sat_Binary_Time, \
	Clock_Step, Leap_Utc, Gps_Utc_Time_Difference, Leap_Sign, \
	Leap_Err, Product_Err, Tot_Size, Sph_Size, Num_Dsd, Dsd_Size, \
	Num_Data_Sets, Specific_Product_Header, Sph_Descriptor, \
	NumMeasurements, NumMieGroups, NumRayleighGroups, \
	NumMieWindResults, NumRayleighWindResults, NumMieProfiles, \
	NumRayleighProfiles, NumAMDprofiles, Intersect_Start_Lat, \
	Intersect_Start_Long, Intersect_Stop_Lat, Intersect_Stop_Long, \
	Sat_Track, List_of_Valid_Mie_Profile_Counts, \
	Valid_Mie_Profile_Count, List_of_Valid_Rayleigh_Profile_Counts, \
	Valid_Rayleigh_Profile_Count, List_of_Invalid_Mie_Profile_Counts, \
	Invalid_Mie_Profile_Count, List_of_Invalid_Rayleigh_Profile_Counts, \
	Invalid_Rayleigh_Profile_Count, Num_Profiles_Surface_Mie, \
	Num_Profiles_Surface_Ray, List_of_Valid_L2B_Mie_Wind_Counts, \
	Valid_L2B_Mie_Wind_Count, List_of_Valid_L2B_Rayleigh_Wind_Counts, \
	Valid_L2B_Rayleigh_Wind_Count, List_of_Invalid_L2B_Mie_Wind_Counts, \
	Invalid_L2B_Mie_Wind_Count, List_of_Invalid_L2B_Rayleigh_Wind_Counts, \
	Invalid_L2B_Rayleigh_Wind_Count, List_of_Dsds, Dsd, Meas_Map_ADS, \
	Mie_Grouping_ADS, Rayleigh_Grouping_Map, Mie_Geolocation_ADS, \
	Rayleigh_Geolocation_ADS, AMD_Product_Confid_Data_ADS, \
	Meas_Product_Confid_Data_ADS, Mie_Wind_Product_Conf_Data_ADS, \
	Rayl_Wind_Prod_Conf_Data_ADS, Mie_Wind_MDS, Rayleigh_Wind_MDS, \
	Mie_Profile_MDS, Rayleigh_Profile_MDS, Aeolus_Level_1B_Product, \
	Aux_Met_Product, Aeolus_RBC, Clim_Product, Cal_Product, \
	Level_2B_Proc_Params
	
def load_dbl_tags(dbl):
	# Open DBL file
	pf = coda.open(dbl)
	
	# Get field names and fetch data according to the document:
	# L-2B/2C_I/O_Data_Definitions
	
	# Main_Product_Header
	Product = coda.fetch(pf, 'mph/Product')
	Proc_Stage = coda.fetch(pf, 'mph/Proc_Stage')
	Ref_Doc = coda.fetch(pf, 'mph/Ref_Doc')
	Acquisition_Station = coda.fetch(pf, 'mph/Acquisition_Station')
	Proc_Center = coda.fetch(pf, 'mph/Proc_Center')
	Proc_Time = coda.fetch(pf, 'mph/Proc_Time')
	Software_Ver = coda.fetch(pf, 'mph/Software_Ver')
	Baseline = coda.fetch(pf, 'mph/Baseline')
	Sensing_Start = coda.fetch(pf, 'mph/Sensing_Start')
	Sensing_Stop = coda.fetch(pf, 'mph/Sensing_Stop')
	Phase = coda.fetch(pf, 'mph/Phase')
	Cycle = coda.fetch(pf, 'mph/Cycle')
	Rel_Orbit = coda.fetch(pf, 'mph/Rel_Orbit')
	Abs_Orbit = coda.fetch(pf, 'mph/Abs_Orbit')
	State_Vector_Time = coda.fetch(pf, 'mph/State_Vector_Time')
	Delta_UT1 = coda.fetch(pf, 'mph/Delta_UT1')
	X_Position = coda.fetch(pf, 'mph/X_Position')
	Y_Position = coda.fetch(pf, 'mph/Y_Position')
	Z_Position = coda.fetch(pf, 'mph/Z_Position')
	X_Velocity = coda.fetch(pf, 'mph/X_Velocity')
	Y_Velocity = coda.fetch(pf, 'mph/Y_Velocity')
	Z_Velocity = coda.fetch(pf, 'mph/Z_Velocity')
	Vector_Source = coda.fetch(pf, 'mph/Vector_Source')
	Utc_Sbt_Time = coda.fetch(pf, 'mph/Utc_Sbt_Time')
	Sat_Binary_Time = coda.fetch(pf, 'mph/Sat_Binary_Time')
	Clock_Step = coda.fetch(pf, 'mph/Clock_Step')
	Leap_Utc = coda.fetch(pf, 'mph/Leap_Utc')
	Gps_Utc_Time_Difference = coda.fetch(pf, 'mph/Gps_Utc_Time_Difference')
	Leap_Sign = coda.fetch(pf, 'mph/Leap_Sign')
	Leap_Err = coda.fetch(pf, 'mph/Leap_Err')
	Product_Err = coda.fetch(pf, 'mph/Product_Err')
	Tot_Size = coda.fetch(pf, 'mph/Tot_Size')
	Sph_Size = coda.fetch(pf, 'mph/Sph_Size')
	Num_Dsd = coda.fetch(pf, 'mph/Num_Dsd')
	Dsd_Size = coda.fetch(pf, 'mph/Dsd_Size')
	Num_Data_Sets = coda.fetch(pf, 'mph/Num_Data_Sets')
	
	# Specific_Product_Header
	Sph_Descriptor = coda.fetch(pf, 'sph/Sph_Descriptor')
	NumMeasurements = coda.fetch(pf, 'sph/NumMeasurements')
	NumMieGroups = coda.fetch(pf, 'sph/NumMieGroups')
	NumRayleighGroups = coda.fetch(pf, 'sph/NumRayleighGroups')
	NumMieWindResults = coda.fetch(pf, 'sph/NumMieWindResults')
	NumRayleighWindResults = coda.fetch(pf, 'sph/NumRayleighWindResults')
	NumMieProfiles = coda.fetch(pf, 'sph/NumMieProfiles')
	NumRayleighProfiles = coda.fetch(pf, 'sph/NumRayleighProfiles')
	NumAMDprofiles = coda.fetch(pf, 'sph/NumAMDprofiles')
	Intersect_Start_Lat = coda.fetch(pf, 'sph/Intersect_Start_Lat')
	Intersect_Start_Long = coda.fetch(pf, 'sph/Intersect_Start_Long')
	Intersect_Stop_Lat = coda.fetch(pf, 'sph/Intersect_Stop_Lat')
	Intersect_Stop_Long = coda.fetch(pf, 'sph/Intersect_Stop_Long')
	Sat_Track = coda.fetch(pf, 'sph/Sat_Track')
	
	#	 Counts
	Valid_Mie_Profile_Count = coda.fetch(pf, 'sph/Valid_Mie_Profile_Count')
	Valid_Rayleigh_Profile_Count = coda.fetch(pf, 'sph/Valid_Rayleigh_Profile_Count')
	Invalid_Mie_Profile_Count = coda.fetch(pf, 'sph/Invalid_Mie_Profile_Count')
	Invalid_Rayleigh_Profile_Count = coda.fetch(pf, 'sph/Invalid_Rayleigh_Profile_Count')
	Num_Profiles_Surface_Mie = coda.fetch(pf, 'sph/Num_Profiles_Surface_Mie')
	Num_Profiles_Surface_Ray = coda.fetch(pf, 'sph/Num_Profiles_Surface_Ray')
	Valid_L2B_Mie_Wind_Count = coda.fetch(pf, 'sph/Valid_L2B_Mie_Wind_Count')
	Valid_L2B_Rayleigh_Wind_Count = coda.fetch(pf, 'sph/Valid_L2B_Rayleigh_Wind_Count')
	Invalid_L2B_Mie_Wind_Count = coda.fetch(pf, 'sph/Invalid_L2B_Mie_Wind_Count')
	Invalid_L2B_Rayleigh_Wind_Count = coda.fetch(pf, 'sph/Invalid_L2B_Rayleigh_Wind_Count')
	
	# Dsds
	Dsd = coda.fetch(pf, 'dsd')
	# print((coda.fetch(pf, 'dsd'), 0)[0][0]) ~ Use this format to retrieve coda records
	
	# Measurement_Maps [These relate the L1B measurements to the L2B wind retrievals]
	Meas_Map = coda.fetch(pf, 'meas_map')
	# ~ Meas_Map_ADS_measurement = coda.fetch(pf, 'meas_map', m)
	Mie_Map_of_L1B_Meas_Used = coda.fetch(pf, 'meas_map', -1, 'mie_map_of_l1b_meas_used')
	"""Access via Mie_Map_of_L1B_Meas_Used[m][n], or:"""
	# ~ Mie_Map_of_L1B_Meas_Used_measurement = coda.fetch(pf, 'meas_map', -1, 'mie_map_of_l1b_meas_used', n)
	Rayleigh_Map_of_L1B_Meas_Used = coda.fetch(pf, 'meas_map', -1, 'rayleigh_map_of_l1b_meas_used')
	"""Access via Rayleigh_Map_of_L1B_Meas_Used[m][n], or:"""
	# ~ Rayleigh_Map_of_L1B_Meas_Used = coda.fetch(pf, 'meas_map', -1, 'rayleigh_map_of_l1b_meas_used', n)
	
	#	 Subrecords:
	# Bin = coda.fetch(pf, 'meas_map', -1, 'mie_map_of_l1b_meas_used', n)[m]
	# Which_L2B_Wind_id = coda.fetch(pf, 'meas_map', -1, 'mie_map_of_l1b_meas_used', n)[m][0]
	# Weight = coda.fetch(pf, 'meas_map', -1, 'mie_map_of_l1b_meas_used', n)[m][1]
	"""[Here n is in range(24) and m is in len(Mie_Map_of_L1B_Meas_Used) = NumMeasurements]"""
	
	# Mie and Rayleigh Grouping
	Mie_Grouping = coda.fetch(pf, 'mie_grouping')
	"""Access via Mie_Grouping[l][PI]"""
	Rayleigh_Grouping = coda.fetch(pf, 'rayleigh_grouping')
	"""Access via Rayleigh_Grouping[l][PI]"""
	"""[Here, l is the length of the grouping array and PI is the python index corresponding
	to the tags in Table 19 in the ICD]"""
	
	# Mie and Rayleigh Geolocation
	Mie_Geolocation = coda.fetch(pf, 'mie_geolocation')
	Rayleigh_Geolocation = coda.fetch(pf, 'rayleigh_geolocation')
	
	# Confidence Data
	AMD_Product_Confid_Data = coda.fetch(pf, 'amd_product_confid_data')
	Meas_Product_Confid_Data = coda.fetch(pf, 'meas_product_confid_data')
	Mie_Wind_Prod_Conf_Data = coda.fetch(pf, 'mie_wind_prod_conf_data')
	Rayleigh_Wind_Prod_Conf_Data = coda.fetch(pf, 'rayleigh_wind_prod_conf_data')
	
	# Mie_HLOS_Wind
	Mie_HLOS_Wind = coda.fetch(pf, 'mie_hloswind')
	Rayleigh_HLOS_Wind = coda.fetch(pf, 'rayleigh_hloswind')
	Mie_Profile = coda.fetch(pf, 'mie_profile')
	Rayleigh_Profile = coda.fetch(pf, 'rayleigh_profile')
	
	return Product, Proc_Stage,	Ref_Doc, Acquisition_Station, \
	Proc_Center, Proc_Time, Software_Ver, Baseline, Sensing_Start, \
	Sensing_Stop, Phase, Cycle, Rel_Orbit, Abs_Orbit, State_Vector_Time, \
	Delta_UT1, X_Position, Y_Position, Z_Position, X_Velocity, Y_Velocity, \
	Z_Velocity, Vector_Source, Utc_Sbt_Time, Sat_Binary_Time, \
	Clock_Step, Leap_Utc, Gps_Utc_Time_Difference, Leap_Sign, \
	Leap_Err, Product_Err, Tot_Size, Sph_Size, Num_Dsd, Dsd_Size, \
	Num_Data_Sets, Sph_Descriptor, NumMeasurements, NumMieGroups, \
	NumRayleighGroups, NumMieWindResults, NumRayleighWindResults, \
	NumMieProfiles,	NumRayleighProfiles, NumAMDprofiles, Intersect_Start_Lat, \
	Intersect_Start_Long, Intersect_Stop_Lat, Intersect_Stop_Long, \
	Sat_Track, Valid_Mie_Profile_Count, Valid_Rayleigh_Profile_Count, \
	Invalid_Mie_Profile_Count, Invalid_Rayleigh_Profile_Count, \
	Num_Profiles_Surface_Mie, Num_Profiles_Surface_Ray, \
	Valid_L2B_Mie_Wind_Count, Valid_L2B_Rayleigh_Wind_Count, \
	Invalid_L2B_Mie_Wind_Count, Invalid_L2B_Rayleigh_Wind_Count, \
	Dsd, Meas_Map, Mie_Map_of_L1B_Meas_Used, Rayleigh_Map_of_L1B_Meas_Used, \
	Mie_Grouping, Rayleigh_Grouping, Mie_Geolocation, Rayleigh_Geolocation, \
	AMD_Product_Confid_Data, Meas_Product_Confid_Data, Mie_Wind_Prod_Conf_Data, \
	Rayleigh_Wind_Prod_Conf_Data, Mie_HLOS_Wind, Rayleigh_HLOS_Wind, Mie_Profile, \
	Rayleigh_Profile
