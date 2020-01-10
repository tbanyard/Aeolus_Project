#!/usr/bin/env python3

# Imports
import os
os.putenv('CODA_DEFINITION', '/opt/anaconda3/envs/virtualenv/share/coda/definitions/AEOLUS-20191015.codadef')
import coda
import numpy

# Opening Files
pf = coda.open('AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.DBL')
pfhdr = coda.open('AE_OPER_ALD_U_N_2B_20191216T070123_20191216T083159_0001.HDR')

# Fetching Data
mie_wind_velocity = coda.fetch(pf, 'mie_hloswind', -1, 'windresult/mie_wind_velocity')
latitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/latitude_cog')
longitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/longitude_cog')
altitude = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/altitude_vcog')
date_time = coda.fetch(pf, 'mie_geolocation', -1, 'windresult_geolocation/datetime_cog')

# Retrieving Field Names
field_names = coda.get_field_names(pf)
print(field_names)
mie_hloswind_field_names = coda.get_field_names(pf, 'mie_hloswind', 0)
print(mie_hloswind_field_names)
mie_geolocation_field_names = coda.get_field_names(pf, 'mie_geolocation', 0)
print(mie_geolocation_field_names)
windresult_geolocation_field_names = coda.get_field_names(pf, 'mie_geolocation', 0, 'windresult_geolocation')
print(windresult_geolocation_field_names)

# Arrays
print(mie_wind_velocity.shape)
print(mie_wind_velocity)
