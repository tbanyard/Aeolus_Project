Aeolus_initial_test_DEPRECATED.py       Initial look at the Aeolus dataset and      plotting of first results [INACTIVE]
Aeolus_timeseries_v1.0.py               Generic program which creates a timeseries using the timeseriesplot function from phdfunctions.py, modified for AEOLUS. [ACTIVE]
Aeolus_testscans_v1.0.py                Program written to produce a graphical plot of the satellite pass. [INACTIVE]
Aeolus_testscans_v1.1.py                Program written to produce a graphical plot of the satellite pass. [INACTIVE]
Aeolus_testscans_v1.2.py                Program written to produce a graphical plot of the satellite pass. [INACTIVE]
Aeolus_testscans_v1.3.py                Program written to produce a graphical plot of the satellite pass. This version is the main program to plot AEOLUS data from DBL files. [ACTIVE]
Aeolus_testscans_v1.3_andesfocus.py     Program written to produce a graphical plot of the satellite pass. This version is an adapted version of the main program, focussing on the andes. [ACTIVE]
Aeolus_ncconvert.py                     Program which reads DBL files and converts to NC files of key parameters. Default year is 2020, default month is January. Uses createAeolusnc, load_rayleigh_data from phdfunctions.py. [ACTIVE]
Aeolus_ncconvert_ubpc-2027.py           Same program as above but for use on ubpc-2027.
Aeolus_ncload.py                        Program which produces timeseries heatmaps from the NC files created by Aeolus_ncconvert.py. Uses timeseriesplot and find_nearest from phdfunctions.py. [INACTIVE]
Aeolus_ncload_v1.1.py                   Program which produces timeseries heatmaps from the NC files created by Aeolus_ncconvert.py. Uses timeseriesplot and find_nearest from phdfunctions.py. [INACTIVE]
Aeolus_ncload_ubpc-2027.py              Same as Aeolus_ncload_v1.1.py but for use on ubpc-2027. [INACTIVE]
Aeolus_ncload_v1.2.py                   Program which produces timeseries heatmaps from the NC files created by Aeolus_ncconvert.py. Uses the groupby function from itertools, and timeseriesplot and find_nearest from phdfunctions.py. [INACTIVE]
Aeolus_ncload_ubpc-2027_v2.py              Same as Aeolus_ncload_v1.2.py but for use on ubpc-2027. [INACTIVE]
Aeolus_ncload_v1.3.py                   Program which produces timeseries heatmaps from the NC files created by Aeolus_ncconvert.py. Uses the groupby function from itertools, and timeseriesplot and find_nearest from phdfunctions.py. [ACTIVE]
Aeolus_ncload_ubpc-2027_v3.py              Same as Aeolus_ncload_v1.3.py but for use on ubpc-2027. [ACTIVE]

OTHER FILES:
Files in the format YYYY-mm-DD_... correspond to data from that date, and produce plots with the same name in the corresponding Plots folder.
