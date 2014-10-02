#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/2/14
###Function: network generation for season level spatial correlation networks for specific age groups (toddler, children, adult, elderly) and all service places.

###Import data: 

###Command Line: python edgelist_seasonweek_age.py
##############################################


### notes ###


### packages/modules ###
import csv
import numpy as np

## local modules ##
import network_generation as ng

### data structures ###
excl_zips = []

### parameters ###
seasons = ng.sp_seasons
kwargs_threshold = ng.cp_threshold_kwargs
kwargs_TSdata_method = ng.cp_TSdata_method_kwargs
threshold_type, value = kwargs_threshold.items()[0]
TSdata, method = kwargs_TSdata_method.values()
agecodes = ng.sp_agegroupcodes

### import list of non-continental zip3s to exclude from analysis ###
noncontzipsin = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/R_export/non_continental_zip3s.csv', 'r')
noncontzipsin.readline()
noncont = csv.reader(noncontzipsin, delimiter=',')
for item in noncont:
	zip3 = str(item[0])
	excl_zips.append(zip3)

### import and process pop data ###
## NB. pop data includes all years and age groups
popfile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/allyears_zip3_popAge.csv'
popfilein = open(popfile, 'r')
# import pop data
d_yearZip3Smallage_pop = ng.import_yearZip3_popAge_data(popfilein)
# process pop data to large age bin 
d_yearZip3Age_pop = ng.bin_smallage_age_pop(d_yearZip3Smallage_pop)

### import and process ILI data ###
for snum in seasons:
	## NB. ILI data includes one flu season and age groups
	basefile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/fluseasonweeks_zip3_incidAge_S%s.csv' %(snum)
	filein = open(basefile,'r')
	d_week_weeknum, d_weekZip3Smallage_ili = ng.import_seasonweekAge_data(filein)
	# process ILI data to large age bin 
	d_weekZip3Age_ili = ng.bin_smallage_age_ili(d_weekZip3Smallage_ili)

	# separate specific age bin for ILI and pop data in preparation for seasonweek_timeseries function
	for ac in agecodes:
		d_weekZip3_iliVisit, d_zip3_pop = ng.prep_age_seasonweek_timeseries(d_weekZip3Age_ili, d_yearZip3Age_pop, ac, snum)
		d_zip3_iliTS, d_zip3_incidTS = ng.seasonweek_timeseries(d_weekZip3_iliVisit, d_zip3_pop)
		d_zip3pair_crosscorr = ng.timeseries_crosscorrelations(d_zip3_incidTS, **kwargs_TSdata_method)
		d_zip3pairFilt_crosscorr, thresh = ng.crosscorrelation_edgelist(d_zip3pair_crosscorr, **kwargs_threshold)

		# create flows list (includes only continental zip3s)
		d_zip3_iliTS_cont, d_zip3_incidTS_cont = ng.exclude_noncontinental_zip3s(d_zip3_iliTS, d_zip3_incidTS, excl_zips)
		d_zip3pair_crosscorr_cont = ng.timeseries_crosscorrelations(d_zip3_incidTS_cont, **kwargs_TSdata_method)
		d_zip3pairFilt_crosscorr_cont, thresh = ng.crosscorrelation_edgelist(d_zip3pair_crosscorr_cont, **kwargs_threshold)

		# write edgelist to file
		edgelist_fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/Py_export/edgelists/seasonweek_age_edgelist_%s%s_%s%s_S%s%s.csv' %(TSdata, method, threshold_type, value, snum, ac)
		ng.write_edgelist(d_zip3pairFilt_crosscorr, edgelist_fname)

		# write flows.csv to file
		flows_fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/jflowmap/input_files/flows/seasonweek_age_flows_%s%s_%s%s_S%s%s.csv' %(TSdata, method, threshold_type, value, snum, ac)
		ng.write_flowscsv_jflowmap(d_zip3pairFilt_crosscorr_cont, flows_fname)