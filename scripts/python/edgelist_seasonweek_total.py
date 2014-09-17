#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 9/12/14
###Function: network generation for season level spatial correlation networks for all ages and service places.

###Import data: 

###Command Line: python 
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

### import list of non-continental zip3s to exclude from analysis ###
noncontzipsin = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/R_export/non_continental_zip3s.csv', 'r')
noncontzipsin.readline()
noncont = csv.reader(noncontzipsin, delimiter=',')
for item in noncont:
	zip3 = str(item[0])
	excl_zips.append(zip3)

### import incidence data ###
for snum in seasons:
	basefile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/fluseasonweeks_zip3_incid_S%s.csv' %(snum)
	filein = open(basefile,'r')
	d_week_weeknum, d_weekZip3_iliVisit, d_zip3_pop = ng.import_seasonweek_data(filein)
	d_zip3_iliTS, d_zip3_incidTS = ng.seasonweek_timeseries(d_weekZip3_iliVisit, d_zip3_pop)
	
	# create zip3 network edgelist (includes all zip3s)
	d_zip3pair_crosscorr = ng.timeseries_crosscorrelations(d_zip3_incidTS)
	d_zip3pairFilt_crosscorr, thresh = ng.crosscorrelation_edgelist(d_zip3pair_crosscorr, **kwargs_threshold)

	# create flows list (includes only continental zip3s)
	d_zip3_iliTS_cont, d_zip3_incidTS_cont = ng.exclude_noncontinental_zip3s(d_zip3_iliTS, d_zip3_incidTS, excl_zips)
	d_zip3pair_crosscorr_cont = ng.timeseries_crosscorrelations(d_zip3_incidTS_cont)
	d_zip3pairFilt_crosscorr_cont, thresh = ng.crosscorrelation_edgelist(d_zip3pair_crosscorr_cont, **kwargs_threshold)

	threshold_type, value = kwargs_threshold.items()[0]
	# write edgelist to file
	edgelist_fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/Py_export/edgelists/seasonweek_total_edgelist_pearson_%s%s_S%s.csv' %(threshold_type, value, snum)
	ng.write_edgelist(d_zip3pairFilt_crosscorr, edgelist_fname)

	# write flows.csv to file
	flows_fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/jflowmap/input_files/flows/seasonweek_total_flows_pearson_%s%s_S%s.csv' %(threshold_type, value, snum)
	ng.write_flowscsv_jflowmap(d_zip3pairFilt_crosscorr_cont, flows_fname)