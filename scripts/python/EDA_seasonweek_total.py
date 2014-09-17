#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 9/11/14
###Function: exploratory data analysis for season level spatial correlation networks for all ages and service places.

###Import data: 

###Command Line: python 
##############################################


### notes ###


### packages/modules ###
import csv
import matplotlib.pyplot as plt
import numpy as np

## local modules ##
import network_generation as ng

### data structures ###

### parameters ###
seasons = range(2,11)

### import data ###
for snum in seasons:
	basefile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/fluseasonweeks_zip3_incid_S%s.csv' %(snum)
	filein = open(basefile,'r')
	d_week_weeknum, d_weekZip3_iliVisit, d_zip3_pop = ng.import_seasonweek_data(filein)
	d_zip3_iliTS, d_zip3_incidTS = ng.seasonweek_timeseries(d_weekZip3_iliVisit, d_zip3_pop)
	d_zip3pair_crosscorr = ng.timeseries_crosscorrelations(d_zip3_incidTS)

### program ###

	# convert to numpy array to use isnan function
	cross_corr = np.array([d_zip3pair_crosscorr[key][0] for key in d_zip3pair_crosscorr])

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	n, bins, patches = ax1.hist(cross_corr[~np.isnan(cross_corr)], bins=50)
	ax1.set_xlabel('pairwise cross-correlation')
	ax1.set_ylabel('number of zip3 pairs')
	plt.savefig('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/hist_pairwiseZip3_pearson_S%s.png' %(snum), transparent=False, bbox_inches='tight', pad_inches=0)
	plt.close()
	# plt.show()


