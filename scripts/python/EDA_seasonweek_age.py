#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/2/14
###Function: exploratory data analysis for season level spatial correlation networks for specific age groups (toddler, children, adult, elderly) and all service places.

###Import data: 

###Command Line: python EDA_seasonweek_age.py
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
seasons = ng.sp_seasons
kwargs_TSdata_method = ng.cp_TSdata_method_kwargs
TSdata, method = kwargs_TSdata_method.values()
agecodes = ng.sp_agegroupcodes

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

### plot distribution of pairwise zip3 cross-correlations ###

		# convert to numpy array to use isnan function
		cross_corr = np.array([d_zip3pair_crosscorr[key] for key in d_zip3pair_crosscorr])
		# draw histogram
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)
		n, bins, patches = ax1.hist(cross_corr[~np.isnan(cross_corr)], bins=50)
		ax1.set_xlabel('pairwise cross-correlation (S%s, %s)' %(snum, ac))
		ax1.set_ylabel('number of zip3 pairs')
		ax1.set_xlim([-0.8,1])
		plt.savefig('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/hist_pairwiseZip3_%s%s_S%s%s.png' %(TSdata, method, snum, ac), transparent=False, bbox_inches='tight', pad_inches=0)
		plt.close()
		# plt.show()

		# print to screen
		print '\tSeason %s, %s complete' %(snum, ac)


