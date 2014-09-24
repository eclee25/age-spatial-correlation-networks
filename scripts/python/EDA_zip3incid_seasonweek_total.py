#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 9/23/14
###Function: exploratory data analysis for season level spatial correlation networks for all ages and service places. Incidence for each zip3 by season, for each HHS region.

###Import data: 

###Command Line: python 
##############################################


### notes ###


### packages/modules ###
import csv
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from scipy import stats

## local modules ##
import network_generation as ng

### data structures ###

### parameters ###
seasons = ng.sp_seasons
regions = ng.sp_HHSregions
fluweeklab = ng.sp_fluweeklabels

### import data ###
zip3RegionFile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/outside_source_data/allpopstat_zip3_season_cl.csv'
filein = open(zip3RegionFile, 'r')
d_zip3_HHS = ng.import_zip3_region_data(filein)

for snum in seasons:
	basefile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/fluseasonweeks_zip3_incid_S%s.csv' %(snum)
	filein = open(basefile,'r')
	d_week_weeknum, d_weekZip3_iliVisit, d_zip3_pop = ng.import_seasonweek_data(filein)
	d_zip3_iliTS, d_zip3_incidTS = ng.seasonweek_timeseries(d_weekZip3_iliVisit, d_zip3_pop)

### program ###

	
	for rnum in regions:

		# list of zip3s in a region
		ls_zip3 = [zip3 for zip3 in d_zip3_HHS if d_zip3_HHS[zip3]==rnum]
		
		# 9/23/14 plot zip3 incidence by HHS region
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)

		for zip3 in ls_zip3:
			ax1.plot([incid * 100000 for incid in d_zip3_incidTS[zip3]], '-', color='0.6')

		ax1.set_xlabel('week number')
		ax1.set_ylabel('incidence per 100,000')
		ax1.set_xticks(range(len(fluweeklab))[::5])
		ax1.set_xticklabels(fluweeklab[::5])
		ax1.set_xlim([0, len(fluweeklab)-1])
		plt.savefig('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/zip3Incidence/zip3Incidence_S%s_R%s.png' %(snum, rnum), transparent=False, bbox_inches='tight', pad_inches=0)
		plt.close()

		# 9/23/14 plot correlations for zip3 in single HHS region
		fig2 = plt.figure()
		ax2 = fig2.add_subplot(111)
		d_z1z2_corr = {}

		for z1, z2 in product(ls_zip3, ls_zip3):
			ls_z1 = d_zip3_incidTS[z1]
			ls_z2 = d_zip3_incidTS[z2]
			d_z1z2_corr[(z1,z2)] = stats.pearsonr(ls_z1, ls_z2)[0]

		ax2.plot([d_z1z2_corr[key] for key in d_z1z2_corr], 'o')
		ax2.set_xlabel('unique pair of zip3s')
		ax2.set_ylabel('pearson r for incidence TS')
		plt.savefig('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/zip3PearsonR/zip3PearsonR_S%s_R%s.png' %(snum, rnum), transparent=False, bbox_inches='tight', pad_inches=0)
		plt.close()
