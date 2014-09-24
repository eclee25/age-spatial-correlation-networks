#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 9/11/14
###Function: functions to generate spatial correlation networks from spatially-oriented ILI data

###Import data: 

###Command Line: python 
##############################################


### notes ###


##############################################
### packages/modules ###
import csv
from collections import defaultdict
from itertools import product
from scipy import stats
from datetime import date
import numpy as np
import networkx as nx

##############################################
### (pre)set parameters ### 
sp_seasons = range(2,11)
sp_HHSregions = range(1,11)
sp_fluweeklabels = range(40,54) # week number labels for plots vs. time
sp_fluweeklabels.extend(range(1,21))

##############################################
### call parameters ### 
# calls preset parameters or otherwise sets frequently changing parameters
cp_threshold_kwargs = {'percent':90} # crosscorrelation_edgelist

##############################################
### import functions ###
##############################################
def import_seasonweek_data(infile):
	''' Import fluseasonweeks_zip3_incid_S#.csv data for spatial correlation networks of weekly ILI incidence time series of the flu season. Returns three dicts: dict_week_weeknum[wk] = weeknum; dict_weekZip3_iliVisit[(wk, zip3)] = (ili, visit); dict_zip3_pop[zip3] = pop
	'''
	main(import_seasonweek_data)
	# init dicts
	dict_week_weeknum = {}
	dict_weekZip3_iliVisit = {}
	dict_zip3_pop = {}
	# data import
	infile.readline() # rm header from mysql
	csvfile = csv.reader(infile, delimiter=',')
	for row in csvfile:
		week, weeknum = row[1], int(row[2])
		zip3, ili, visit, pop = str(row[3]), float(row[4]), int(row[5]), int(row[6])
		wk = date(int(week[:4]), int(week[5:7]), int(week[8:]))
		# assign data to dicts
		dict_week_weeknum[wk] = weeknum
		dict_weekZip3_iliVisit[(wk, zip3)] = (ili, visit)
		dict_zip3_pop[zip3] = pop

	return dict_week_weeknum, dict_weekZip3_iliVisit, dict_zip3_pop

##############################################
def import_network(infile):
	''' Import edgelist into graph object. Returns graph object.
	'''
	main(import_network)
	# create empty network
	Graph = nx.Graph()
	for edge in infile:
		edge_ls = edge.strip().split(',')
		Graph.add_edge(*edge_ls)

	return Graph

##############################################
def import_zip3_region_data(infile):
	''' Import crosswalk of zip3 and HHS region from explore/ project (R_export). Returns one dict: dict_zip3_HHS[zip3] = HHS region number.
	'''
	main(import_zip3_region_data)
	# init dict
	dict_zip3_HHS = {}
	# data import
	infile.readline() # rm header from R export
	csvfile = csv.reader(infile, delimiter=',')
	for row in csvfile:
		zip3, HHSreg = str(row[0]), int(row[8])
		# assign data to dict
		dict_zip3_HHS[zip3] = HHSreg

	return dict_zip3_HHS

##############################################
## processing functions ##
##############################################
def seasonweek_timeseries(dict_weekZip3_iliVisit, dict_zip3_pop):
	''' Process original data dicts to time series dicts. Returns two dicts: dict_zip3_iliTS[zip3] = iliTS; dict_zip3_incidTS[zip3] = incidTS
	'''
	main(seasonweek_timeseries)
	# init dicts
	dict_zip3_iliTS = defaultdict(list)
	dict_zip3_incidTS = defaultdict(list)
	# sorted unique week and zip3 lists
	list_week, list_zip3 = sortedlists_week_zip3(dict_weekZip3_iliVisit)
	# assign ili weekly time series as list to zip3
	for zip3 in list_zip3:
		iliTS = [dict_weekZip3_iliVisit.get((week,zip3),(0,0))[0] for week in list_week] # some zip3s missing data, 0 as placeholder
		incidTS = [ili/dict_zip3_pop[zip3] for ili in iliTS]
		# assign data to dict
		dict_zip3_iliTS[zip3] = iliTS
		dict_zip3_incidTS[zip3] = incidTS

	return dict_zip3_iliTS, dict_zip3_incidTS

##############################################
def sortedlists_week_zip3(dict_weekZip3_iliVisit):
	''' Produce sorted week and zip3 lists for use in other functions. Returns two lists: list_week; list_zip3.
	'''
	main(sortedlists_week_zip3)
	list_week = sorted(list(set([week for week,zip3 in dict_weekZip3_iliVisit])))
	list_zip3 = sorted(list(set([zip3 for week,zip3 in dict_weekZip3_iliVisit])))
	return list_week, list_zip3

##############################################
def timeseries_crosscorrelations(dict_zip3_incidTS):
	''' Calculate cross-correlations between each pair of zip3s for their flu season incidence time series. Pearson's r is the default method for calculating time series correlations. Returns one dict: dict_zip3pair_crosscorr[(zip_1, zip_2)] = corr_coefficient
	'''
	main(timeseries_crosscorrelations)
	# init dicts
	dict_zip3pair_crosscorr = {}
	# sorted unique zip3 list
	list_zip3 = sorted([zip3 for zip3 in dict_zip3_incidTS])
	

	# calculate pairwise correlations between time series
	ct=0
	for zip_1, zip_2 in product(list_zip3, list_zip3):
		# skip dict entry if swapped tuple pair is already in dict
		if (zip_2, zip_1) in dict_zip3pair_crosscorr:
			continue
			# print zip_1, zip_2, 'pair already in dict'
		# skip dict entry if zip3 is equal
		elif zip_1 == zip_2:
			continue
		else:
			ct+=1
			list_zip_1 = dict_zip3_incidTS[zip_1]
			list_zip_2 = dict_zip3_incidTS[zip_2]
			corr_coefficient = stats.pearsonr(list_zip_1, list_zip_2)
			dict_zip3pair_crosscorr[(zip_1, zip_2)] = corr_coefficient
	print ct, 'zip3 pair correlations calculated'

	return dict_zip3pair_crosscorr

##############################################
def crosscorrelation_edgelist(dict_zip3pair_crosscorr, **kwargs):
	''' Create zip3 network edgelist based on cross-correlations dictionary and a threshold. Returns one dict with edgelist as keys and cross-correlations as values: dict_zip3pair_filtered[(zip_1, zip_2)] = corr_coefficient; returns threshold value that defines edge existence between zip3s.
	'''
	main(crosscorrelation_edgelist)
	# list of all correlations
	all_corrs = [dict_zip3pair_crosscorr[key][0] for key in dict_zip3pair_crosscorr]
	threshold = kwargs.get('threshold')
	threshold = np.percentile(all_corrs, kwargs.get('percent'))
	# correlation edgelist is filtered by threshold
	dict_zip3pairFilt_crosscorr = dict((key, dict_zip3pair_crosscorr[key]) for key in dict_zip3pair_crosscorr if dict_zip3pair_crosscorr[key][0] >= threshold)

	print 'correlation threshold:', threshold
	print 'number of zip3 pairs:', len(dict_zip3pairFilt_crosscorr)

	return dict_zip3pairFilt_crosscorr, threshold

##############################################
def exclude_noncontinental_zip3s(dict_zip3_iliTS, dict_zip3_incidTS, exclude_list):
	''' Remove non-continental zip3s from dict_zip3_iliTS and dict_zip3_incidTS for cleaner jflowmap flows files.
	'''
	main(exclude_noncontinental_zip3s)
	# remove non-continental zip3s from ili dict
	dict_zip3_iliTS_continental = dict((key, dict_zip3_iliTS[key]) for key in dict_zip3_iliTS if key not in exclude_list)
	# remove non-continental zip3s from incid dict
	dict_zip3_incidTS_continental = dict((key, dict_zip3_incidTS[key]) for key in dict_zip3_incidTS if key not in exclude_list)

	print len(dict_zip3_incidTS_continental)

	return dict_zip3_iliTS_continental, dict_zip3_incidTS_continental


##############################################
## write functions ##
##############################################
def write_flowscsv_jflowmap(dict_zip3pairFilt_crosscorr, filename):
	''' Write 'flows.csv' file for jflowmap implementation, with columns 'origin', 'destination', and 'correlation'. The origin and destination zip3s are arbitrary; the correlations are not actually directed. Trial use of jflowmap.
	'''
	main(write_flowscsv_jflowmap)
	with open(filename, 'w+') as fwriter:
		fwriter.write("%s,%s,%s\n" %('Origin','Dest','Correlation'))
		for key, value in dict_zip3pairFilt_crosscorr.items():
			fwriter.write("%s,%s,%s\n" %(key[0], key[1], value[0]))


##############################################
def write_edgelist(dict_zip3pairFilt_crosscorr, filename):
	''' Write zip3 pairs as an edgelist file.
	'''
	main(write_edgelist)
	with open(filename, 'w+') as fwriter:
		for key in dict_zip3pairFilt_crosscorr:
			fwriter.write("%s,%s\n" %(key[0], key[1]))



##############################################
##############################################
# footer

def main(function):
	print 'Running', __name__, function.__name__

if __name__ == '__main__':
	print 'Executed from the command line'
	main()