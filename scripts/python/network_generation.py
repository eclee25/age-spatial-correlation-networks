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
sp_agegroupcodes = ['T', 'C', 'A', 'E']

# keywords (function: crosscorrelation_edgelist)
sp_threshold_kwargs1 = {'percent':90}
sp_threshold_kwargs2 = {'threshold':0.8}
sp_threshold_kwargs3 = {'percent':99}
# keywords (function: timeseries_crosscorrelations)
sp_TSdata_method_kwargs1 = {'TS data': 'raw', 'method':'Pearson'}
sp_TSdata_method_kwargs2 = {'TS data': 'fd', 'method':'Pearson'} # fd = first differenced
sp_TSdata_method_kwargs3 = {'TS data': 'zn', 'method':'Pearson'} # zn = z-normalized


##############################################
### call parameters ### 
# calls preset parameters or otherwise sets frequently changing parameters

# keywords (function: crosscorrelation_edgelist)
cp_threshold_kwargs = sp_threshold_kwargs3
# keywords (function: timeseries_crosscorrelations)
cp_TSdata_method_kwargs = sp_TSdata_method_kwargs2


##############################################
### import functions ###
##############################################
def import_yearZip3_popAge_data(infile):
	''' Import allyears_zip3_popAge.csv for small age-specific population by calendar year by zip3. Returns one dict: dict_yearZip3Smallage_pop[(year, zip3, smallage)] = population size
	'''
	main(import_yearZip3_popAge_data)
	# init dict
	dict_yearZip3Smallage_pop = {}
	# data import
	infile.readline() # rm header
	csvfile = csv.reader(infile, delimiter=',')
	for row in csvfile:
		year, zip3, smallage = int(row[0]), str(row[1]), str(row[2])
		pop = float(row[3])
		dict_yearZip3Smallage_pop[(year, zip3, smallage)] = pop

	return dict_yearZip3Smallage_pop

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
def import_seasonweekAge_data(infile):
	''' Import fluseasonweeks_zip3_incidAge_S#.csv data for spatial correlation networks of weekly age-specific ILI incidence time series of the flu season. Small age groups and total population are pulled from the dataset. Returns two dicts: dict_week_weeknum[wk] = weeknum; dict_weekZip3_Smallage[(wk, zip3, small age)] = ili
	'''
	main(import_seasonweekAge_data)
	# init dicts
	dict_week_weeknum = {}
	dict_weekZip3Smallage_ili = {}
	# data import
	infile.readline() # rm header from mysql
	csvfile = csv.reader(infile, delimiter=',')
	for row in csvfile:
		week, weeknum = row[1], int(row[2])
		zip3, smallage, ili = str(row[3]), str(row[4]), float(row[5])
		wk = date(int(week[:4]), int(week[5:7]), int(week[8:]))
		# assign data to dicts
		dict_week_weeknum[wk] = weeknum
		dict_weekZip3Smallage_ili[(wk, zip3, smallage)] = ili

	return dict_week_weeknum, dict_weekZip3Smallage_ili

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
def bin_smallage_age_ili(dict_weekZip3Smallage_ili):
	''' Bin ILI data from small age group level to toddler (<=4), children (5-19), adult (20-59), elderly (60+) age groups. Returns one dict dict_weekZip3Age_ili[(wk, zip3, age)] = ili
	'''
	main(bin_smallage_age_ili)
	# create dict_smallAge_age['agecode'] = [list of small age strings]
	dict_smallAge_age = createDict_smallage_age()
	# generate list of unique weeks
	unique_weeks = list(set([key[0] for key in dict_weekZip3Smallage_ili]))
	# generate list of unique zip3s
	unique_zip3s = list(set([key[1] for key in dict_weekZip3Smallage_ili]))
	# initialize binned dict: dict_weekZip3Age_ili[(wk, zip3, age)] = ili
	dict_weekZip3Age_ili = {}
	# iterate through ILI dict, binning smallAge to age
	for wk, zip3 in product(unique_weeks, unique_zip3s):
		# bin toddlers
		T = [dict_weekZip3Smallage_ili.get((wk, zip3, smallage), 0) for smallage in dict_smallAge_age['T']]
		# bin children
		C = [dict_weekZip3Smallage_ili.get((wk, zip3, smallage), 0) for smallage in dict_smallAge_age['C']]	
		# bin adults
		A = [dict_weekZip3Smallage_ili.get((wk, zip3, smallage), 0) for smallage in dict_smallAge_age['A']]	
		# bin elderly
		E = [dict_weekZip3Smallage_ili.get((wk, zip3, smallage), 0) for smallage in dict_smallAge_age['E']]
		# add to dict only if all ILIs are non-zero
		if sum(T)>0 and sum(C)>0 and sum(A)>0 and sum(E)>0:
			dict_weekZip3Age_ili[(wk, zip3, 'T')] = sum(T)
			dict_weekZip3Age_ili[(wk, zip3, 'C')] = sum(C)
			dict_weekZip3Age_ili[(wk, zip3, 'A')] = sum(A)
			dict_weekZip3Age_ili[(wk, zip3, 'E')] = sum(E)

	return dict_weekZip3Age_ili

##############################################
def bin_smallage_age_pop(dict_yearZip3Smallage_pop):
	''' Bin popstat data from small age group level to toddler (<=4), children (5-19), adult (20-59), elderly (60+) age groups. Returns one dict dict_yearZip3Age_pop[(yr, zip3, age)] = pop
	'''
	main(bin_smallage_age_pop)
	# create dict_smallAge_age['agecode'] = [list of small age strings]
	dict_smallAge_age = createDict_smallage_age()
	# initialize dict
	dict_yearZip3Age_pop = {}
	# generate list of unique years
	unique_years = list(set([key[0] for key in dict_yearZip3Smallage_pop]))
	# generate list of unique zip3s
	unique_zip3s = list(set([key[1] for key in dict_yearZip3Smallage_pop]))
	# iterate through pop dict for each year-zip3 combo
	for yr, zip3 in product(unique_years, unique_zip3s):
		# bin toddler
		T = [dict_yearZip3Smallage_pop.get((yr, zip3, smallage), 0) for smallage in dict_smallAge_age['T']]
		# bin children
		C = [dict_yearZip3Smallage_pop.get((yr, zip3, smallage), 0) for smallage in dict_smallAge_age['C']]
		# bin adult
		A = [dict_yearZip3Smallage_pop.get((yr, zip3, smallage), 0) for smallage in dict_smallAge_age['A']]
		# bin elderly
		E = [dict_yearZip3Smallage_pop.get((yr, zip3, smallage), 0) for smallage in dict_smallAge_age['E']]
		# add to dict only if all populations are non-zero
		if sum(T)>0 and sum(C)>0 and sum(A)>0 and sum(E)>0:
			dict_yearZip3Age_pop[(yr, zip3, 'T')] = sum(T)
			dict_yearZip3Age_pop[(yr, zip3, 'C')] = sum(C)
			dict_yearZip3Age_pop[(yr, zip3, 'A')] = sum(A)
			dict_yearZip3Age_pop[(yr, zip3, 'E')] = sum(E)

	return dict_yearZip3Age_pop

##############################################
def prep_age_seasonweek_timeseries(dict_weekZip3Age_ili, dict_yearZip3Age_pop, agegroupCode, season):
	''' Prepare age-specific ILI by wk-zip3 and age-specific pop by yr-zip3 for use in seasonweek_timeseries function. This function subsets input dictionaries by input age group code. Returns two dicts: dict_weekZip3_iliVisit[(wk, zip3)] = (ili, dummyvalue), dict_zip3_pop[zip3] = pop
	'''
	main(prep_age_seasonweek_timeseries)
	# initialize dicts
	dict_weekZip3_iliVisit, dict_zip3_pop = {}, {}
	# subset pop dict based on equality between year suffix and season number, and agegroupcode
	dict_zip3_pop = dict((key[1], dict_yearZip3Age_pop[key]) for key in dict_yearZip3Age_pop if int(str(key[0])[-2:])==season and key[2]==agegroupCode)
	# subset ILI dict based on agegroupcode, 10/2/14 add dummy 0 value for "# visits"
	# rm from ILI dict if pop not in dict_zip3_pop
	dict_weekZip3_iliVisit = dict(((key[0], key[1]), (dict_weekZip3Age_ili[key], 0)) for key in dict_weekZip3Age_ili if key[2]==agegroupCode and key[1] in dict_zip3_pop)

	return dict_weekZip3_iliVisit, dict_zip3_pop

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
		# 10/2/14 calculate iliTS and incidTS only if popsize exists (see condition for dict_weekZip3_iliVisit in prep_age_seasonweek_timeseries)
		# NB. some zip3s missing data in a given week, 0 as placeholder
		iliTS = [dict_weekZip3_iliVisit.get((week,zip3),(0,0))[0] for week in list_week] 
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
def timeseries_crosscorrelations(dict_zip3_incidTS, **kwargs):
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
			TSdata, method = kwargs.get('TS data'), kwargs.get('method')
			# 9/29/14 kwargs control flow
			if TSdata == 'raw':
				list_zip_1 = dict_zip3_incidTS[zip_1]
				list_zip_2 = dict_zip3_incidTS[zip_2]
			elif TSdata == 'fd':
				list_zip_1 = np.diff(dict_zip3_incidTS[zip_1]) # first diff
				list_zip_2 = np.diff(dict_zip3_incidTS[zip_2])
			elif TSdata == 'zn':
				mn1, sd1 = np.mean(dict_zip3_incidTS[zip_1]), np.std(dict_zip3_incidTS[zip_1])
				mn2, sd2 = np.mean(dict_zip3_incidTS[zip_2]), np.std(dict_zip3_incidTS[zip_2])
				list_zip_1 = [(incid-mn1)/sd1 for incid in dict_zip3_incidTS[zip_1]] # z-normalized
				list_zip_2 = [(incid-mn2)/sd2 for incid in dict_zip3_incidTS[zip_2]]
			if method == 'Pearson':
				corr_coefficient = stats.pearsonr(list_zip_1, list_zip_2)[0]
			elif method == 'Sax':
				continue
			# assign correlation to dictionary 
			dict_zip3pair_crosscorr[(zip_1, zip_2)] = corr_coefficient
	print ct, 'zip3 pair correlations calculated'
	print 'with %s data and %s method' %(TSdata, method)

	return dict_zip3pair_crosscorr

##############################################
def crosscorrelation_edgelist(dict_zip3pair_crosscorr, **kwargs):
	''' Create zip3 network edgelist based on cross-correlations dictionary and a threshold. Returns one dict with edgelist as keys and cross-correlations as values: dict_zip3pair_filtered[(zip_1, zip_2)] = corr_coefficient; returns threshold value that defines edge existence between zip3s.
	'''
	main(crosscorrelation_edgelist)
	# list of all correlations
	all_corrs = [dict_zip3pair_crosscorr[key] for key in dict_zip3pair_crosscorr]
	# 9/29/14 kwargs control flow
	if kwargs.get('threshold'):
		threshold = kwargs.get('threshold')
	elif kwargs.get('percent'):
		threshold = np.percentile(all_corrs, kwargs.get('percent'))
	# correlation edgelist is filtered by threshold
	dict_zip3pairFilt_crosscorr = dict((key, dict_zip3pair_crosscorr[key]) for key in dict_zip3pair_crosscorr if dict_zip3pair_crosscorr[key] >= threshold)

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
## fixed calls ##
##############################################
def createDict_smallage_age():
	''' Initialize dictionary of toddlers, children, adults, and elderly to list of all small age groups in that age bin. Returns one dict:
	'''
	main(createDict_smallage_age)
	# initialize dict
	dict_smallAge_age = defaultdict(list)
	# define larger age groups with lists of small age group strings
	dict_smallAge_age['T'] = ['<2 years', '2-4 years']
	dict_smallAge_age['C'] = ['5-9 years', '10-14 years', '15-19 years']
	dict_smallAge_age['A'] = ['20-29 years', '30-39 years', '40-49 years', '50-59 years']
	dict_smallAge_age['E'] = ['60-69 years', '70-79 years', '80+ years']

	return dict_smallAge_age

##############################################



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
			fwriter.write("%s,%s,%s\n" %(key[0], key[1], value))


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