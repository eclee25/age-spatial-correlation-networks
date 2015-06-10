#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/2/14
###Function: Export networks as graphml format for import into gephi. Node attributes lat and lng should be included.

###Import data: 

###Command Line: python graphml_seasonweek_age.py
##############################################


### notes ###


### packages/modules ###
import csv
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

## local modules ##
import network_generation as ng
### data structures ###


### parameters ###
seasons = ng.sp_seasons
kwargs_threshold = ng.cp_threshold_kwargs
kwargs_TSdata_method = ng.cp_TSdata_method_kwargs
threshold_type, value = kwargs_threshold.items()[0]
TSdata, method = kwargs_TSdata_method.values()
agecodes = ng.sp_agegroupcodes


### import lat/long coordinates ###
# format: zip3, lat, long
llfile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/outside_source_data/zip3_latlong.csv'
llfilein = open(llfile, 'r')
d_zip3_lat, d_zip3_lng = ng.import_zip3_latLng(llfilein)

### import list of non-continental zip3s to exclude from visualization ###
noncontzipsin = open('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/R_export/non_continental_zip3s.csv', 'r')
excl_zips = ng.import_exclude_zip3s(noncontzipsin)

### import edgelist ###
for snum in seasons:
	# import age-specific edgelist for a given season
	for ac in agecodes:
		base_edgelist = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/Py_export/edgelists/seasonweek_age_edgelist_%s%s_%s%s_S%s%s.csv' %(TSdata, method, threshold_type, value, snum, ac)
		graphfile = open(base_edgelist, 'r')
		G = ng.import_network_continentalUS(graphfile, excl_zips)
		# remove zip3 from lat/long dicts if not in graph
		d_zip3f_lat, d_zip3f_lng = ng.filter_latlong_dicts(d_zip3_lat, d_zip3_lng, G)
		# set latitude node attribute
		nx.set_node_attributes(G, 'lat', d_zip3f_lat)
		# set longitude node attribute
		nx.set_node_attributes(G, 'lng', d_zip3f_lng)
		# write to file
		savefile = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/gephi/input_files/seasonweek_age_edgelist_%s%s_%s%s_S%s%s.graphml' %(TSdata, method, threshold_type, value, snum, ac)
		nx.write_graphml(G, savefile)