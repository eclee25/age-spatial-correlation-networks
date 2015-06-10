#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 10/2/14
###Function: Plot degree distributions for season week age-specific plots

###Import data: 

###Command Line: python EDAnx_seasonweek_age.py
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

### import data ###
for snum in seasons:
	# import age-specific edgelist for a given season
	for ac in agecodes:
		base_edgelist = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/Py_export/edgelists/seasonweek_age_edgelist_%s%s_%s%s_S%s%s.csv' %(TSdata, method, threshold_type, value, snum, ac)
		graphfile = open(base_edgelist,'r')
		G = ng.import_network(graphfile)

### program ###
		nodes, edges = G.number_of_nodes(), G.number_of_edges()
		print "network size", nodes
		mean_degree = np.mean(edges/nodes)
		print 'mean deg', mean_degree

		# plot degree distribution
		degrees = G.degree()
		values = sorted(set(degrees.values()))
		hist = [degrees.values().count(x) for x in values]
		plt.plot(values, hist, 'ko-')
		plt.axvline(x=mean_degree, color = 'r')
		plt.xlabel('number of correlated zip3s')
		plt.ylabel('zip3 count')
		plt.xlim([0,100])
		# plt.ylim([0,45])
		plt.title('Nodes: %s, Uq Edges: %s, Mean Degree: %s, Age: %s' %(nodes, edges, mean_degree, ac))

		# save figure
		fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/degreeDist_seasonweek_age_edgelist_%s%s_%s%s_S%s%s.png' %(TSdata, method, threshold_type, value, snum, ac)
		plt.savefig(fname, transparent=False, bbox_inches='tight', pad_inches=0)
		plt.close()
