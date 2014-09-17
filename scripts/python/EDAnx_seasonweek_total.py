#!/usr/bin/python

##############################################
###Python template
###Author: Elizabeth Lee
###Date: 9/12/14
###Function: Plot degree distributions for season week total plots

###Import data: 

###Command Line: python 
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
seasons = range(2,11)
kwargs_threshold = {'percent':90}

### import data ###
threshold_type, value = kwargs_threshold.items()[0]

for snum in seasons:
	base_edgelist = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/Py_export/edgelists/seasonweek_total_edgelist_pearson_%s%s_S%s.csv' %(threshold_type, value, snum)
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
	plt.xlim([0,375])
	plt.ylim([0,45])
	plt.title('Nodes: %s, Uq Edges: %s, Mean Degree: %s' %(nodes, edges, mean_degree))

	# save figure
	fname = '/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/plot_outputs/degreeDist_seasonweek_total_edgelist_pearson_%s%s_S%s.png' %(threshold_type, value, snum)
	plt.savefig(fname, transparent=False, bbox_inches='tight', pad_inches=0)
	plt.close()
