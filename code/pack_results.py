#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Packges the results folder into a csv for each tumor
#This gives a csv with filenames as the columns, genes as the rows, and depth as the values
#Note: This was written for Python 3
import pandas as pd
import glob

result_dir = 'results/'
files = glob.iglob(result_dir + '*_stats.txt')
gene_dict = dict()
#Build dictionary containing list of genes for each tumor
for f in files:
	ind_0 = f.index('/') + 1 
	ind_1 = f.index('_', ind_0)
	gene = f[ind_0:ind_1] 
	ind_1 = ind_1 + 1
	ind_2 = f.index('_', ind_1)
	tumor = f[ind_1:ind_2]
	gene_list = gene_dict.get(tumor, list())
	gene_list.append(gene)
	gene_dict[tumor] = gene_list

tumors = gene_dict.keys()
for t in tumors:
	gene_list = gene_dict[t]
	depth_list = list()
	for g in gene_list:
		filename = '%s%s_%s_stats.txt' % (result_dir, g, t)
		samples = list()
		coverage_list = list()
		#Build up Pandas series to eventually export to CSV
		for line in open(filename):
			spl = line.split(';')
			sample = spl[0].strip()
			depth = spl[2].strip()
			samples.append(sample)
			coverage_list.append(depth)
		s_depth = pd.Series(coverage_list, name=g, index=samples)
		depth_list.append(s_depth)
	df = pd.concat(depth_list, axis=1).T
	df.to_csv('depth_results/%s_depth.csv' % t)
