#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Determine read depth in canonical exon CDS regions based on samples from TCGA-UCS and output any problematic gene/exons

import glob
import pandas as pd
import os.path
import pickle
import numpy as np

cds_dict = pickle.load(open('tcga_cds_dict.p', 'rb'))
output_file = 'UCS_cds_undercovered.txt'

#These files must be requested due to file size constraints
g = glob.iglob('UCS_all/depth/*_base_depth.txt')

def cds_index(cds):
	'''
	Generate a Pandas Index based on the input CDS coordinates
	'''
	cds_s_l = cds[1]
	cds_e_l = cds[2]
	if cds_s_l[0] == cds_e_l[0]:
		return pd.Index([])
	np_index = np.array([], dtype=int)
	for s, e in zip(cds_s_l, cds_e_l):
		arr = np.arange(int(s) + 1, int(e) + 1)
		np_index = np.concatenate((np_index, arr))
	index = pd.Index(np_index)
	index = index.drop_duplicates()
	return index
	
frac_dict = dict()

all_exons_bad = 0
for f in g:
	df = pd.read_csv(f, index_col=0)
	bname = os.path.basename(f)
	gene = bname.replace('_UCS_base_depth.txt', '')
	cds = cds_dict[gene][0] #Only look at first one for right now
	cds_s_l = cds[1]
	cds_e_l = cds[2]
	if cds_s_l[0] == cds_e_l[0]: #No CDS
		continue
	for s, e, n in zip(cds_s_l, cds_e_l, cds[7]):
		index = pd.Index(np.arange(int(s) + 1, int(e) + 1))	
		df_p = df.loc[index]
		df_m = df_p.mean(axis=1) #Mean across all samples in UCS
		df_t = df_m > 20 #Threshold per base to determine good/bad read depth
		size = df_t.size #Exon length
		passed = np.sum(df_t) #Number of bases with sufficient depth
		frac = float(passed) / float(size)
		val = frac_dict.get(gene, [0, 0, list()]) 
		val[0] += passed
		val[1] += size
		val[2].append((n, frac))
		frac_dict[gene] = val

keys = list(frac_dict.keys())
keys.sort()

with open(output_file, 'w') as writer:
	writer.write('Gene Name\tFraction of Bases >20 reads\t# Exons\tExon #: Fraction of Bases >20 reads in exon\n')
	for k in keys:
		v = frac_dict[k]
		t_frac = float(v[0]) / float(v[1])
		if t_frac < 0.8: #Only write out info if gene overall has less than ideal coverage
			writer.write('{:s}\t{:0.3f}\t'.format(k, t_frac)) #Write gene name
			cds = cds_dict[k]
			exon_num = cds[0][8]
			writer.write(str(exon_num) + '\t')
			exons = v[2]
			exons.sort(key=lambda x:x[0])
			output = ''
			#Write each exon which also falls below the coverage threshold
			for exon in exons:
				x_frac = exon[1]
				if x_frac < 0.8:
					x_n = exon[0]
					output += '{:d}:{:0.3f},'.format(x_n, x_frac)
			output = output[:-1] + '\n'
			writer.write(output)
