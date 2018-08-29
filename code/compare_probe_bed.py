#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Determine overlap between canonical exon CDS regions with a provided capture kit's probe target BED file

import pickle
import numpy as np
import pandas as pd

#Sorted probe target BED file. If file is not sorted, this script will not work properly.
design_file = ''
output_file = ''

if design_file == '':
	print('No comparison BED file was provided')
	exit()

if output_file == '':
	output_file = design_file[:design_file.rindex('.')] + '.undercovered.txt'

def num_sorter_key(value):
	'''
	Comparison function for proper sorting of exons analogous to IGV
	'''
	chrom = value[0]
	chrom_n = chrom.replace('chr', '')
	if chrom_n == 'X':
		return 23
	elif chrom_n == 'Y':
		return 24
	else:
		return(int(chrom_n))

def compare_chrom(chr_a, chr_b):
	'''
	Comparison function for chromosomes to use integers instead of strings
	'''
	val_a = chr_a.replace('chr', '')
	val_b = chr_b.replace('chr', '')
	if val_a == 'X':
		val_a = 23
	elif val_a == 'Y':
		val_a = 24
	else:
		val_a = int(val_a)
		
	if val_b == 'X':
		val_b = 23
	elif val_b == 'Y':
		val_b = 24
	else:
		val_b = int(val_b)
	if val_a == val_b:
		return 0
	elif val_a > val_b:
		return 1
	else: #val_a < val_b
		return -1

def overlap_s(exon_series, s_2, e_2):
	'''
	Check for overlap of two sequences of coordinates encoded as a Pandas Series
	'''
	if (exon_series.index[0] < e_2) and (exon_series.index[-1] > s_2):
		return True
	else:
		return False

#Association between TCGA gene name and CDS information
cds_dict = pickle.load(open('ref/tcga_cds_dict.p', 'rb'))
#Coordinate-sort exons
exon_list = list() #Exon list, with CDS only
for key, value in cds_dict.items():
	cds = value[0] #Only look at first one for now
	chrom = cds[0]
	if '_' in chrom:
		#In scaffold, skip as design files don't contain scaffolds
		continue
	#Lists for exon start/end coordinates
	x_starts = cds[1]
	x_ends = cds[2]
	if x_starts[0] == x_ends[0]:
		continue #No CDS
	for i in range(len(x_starts)): #UCSC coordinates for lower index are off by 1
		x_s = x_starts[i]
		x_s_i = int(x_s) + 1
		x_starts[i] = str(x_s_i)
	strand = cds[5]
	x_n = cds[7]

	for s, e, n in zip(x_starts, x_ends, x_n):
		#Pandas Series have coordinates as the index, associated value is a boolean to say whether or not probes cover that coordinate
		exon_series = pd.Series(False, index=np.arange(int(s), int(e) + 1)) #Need to include higher index from UCSC
		exon_list.append([chrom, exon_series, key, strand, n])
		
#Sort by coordinate followed by chromosome
exon_list.sort(key=lambda x:x[1].index[0])
exon_list.sort(key=num_sorter_key)

it = iter(exon_list)
current = next(it)
stop = False

#Simultaneously go through both sorted BED file and sorted exon CDS list
for line in open(design_file):
	#Don't need to skip first line
	spl = line.split('\t')
	chrom = spl[0]

	if '_' in chrom:
		#Skip scaffolds
		continue
	
	x_start = int(spl[1])
	x_end = int(spl[2])
	if compare_chrom(chrom, current[0]) == -1:
		#UCSC exon list is ahead of the probe
		continue
	
	while compare_chrom(chrom, current[0]) == 1:
		#Probe is ahead of the UCSC exon list
		current = next(it)
	
	while (chrom == current[0]) and ( current[1].index[0] < x_end ): 
		if overlap_s(current[1], x_start, x_end):
			x_range = np.arange(x_start, x_end)
			index = current[1].index
			overlap_ind = index.isin(x_range)
			current[1][overlap_ind] = True #The overlapped coordinates have been visited
			if x_end < current[1].index[-1]: #Probe is lagging behind, don't move to next exon yet
				break #Move to next probe

		#Otherwise, move to next coordinate in UCSC exon list
		try:
			current = next(it)
		except StopIteration:
			#Hit the end of the UCSC exon list before the end of the probe BED file
			stop = True
			break
		
	if stop:
		break

frac_dict = dict()
for exon in exon_list:
	covered_num = np.sum(exon[1]) #Number of bases in this exon with probe coverage
	len_exon = len(exon[1].index)
	gene = exon[2]
	x_n = exon[4] #Exon number
	val = frac_dict.get(gene, [0, 0, list()])
	#Store fraction of covered bases in this exon
	val[0] += covered_num
	val[1] += len_exon
	val[2].append((x_n, float(covered_num)/float(len_exon)))
	frac_dict[gene] = val

keys = sorted(list(frac_dict.keys()))

with open(output_file, 'w') as writer:
	writer.write('Gene Name\tFraction of bases with probe coverage\t# Exons\tExon #: Fraction of bases with probe coverage in exon\n')
	for k in keys:
		v = frac_dict[k]
		t_frac = float(v[0]) / float(v[1])
		if t_frac < 0.8:
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
