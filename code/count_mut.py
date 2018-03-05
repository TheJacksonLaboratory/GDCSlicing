#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Counts the number of patients with mutations in each gene for each tumor and outputs into a table
#Note: This was written for Python 3
import csv
import glob
import pandas as pd
import pickle

maf_files = glob.iglob('*_filtered.csv') #Use filtered files, otherwise will fail
gene_list = pickle.load(open('gene_list.p', 'rb'))
count_dict = dict()
size_dict = dict()

for f in maf_files:
	tumor = f[:f.index('_')]
	reader = open(f)
	line = reader.readline()
	#This is why the number of samples is preserved in filtered CSV
	if 'samples' in line:
		sample_num = int(line.split(' ')[1])
		size_dict[tumor] = sample_num
			
	#Should be past comment lines at this point
	csv_reader = csv.DictReader(reader)
	
	#Only really care about a couple of columns: Hugo_Symbol, maybe Entrez_Gene_Id, Variant_Classification, Variant_Type, Tumor_Sample_Barcode
	#Only count 1 mutation per patient, filter out silent mutation (Variant_Classification = 'Silent')
	#This could easily change if needed
	mut_dict = dict()
	
	for line in csv_reader:
		var_class = line['Variant_Classification']
		#Ignore silent mutations for now, should have already been filtered out previously
		if var_class == 'Silent':
			continue
		gene = line['Hugo_Symbol']
		var_type = line['Variant_Type']
		barcode = line['Tumor_Sample_Barcode']
		
		patients = mut_dict.get(gene, set())
		#Only consider 1 mutation per patient
		if barcode not in patients:
			patients.add(barcode)
		mut_dict[gene] = patients
	
	#This will contain the patients with a mutation in each gene, organized by tumor type
	count_dict[tumor] = mut_dict

tumors = sorted(count_dict.keys())
gene_set = set()
sizes = pd.Series(size_dict)

#Where we actually tabulate the number of patients with mutations per gene per tumor
df = pd.DataFrame(0, index=gene_list, columns=tumors, dtype='int')

for t in tumors:
	mut_dict = count_dict[t]
	for g in gene_list:
		patients = mut_dict.get(g, set())
		df[t][g] = len(patients)
		
pickle.dump(df, open('mut_table.p', 'wb'))
pickle.dump(df/sizes, open('mut_frac.p', 'wb'))
pickle.dump(sizes, open('sample_sizes.p', 'wb'))
