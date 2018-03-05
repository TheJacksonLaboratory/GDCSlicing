#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Filters MAF files to remove silent mutations and other columns we really don't care about later
#Note: This was written for Python 3
import glob
import csv
import pickle

#These MAF files should be downloaded from GDC
maf_files = glob.iglob('maf/*.maf')
gene_set = set()

for f in maf_files:
	ind = f.index('.') + 1
	tumor = f[ind:f.index('.', ind)]
	reader = open(f)
	writer = open(tumor + '_filtered.csv', 'w')
	#Preserve information about number of samples for simplicity
	for l in range(5):
		line = reader.readline()
		if 'samples' in line:
			writer.write(line)
			
	#Should be past comment lines at this point
	csv_reader = csv.DictReader(reader, delimiter='\t')
	writer.write('Hugo_Symbol,Variant_Classification,Variant_Type,Tumor_Sample_Barcode\n')
	
	for line in csv_reader:
		var_class = line['Variant_Classification']
		#Ignore silent mutations for now
		if var_class == 'Silent':
			continue
		gene = line['Hugo_Symbol']
		var_type = line['Variant_Type']
		barcode = line['Tumor_Sample_Barcode']
		gene_set.add(gene)
		
		writer.write('%s,%s,%s,%s\n' % (gene, var_class, var_type, barcode))
	
	writer.close()
		
gene_list = sorted(list(gene_set))
pickle.dump(gene_list, open('gene_list.p', 'wb'))
