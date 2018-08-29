#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Parses UCSC Table Browser files to retrieve canonical exon CDS coordinates

#Modified to do CDS, not exons
import pickle

#UCSC Table knownGene, track Gencode v24, assembly hg38
gene_file = 'ref/hg38_ucsc_genes.txt'
#UCSC Table knownCanonical, track Gencode v24, assembly hg38
canonical = 'ref/hg38_canonical_ensg.txt'
#Gathered all the ENST ids from TCGA MAF files previously
gene_info = pickle.load(open('ref/tcga_gene_info.p', 'rb'))

gene_set = set()
info_dict = dict()
#Each element is ( Gene symbol, ENSG ID, ENST ID, HGNC ID)
#Associate ENSG IDs with TCGA gene symbols
for g in gene_info:
	info_dict[g[1]] = g[0]
	gene_set.add(g[0])

gene_dict = dict()

canonical_file = open(canonical)
#Skip header line
canonical_file.readline()
for line in canonical_file:
	spl = line.split('\t')
	chromosome = spl[0]
	#Ignore scaffolds
	if '_' in chromosome:
		continue
	ensg = spl[5]
	ensg = ensg[:ensg.index('.')]
	gene = info_dict.get(ensg, '')
	if gene != '':
		g_id = spl[4] #UCSC Gene ID
		gene_dict[g_id] = gene

cds_dict = dict()
reader = open(gene_file)
#Skip header line
reader.readline()

def between(cds_coord, x_s, x_e):
	'''
Check if the provided CDS coord in UCSC gene file is between the given exon start and end coordinates. 
Need to check so you can truncate exons where CDS starts/ends
	'''
	x_s = int(x_s)
	x_e = int(x_e)
	cds_coord = int(cds_coord)
	if (cds_coord > x_s) and (cds_coord < x_e):
		return True
	else:
		return False

for line in reader:
	spl = line.split('\t')
	u_id = spl[0]
	gene = gene_dict.get(u_id, '')
	if gene != '':
		chrom = spl[1]
		strand = spl[2]
		x_start = spl[8][:-1]#Remove trailing comma
		x_end = spl[9][:-1]#Remove trailing comma
		if ',' in x_start:
			x_start = x_start.split(',')
			x_end = x_end.split(',')
		else:
			x_start = [x_start]
			x_end = [x_end]

		cds_s = spl[5]
		cds_e = spl[6]

		cds_s_l = list()
		cds_e_l = list()
		cds_list = cds_dict.get(gene, list())
		exon = 0 #Keep track of which exon is which for downstream analysis outputs
		if cds_s == cds_e: #No CDS
			cds_list.append((chrom, [cds_s], [cds_e], gene, '0', strand, u_id, [-1], 0))
		else:
			n_exons = len(x_start)
			exon_num_list = list()
			for s, e in zip(x_start, x_end):
				if (int(cds_s) > int(e)) or (int(cds_e) < int(s)): #Exon is all non-coding
					exon += 1
					continue #Skip exon since it is entirely a UTR
				if between(cds_s, s, e):
					s = cds_s #If CDS starts in this exon, truncate exon coordinates
				if between(cds_e, s, e):
					e = cds_e #IF CDS ends in this exon, truncate exon coordinates
				cds_s_l.append(s)
				cds_e_l.append(e)
				exon += 1
				exon_out = exon
				if strand == '-': #Flip exon count for - strand
					exon_out = n_exons - exon + 1
				exon_num_list.append(exon_out)
			#Output similar to a BED file, but with added UCSC ID and exon number associated with information)
			cds_list.append((chrom, cds_s_l, cds_e_l, gene, '0', strand, u_id, exon_num_list, n_exons))
		cds_dict[gene] = cds_list
		
pickle.dump(cds_dict, open('ref/tcga_cds_dict.p', 'wb'))
