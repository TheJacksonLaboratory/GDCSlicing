#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Converts UCSC information on canonical exons to BED files for exon coverage calculations
#Note: This was written for Python 3

#These two files are from the UCSC Table Browser
exon_bed = 'known_genes.bed'
canonical = 'known_canonical.txt'

#This file contains the genes you are interested in created BED files for
genes = 'genes.txt'

gene_set = set([x.strip() for x in open(genes)])
id_set = set()
gene_dict = dict()

#Store UCSC IDs of canonical exons only
canonical_file = open(canonical)
canonical_file.readline() #Ignore header
for line in canonical_file:
	spl = line.split('\t')
	chromosome = spl[0]
	#Ignore scaffolds as TCGA BAM files won't align there
	if '_' in chromosome:
		break
	gene = spl[4].strip()
	if gene in gene_set:
		ucsc_id = spl[3]
		id_set.add(ucsc_id)
		gene_dict[ucsc_id] = gene

exon_dict = dict()

#Pull exon coordinates for canonical exons
for line in open(exon_bed):
	spl = line.split('\t')
	ucsc_id = spl[3]
	ucsc_id = ucsc_id[:ucsc_id.index('_')]
	#Ignore coordinates that are not from the canonical exons
	#Ignore coordinates that are not from genes we care about
	if ucsc_id in id_set:
		gene = gene_dict[ucsc_id]
		exon_list = exon_dict.get(gene, list())
		chrm = spl[0]
		exon_start = spl[1]
		exon_end = spl[2]
		exon_list.append((chrm, exon_start, exon_end))
		exon_dict[gene] = exon_list

#Write out BED files
for gene in exon_dict.keys():
	exon_list = exon_dict[gene]
	writer = open('bed/' + gene + '.bed','w')
	for exon in exon_list:
		writer.write('%s\t%s\t%s\n' % exon)
	writer.close()
