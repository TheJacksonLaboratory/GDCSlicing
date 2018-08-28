#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Parses UCSC Table Browser files to retrieve canonical exon coordinates and write to individual BED files

import pickle

#UCSC Table knownGene, track Gencode v24, assembly hg38
ucsc_gene = 'hg38_ucsc_genes.txt'
#UCSC Table knownCanonical, track Gencode v24, assembly hg38
canonical = 'hg38_canonical_ensg.txt'
#UCSC files are in Gencode 24, get metadata
gencode = 'gencode.v24.metadata.HGNC'
#Gathered all the ENST ids from TCGA MAF files previously
gene_info = pickle.load(open('tcga_gene_info.p', 'rb'))

enst_set = set()
#Each element is ( Gene symbol, ENSG ID, ENST ID, HGNC ID)
for g in gene_info:
	enst_set.add(g[2])

id_set = set()
canonical_file = open(canonical)
#Skip header line
canonical_file.readline()
for line in canonical_file:
	spl = line.strip().split('\t')
	chromosome = spl[0]
	#Ignore scaffolds
	if '_' in chromosome:
		continue
	g_id = spl[4] #UCSC gene ID
	id_set.add(g_id)

gene_file = open(ucsc_gene)
#Skip header line
gene_file.readline()
exon_dict = dict()
for line in gene_file:
	spl = line.strip().split('\t')
	g_id = spl[0] #UCSC gene ID
	if g_id not in id_set: #Check if this gene is the canonical transcript
		continue
	enst = spl[11]
	#TCGA doesn't have ID version number, so ignore it (not safe but not many other workarounds)
	if enst[:enst.rindex('.')] not in enst_set:
		continue
	chrom = spl[1]
	strand = spl[2]
	exon_start = spl[8][:-1] #Remove trailing comma
	exon_end = spl[9][:-1]
	#Associate gene coordinate information with ENST
	exon_dict[enst] = [chrom, exon_start, exon_end, strand]

#Associate UCSC gene coordinate info with the Gencode v24 name
for line in open(gencode):
	spl = line.strip().split('\t')
	enst = spl[0]
	d = exon_dict.get(enst, None)
	if d != None:
		chrom = d[0]
		exon_s = d[1]
		exon_e = d[2]
		strand = d[3]
		gene = spl[1] #Gene name from Gencode v24
		#Write to BED file
		writer = open('tcga_bed_v24/' + gene + '.bed','w')
		for s, e in zip(exon_s, exon_e):
			writer.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, s, e, gene, 0, strand))
		writer.close()
