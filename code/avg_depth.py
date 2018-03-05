#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Calculate average exonic coverage based on the output of samtools depth
#Note: This was written for Python 3
import sys

#Pass in temporary file containing depth at each exon position
in_file = sys.argv[-1]

counter = 0
total_coverage = 0
for line in open(in_file):
	coverage = line.split('\t')[2]
	counter += 1
	total_coverage += int(coverage)
avg_coverage = float(total_coverage) / float(counter)
print(avg_coverage)
