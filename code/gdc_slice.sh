#!/bin/bash
#Author: Victor Wang
#Affiliation: The Jackson Laboratory for Genomic Medicine
#Affiliation: The University of Connecticut School of Medicine, Department of Genetics and Genomic Medicine
#Desc: Bash script for using GDC's BAM slicing API. A PBS Torque version is also available.

##########################
# USER SPECIFIED OPTIONS #
##########################

TEMP_LOC=temp
#CSV formatted file with GENE in col 1, TUMORS in following columns. See Supp Table 2 for proper format.
ANALYZE_FILE=
#GDC Token for downloading controlled data
TOKENFILE=

#Tests to determine if environment is set up properly

if [ -z "$ANALYZE_FILE" ]
then
	>&2 echo "No file was provided to determine genes/tumors to download"
	exit 1
fi

if [ -z "$TOKENFILE" ]
then
	>&2 echo "No token was provided to download controlled data from GDC"
	exit 2
fi

#Text files which contain filenames to download for a given cohort
if [ ! -d "filenames/" ]
then
	>&2 echo "No filenames directory found, populate filenames directory"
	exit 3
fi

if [ ! -d "bed/" ]
then
	>&2 echo "No exon coordinate directory found, provide gene-by-gene bed files"
	exit 4
fi

if [ -d "$( whereis samtools )" ]
then 
	>&2 echo "samtools is required for this script to work properly. Please download samtools"
	exit 5
fi

if [ -d "$( whereis python )" ]
then
	>&2 echo "Python is required for this script to work properly."
	exit 6
fi

#End environment test

if [ ! -d "missing/" ]
then
	mkdir missing
fi

if [ ! -d "results/" ]
	mkdir results
fi

#USER MUST PROVIDE PATH TO GDC TOKEN as $TOKENFILE
TOKEN=$( < $TOKENFILE )

#USER MUST PROVIDE PATH TO ANALYSIS FILE WHERE FIRST COLUMN IS GENE
#AND THE FOLLOWING COLUMNS ARE WHICH TUMORS TO DOWNLOAD FROM
LINE=$( sed -n "${PBS_ARRAYID}p" $ANALYZE_FILE )
GENE=$( echo $LINE | cut -f 1 -d , ) #Gene in first column
TUMORS=$( echo $LINE | cut -f 2- -d , ) #Tumors to look at in following columns
TUMORS=${TUMORS//,/ }

###################################################
#ALTERNATIVELY YOU CAN SIMPLY PROVIDE A FILE WITH COHORTS
#TO DOWNLOAD FROM LINE BY LINE (replace all_tumors.txt) AND
#ANALYZE_FILE CAN CONTAIN A LIST OF GENES
###################################################
#LINE=$( sed -n "${PBS_ARRAYID}p" $ANALYZE_FILE )
#TUMORS=$(<all_tumors.txt)

#Build out one long URL and series of outputs to maintain a single connection to
#GDC instead of opening up multiple connections which is more prone to failure
URL=""
OUT_STR=""
for TUMOR in $TUMORS #Iterate through all tumors
do
	NAME_FILE=filenames/${TUMOR}_filenames.txt #File contains all the bam files per tumor to look at
	TEMP_DIR=$TEMP_LOC/${TUMOR}_${GENE}
	mkdir -p $TEMP_DIR
	#BASE_URL will be used to insert the appropriate file
	BASE_URL=https://api.gdc.cancer.gov/slicing/view/FOO?gencode=$GENE
	for FILENAME in $( < $NAME_FILE ) #Go through all BAM files per tumor
	do
		URL="$URL ${BASE_URL/FOO/$FILENAME}"
		OUT_STR="$OUT_STR -o $TEMP_DIR/${FILENAME}_slice.bam"
	done
done

curl -s -S -v --header "X-Auth-Token:$TOKEN" $OUT_STR $URL

#If any files failed to download for whatever reason, want to be able to rerun
MISSING=missing/${GENE}_failed.txt
echo -n "" > $MISSING

for TUMOR in $TUMORS
do
	NAME_FILE=filenames/${TUMOR}_filenames.txt #File contains all the bam files per tumor to look at
	OUT_FILE=results/${GENE}_${TUMOR}_stats.txt
	echo -n "" > $OUT_FILE #Clear file in case it's already there from a previous run (for some reason)
	TEMP_DEPTH=$TEMP_LOC/${GENE}_${TUMOR}_depth.txt
	TEMP_DIR=$TEMP_LOC/${TUMOR}_${GENE}
	
	for FILENAME in $( < $NAME_FILE )
	do
		TEMP_BAM=$TEMP_DIR/${FILENAME}_slice.bam
		#Check if downloaded BAM is intact
		if ! samtools quickcheck $TEMP_BAM 
		then
			echo "$TUMOR,$FILENAME" >> $MISSING
			continue
		fi
		FILTER_BAM=$TEMP_BAM.rmdup.hq
		#Filter out duplicate and low mapping quality reads
	 	samtools rmdup $TEMP_BAM - | samtools view -hq 30 - > $FILTER_BAM
		#Count reads for a sanity check
		COUNT=$( samtools view -c $FILTER_BAM )
		if [ "0" -eq "$COUNT" ]
		then
			DEPTH=0	
		else
			#Calculate depth of exons
			samtools depth -a -b bed/${GENE}.bed $FILTER_BAM > $TEMP_DEPTH
			DEPTH=$( python avg_depth.py $TEMP_DEPTH )
		fi
		echo "$FILENAME ; $COUNT ; $DEPTH" >> $OUT_FILE
	done
	if [ ! -s "$MISSING" ]
	then
		rm $MISSING
	fi
	rm $TEMP_DEPTH
	rm -r $TEMP_DIR
done
