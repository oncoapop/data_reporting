#!/bin/sh

# Script to check coverage of genome (short)
# custom genome in fasta

# Directory
dir="/home/dyap/Projects/PrimerDesign/manual"

# custom genome
name="gp140" 
gen=$dir"/"$name".fa"

# Genome size
size=`cat $gen | grep ">" -A1 | grep -v ">" | wc -c`

	if [[ $size > "100000000" ]];
		then echo "Genome too large..."
		exit 1;
	fi

coverage=$size

# manifest QC file
QCfile=$dir"/"$name".AmpliconManifest_QC"

starange=`cat $QCfile | awk -F"\t" '{print $8}' | sort -u`

	nosta=`echo $starange | wc -l`
	for left in $starange
		do
		start=`echo $left | awk -F"=" '{print $2}' | awk -F"-" '{print $1}'`
		end=`echo $left | awk -F"=" '{print $2}' | awk -F"-" '{print $2}'`
		
		done
	
endrange=`cat $QCfile | awk -F"\t" '{print $9}' | sort -u`

