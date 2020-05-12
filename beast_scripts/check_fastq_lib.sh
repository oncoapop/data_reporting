#!/bin/sh

# Script that will check the fastq files for PCR primers 

# reads lib location from
library="/home/dyap/Projects/Takeda_T3/primerQC/lib"

for id in `cat $library | awk -F"\t" '{print $1}'`
	do
	lib=`grep $id $library | awk -F"\t" '{print $2}'`	
        path="/share/lustre/archive/"$id"/illumina_wtss/"$lib"/sequence"

	cd $path
	for file in `ls *2.fastq` 
		do
		echo "File query= "$file
		done

	done
