#!/bin/sh
# Script to convert lookup unknown barcodes from FLD list
# Version 1.0
# by Damian Yap, Molecular Oncology, BCCRC
# on beast.cluster.bccrc.ca
# 10 Jun 2013

# This is the path and name of the file that we want to identify the barcodes in
qpath="/home/dyap/Projects/MiSeq_Data/MiSeq_QC"
query=$qpath/"Matched_BC2.csv"

# Path and name of file that contains all the barcode mappings
inpath="/home/dyap/Projects/MiSeq_Data/MiSeq_QC"
source=$inpath/"FLD_barcodes.txt"

# Output file
out=$qpath/"Matched_Identified_BC2.csv"

echo "Barcode,Name,FLD identity,Count,," > $out

# Read the barcodes from one file and match them with source mapping file
for i in `cat $query  | awk -F, '{print $3}'`
	do
		FLD=`grep $i $source | awk -F"\t" '{print $1}'`
		BC=$i
		NAME=`grep $i $query  | awk -F, '{print $2}'`
		COUNT=`grep $i $query  | awk -F, '{print $6}'`
		echo $BC"," $NAME"," $FLD"," $COUNT",," >> $out
	done

exit 

