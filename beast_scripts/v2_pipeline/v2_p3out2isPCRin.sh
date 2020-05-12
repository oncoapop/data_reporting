#!/bin/sh

# This Script was written by Damian Yap (Aug 2013)
# WSOP2013-001 version 4.0

# Script to generate isPCR input from primer3 output
# $Project and $type is exported from previous script

# Project Directory
#Project="TOV"
#datestamp=`date "+%Y%m%d%N"`

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

cd $dir
ls

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here 
# type="SNV"
# name=$Project"-"$type

echo $name

# Source and Output directories where Barcoded files stored
sourcedir=$posdir
outdir=$p3dir

# Part of the pipeline, use default output of fasta2primer3.sh
outfile=$tmp"/"$name"_p3_output"
infile=$tmp"/"$filename

# Final output
primerlist=$tmp"/"$name"_primerlist.txt"

echo "File from primer3 output to process: "$outfile
echo Output to this directory $outdir

if [ -f $primerlist ];
	then 
	echo $primerlist" will be overwritten. Press Return to continue, Ctrl-C to exit."
	read ans
	rm $primerlist
fi


wkfile=$outfile

echo "Processing..."
for i in `grep "SEQUENCE_ID=" $wkfile`
	do
	n=`grep -A12 $i $wkfile | grep "PRIMER_PAIR_NUM_RETURNED=" | awk -F"=" '{print $2}'`
		if [[ $n =~ "0" ]]; 
			then continue
		fi
	for j in 0 1 2 3 4
		do
	left=`grep -A140 $i $wkfile | grep -m1 "PRIMER_LEFT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
	right=`grep -A140 $i $wkfile | grep -m1 "PRIMER_RIGHT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
	size=`grep -A160 $i $wkfile | grep -m1 "PRIMER_PAIR_"$j"_PRODUCT_SIZE" | awk -F"=" '{print $2}'`
	snv=`grep -A10 $i $wkfile | grep -m1 "P3_COMMENT=" | awk -F"=" '{print $2}'`
	id=`echo $i | awk -F"=" '{print $2}'`

		if [ -z "$right" -o -z "$left"  ]; 
			then continue
		fi


		echo $id,$snv,$left,$right,$size >> $primerlist
		done

	echo "."
	done

cp $tmp"/"$name"_primerlist.txt" $p3dir"/"$name"_primerlist.txt"
cat $p3dir"/"$name"_primerlist.txt" | awk -F, '{print $2,$3,$4}' > $p3dir"/"$name"_isPCR-input"

echo "Primerlist generated for in-silico PCR..."

exit;
