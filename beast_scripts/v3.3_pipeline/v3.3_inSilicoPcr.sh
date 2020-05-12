#!/bin/bash

# This Script was written by Damian Yap (Aug 2013)
# WSOP2013-001 version 4.0

# Script to generate isPCR input from primer3 output
# $Project and $type is exported from previous script
# $Project (project name) is exported
# $name (sample name)is exported
# $dir (~/Projects/Pipeline)  is exported

# Project Directory

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
outfile=$p3dir"/"$name"_primer3_output.txt"

# Final output
primerlist=$p3dir"/"$name"_primerlist.txt"

echo "File from primer3 output to process: "$outfile
echo Output to this directory $outdir

if [ -f $primerlist ];
	then
	echo $primerlist" will be overwritten. Press Return to continue, Ctrl-C to exit."
	# read ans
	rm $primerlist
fi


wkfile=$outfile

echo "Processing..."
for i in `grep "SEQUENCE_ID=" $wkfile`
	do
	echo $i
	n=`grep -A12 $i $wkfile | grep "PRIMER_PAIR_NUM_RETURNED=" | awk -F"=" '{print $2}'`
		if [[ "$n" =~ "ok 0" ]];
			then continue
		fi
	for j in 0 1 2 3 4
		do
	left=`grep -A140 $i $wkfile | grep -m1 "PRIMER_LEFT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
	right=`grep -A140 $i $wkfile | grep -m1 "PRIMER_RIGHT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
	size=`grep -A160 $i $wkfile | grep -m1 "PRIMER_PAIR_"$j"_PRODUCT_SIZE" | awk -F"=" '{print $2}'`
	id=`echo $i | awk -F"=" '{print $2}'`
	snv=`echo $id | awk -F"_" '{print $NF}'`
	chr=`echo $id | awk -F"_" '{print $(NF-1)}'`
		if [ -z "$right" -o -z "$left"  ];
			then continue
		fi


		echo $id","$chr","$snv","$left","$right","$size >> $primerlist
		done

	echo -ne '#\r'
	done

	echo -ne '\n'

cat $p3dir"/"$name"_primerlist.txt" | awk -F, '{print $1,$4,$5}' > $p3dir"/"$name"_isPCR-input"

# Module to check the primers by in silico PCR

suffix="_isPCR-input"

p3dir=$dir"/primer3"

# isPCR is on beast at
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
#database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"
#database="/home/dyap/Projects/PrimerDesign/manual/splice1.2bit"
		if [ -z "$database"  ];
			then 
				echo $database
				exit1;
		fi

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

# output format
output=fa      # fasta format (default)
#output=bed      # bed format (tab-delimited; Fields: chrom/start/end/name/score/strand)
#output=psl     # blat format
outfilesuffix="_isPCR-output."$output

# Name of the input file
inputfile=$p3dir"/"$name$suffix

# Name of the output file
outputfile=$p3dir"/"$name$outfilesuffix

cat $inputfile
echo $outputfile

echo "Performing in-silico PCR using primers on hg19.... (This takes at least 7 min)"

while :;do echo -n .;sleep 1;done &
$command $database $flip "-out="$output $inputfile $outputfile
kill $!; trap 'kill $!' SIGTERM

echo "In-silico PCR is completed."

exit;
