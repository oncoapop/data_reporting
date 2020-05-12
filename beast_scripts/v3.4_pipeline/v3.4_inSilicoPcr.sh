#!/bin/bash

# This Script was written by Damian Yap (Aug 2013)
# WSOP2013-001 version 4.0
# upgraded to specifically take primer3 records one by one (Feb 2017)
# This does PCR on the genomic DNA as a control for RT-PCR

# Script to generate isPCR input from primer3 output
# $Project and $type is exported from previous script
# $Project (project name) is exported
# $name (sample name)is exported
# $dir (~/Projects/Pipeline)  is exported

############## Manual run ##################
Project="eIF4A3"
name="eIF4A3_NMD"
dir="/home/dyap/Projects/PrimerDesign/"$Project
############################################

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
primerlist=$p3dir"/"$name"_primers_genomic.txt"

echo "File from primer3 output to process: "$outfile
echo Output to this directory $outdir

if [ -f $primerlist ];
	then
	echo $primerlist" will be overwritten. Press Return to continue, Ctrl-C to exit."
	# read ans
	rm $primerlist
fi


# Module to check the primers by in silico PCR

suffix="_isPCR-input"
p3dir=$dir"/primer3"

# isPCR is on beast at
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
# Genomic database
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"
# Splice custom transcriptomic database
# database="/home/dyap/Projects/PrimerDesign/manual/splice1.2bit"
# Genocode transcripts of all genes in fasta v19
# database="/home/dyap/dyap_temp/genomes/gencode.v19.pc_transcripts.2bit"

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
outfilesuffix="_isPCR-genomic."$output

# Name of the input file
inputfile=$p3dir"/"$name$suffix

# Name of the output file
outputfile=$p3dir"/"$name$outfilesuffix

cat $inputfile
echo "+++++++++++++++"
echo $outputfile
echo "+++++++++++++++"

echo "Performing in-silico PCR using primers on hg19.... (This takes at least 7 min for whole genome)"

while :;do echo -n .;sleep 1;done &
$command $database $flip "-out="$output $inputfile $outputfile
kill $!; trap 'kill $!' SIGTERM

echo "In-silico PCR is completed."

exit;
