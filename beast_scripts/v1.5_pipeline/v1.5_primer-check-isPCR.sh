#!/bin/sh
# Script to check the primers by in silico PCR

# IF run as pipeline (comment this out)
# Project="TNBC"
# type="SNV"
# name=$Project"-"$type

suffix="_isPCR-input"
p3dir=$dir"/primer3"

outfilesuffix="_isPCR-output.txt"

# isPCR is on beast at 
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

# output format
output=fa 	# fasta format (default)
#output=bed 	# bed format (tab-delimited; Fields: chrom/start/end/name/score/strand)
#output=psl	# blat format 

# Name of the input file
inputfile=$p3dir"/"$name$suffix

# Name of the output file
outputfile=$p3dir"/"$name$outfilesuffix

cat $inputfile
echo $outputfile
$command $database $flip "-out="$output $inputfile $outputfile




