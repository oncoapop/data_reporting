#!/bin/sh
# Script to check the primers by in silico PCR

# IF run as pipeline (comment this out)
Project="manual"
name="PPoI"

suffix="_isPCR-input"
dir="/home/dyap/Projects/PrimerDesign/"$Project"/"$name

outfilesuffix="_isPCR-output.txt"

# isPCR is on beast at 
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

# Max size of the amplicon that we allow in isPCR
size=600

# output format
output=fa 	# fasta format (default)
#output=bed 	# bed format (tab-delimited; Fields: chrom/start/end/name/score/strand)
#output=psl	# blat format 

# Name of the input file
inputfile=$dir"/"$name$suffix

# Name of the output file
outputfile=$dir"/"$name$outfilesuffix

cat $inputfile
echo $outputfile
$command "-maxSize="$size $database $flip "-out="$output $inputfile $outputfile




