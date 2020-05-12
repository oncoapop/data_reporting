#!/bin/bash

# This Script was written by Damian Yap (Aug 2013)
# WSOP2013-001 version 4.0

# This script can inly be run AFTER the pipeline is run
# This just creates a bed file for viewing on UCSC

name="DG1136g"
dir="/home/dyap/Projects/PrimerDesign/TITAN-SS"

p3dir=$dir"/primer3"

cat $p3dir"/"$name"_primerlist.txt" | awk -F, '{print $1,$4,$5}' > $p3dir"/"$name"_isPCR-input"

# Module to check the primers by in silico PCR

suffix="_isPCR-input"


# isPCR is on beast at
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

# output format
#output=fa      # fasta format (default)
output=bed      # bed format (tab-delimited; Fields: chrom/start/end/name/score/strand)
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
