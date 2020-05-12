#!/bin/sh
# Script to check the primers by in silico PCR

# IF run as pipeline (comment this out)
Project="MitoVar"
#name="DAH55_56_Shared"
name="mito"

dir="/home/dyap/Projects/Single_Cell/"$Project
#dir="/home/dyap/Projects/PrimerDesign/manual"
outfilesuffix="_isPCR-output.txt"
suffix="_isPCR-input"

# Adaptors to be removed for isPCR
# Illlumina Adaptors (5'->3')
fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

# Forward adaptor for Fluidigm
# fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
# ra="TACGGTAGCAGAGACTTGGTCT"

# Automated removal of adaptors
sed "s/$fa//g" $dir"/"$name"_primers.txt" > $dir"/temp1.tmp"
sed "s/$ra//g" $dir"/temp1.tmp"| tr "," "\t" > $dir"/"$name$suffix
 
# isPCR is on beast at 
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
#database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"
#database="/home/dyap/Projects/PrimerDesign/manual/gp140.2bit"
database="/home/dyap/Projects/Single_Cell/MitoVar/mitochondria_complete.2bit"

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

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
$command $database $flip "-out="$output $inputfile $outputfile

rm -f $dir"/*.tmp"



