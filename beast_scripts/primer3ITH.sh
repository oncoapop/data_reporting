#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# Modified for ITH Pipeline (v2.0)

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Working Directory
dir="/home/dyap/Projects/ITH/primer3"

clear 
cd $dir
ls

# Source and Output directories where Barcoded files stored
sourcedir=$dir
outdir=$dir

# Part of the pipeline, use default output of fasta2primer3.sh
filename="primer3_input.txt"

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here 
# name=Indel

outfile=$dir"/"$name"_p3_output"
infile=$dir"/"$filename
settings=$dir"/ITH-primer3_settings.txt"
echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

primer3_core -output=$outfile -p3_settings_file=$settings  <  $infile

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile

grep -B26 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

cat $outfile > $outfile".txt"

###################################
# reiteration module goes in here #
###################################

# Part of the redesign pipeline
rename=$name"_p3_redesign.txt"

#cat $outfile.failed | sed 's/PRIMER_PRODUCT_SIZE_RANGE=260-280\ 250-300/PRIMER_PRODUCT_SIZE_RANGE=240-340\ 
#220-380\ 200-420/'| sed 's/PRIMER_MAX_SIZE=25/PRIMER_MAX_SIZE=32/' | sed 
#'s/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | grep -A19 "PRIMER_SEQUENCE_ID=" | sed 's/--/=/' > $rename

cat $outfile.failed | sed 's/PRIMER_MAX_SIZE=25/PRIMER_MAX_SIZE=32/' | sed 's/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | grep -A19 "PRIMER_SEQUENCE_ID=" | sed 's/--/=/' > $rename

echo "=" >> $rename

reoutfile=$dir/$name"_p3_output2"
primer3_core -output=$reoutfile -p3_settings_file=$settings  < $rename


clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $rename

echo "Number of outputs in " $reoutfile
grep "^=" -c $reoutfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $reoutfile

grep -B26 "PAIR_NUM_RETURNED=0" $reoutfile > $reoutfile.failed

######################################

cat $reoutfile >> $outfile".txt"


exit;

