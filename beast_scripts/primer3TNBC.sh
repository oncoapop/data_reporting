#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# Modified for TNBC Pipeline (v2.0)

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Working Directory
dir="/home/dyap/Projects/TNBC/primer3"

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

outfile=$dir/$name"_p3_output"
infile=$dir/$filename

echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

primer3_core -output=$outfile < $infile

echo "Number of sequences in input file" > $outfile.screen
grep "SEQUENCE_TEMPLATE=" -c $infile >> $outfile.screen

echo "Number of outputs in " $outfile >> $outfile.screen
grep "^=" -c $outfile >> $outfile.screen

echo ================================

echo "Number of failed sequences where there are no primers" >> $outfile.screen
grep -c "PAIR_NUM_RETURNED=0" $outfile >> $outfile.screen

grep -B30 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

###################################
# reiteration module goes in here #
###################################

# View file not required any more (takes too much time)
# view=$outfile"_view.txt"
# primer3_core -output=$view -format_output < $infile

mv $outfile $outfile".txt"
