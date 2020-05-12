#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# Modified to generate Multiplex design primers

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/SNV"

clear 
cd $dir
ls

# Source and Output directories where Barcoded files stored
sourcedir=$dir
outdir=$dir

# Part of the pipeline, use default output of fasta2primer3.sh
filename=primer3_input.txt

rm -f nohup.out 

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here 
# name=SA029

outfile=$dir/$name"_p3_output"
infile=$dir/$filename

echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

primer3_core -output=$outfile -p3_settings_file=/home/dyap/Scripts/primer3_settings.txt < $infile

echo "Press enter to continue..."
read ans

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


echo "Press Enter to show..."

read ans
grep -B24 "PAIR_NUM_RETURNED=0" $outfile | more

echo Press Return to continue...
read ans

view=$outfile"_view.txt"

primer3_core -output=$view -format_output -p3_settings_file=/home/dyap/Scripts/primer3_settings.txt < $infile

echo Formatting for viewing... 
echo Press return to continue...

read ans

mv $outfile $outfile".txt"
