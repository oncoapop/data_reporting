#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Working Directory
dir="/home/dyap/Projects/Tumour_Evol/"

clear 
cd $dir
ls

# Source and Output directories where Barcoded files stored
sourcedir=$dir

echo "Input the name of the file that you want to run (without .ext but must be .txt):"
read filename

rm -f nohup.out 
outdir=Primer3_outputs
outfile=$outdir/$filename"_primer3"
infile=$dir$filename".txt"

echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir


nohup primer3_core < $infile

read ans

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c nohup.out

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" nohup.out 


echo "Press Enter to show..."

read ans
grep -B24 "PAIR_NUM_RETURNED=0" nohup.out | more

echo Press Return to continue...
read ans

nohup primer3_core -format_output < $infile

echo Formatting for viewing... 
echo Press return to continue...

read ans
view=$outfile"_view.txt"
cp nohup.out $view


