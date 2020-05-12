#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# Modified to generate Multiplex design primers

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Working Directory
dir="/share/lustre/backup/dyap/Projects/PrimerDesign/manual"
#dir="/home/dyap/Projects/PrimerDesign/manual"

clear 
cd $dir
ls

# Source and Output directories where Barcoded files stored
sourcedir="/home/dyap/dyap_temp/Pipeline/pipelines/primer_design_pipeline_v1.0/SA535-2/SA535-2_primer_design_pipeline_beast/outputs/intermediate"
outdir=$dir

# Part of the pipeline, use default output of fasta2primer3.sh
#filename=primer3_input.txt

# use of a manual input file TASK_2 (from primer design pipeline)
filename="TASK_2_primer3Input.txt"

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here 
name=SA535-BRCA

outfile=$dir/$name"_p3_output"
settings="/home/dyap/Projects/PrimerDesign/settings/manual_p3_settings.txt"
infile=$sourcedir/$filename

echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

primer3_core -output=$outfile -p3_settings_file=$settings < $infile

echo "Press enter to continue..."
read ans

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "SEQUENCE_TEMPLATE=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


echo "Press Enter to show..."

read ans
grep -B24 "PAIR_NUM_RETURNED=0" $outfile | more 
grep -B24 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

echo Press Return to continue...
read ans

view=$outfile"_view.txt"

primer3_core -output=$view -format_output < $infile

echo Formatting for viewing... 
echo Press return to continue...

read ans

mv $outfile $outfile".txt"
