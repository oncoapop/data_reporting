#!/bin/sh

# This Script was written by Damian Yap (Aug 2013)

# Primer3 2.3.5 is installed on ma.sdf
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here
Project="Single_Cell"
sample="ITH-Case2"
name="20130930221649"
Miseq=150

echo "Name of project: "$name
echo "If this is correct, press ENTER or else Ctrl-C to exit..."

# Project Directory
dir="/meta/o/oncoapop/Projects/Pipeline/"$Project"/"$sample

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

tmp="/meta/o/oncoapop/temp"

if [ -e $tmp ];
        then
        echo "tmp directory is ok"
        else
        mkdir $tmp
        chmod og-wrx $tmp
fi

# Source and Output directories where Barcoded files stored
sourcedir=$posdir
outdir=$tmp
settings="/meta/o/oncoapop/Projects/Pipeline/settings/"$p3set"_pipeline_p3_settings.txt"
raw=$posdir"/"$name"-hg19.csv"
uniq=$posdir"/"$name"_uniqpos.txt"

# Part of the pipeline, use default output of R-script GetRanges_TOV_SNV.R
in1=$p3dir"/"$name"_p3_design.txt"
primerlist=$outdir"/"$name"_primerlist.txt"
outfile=$outdir"/"$name"_p3_output"
infile=$outdir"/"$name"_primer3_input.txt"


export dir=$dir
export tmp=$tmp
export Project=$Project
export type=$type
export raw=$raw
export name=$name
export reads=$Miseq
export sample=$sample

echo "####################################################"
echo "Generating TEST order file..."

~/Scripts/v2_pipeline/v2_isPCRweb_primer_summary.sh

echo "Files transferred to html path for viewing..."

exit;

