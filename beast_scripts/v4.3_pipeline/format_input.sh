#!/bin/sh

# Script to format the input file 

dir="/home/dyap/Projects/PrimerDesign/cellmix3"
cd $dir

outdir=$dir"/positions"
mkdir $outdir

cat targets_54.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_54.tsv"
cat targets_55.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_55.tsv"
cat targets_Shared.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_Shared.tsv"



