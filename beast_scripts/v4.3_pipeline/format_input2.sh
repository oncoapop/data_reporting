#!/bin/sh

# Script to format the input file 

dir="/home/dyap/Projects/PrimerDesign/cellmix4"
cd $dir

outdir=$dir"/positions"
mkdir $outdir

cat targets_HCT.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_HCT.tsv"
cat targets_hTERT.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_hTERT.tsv"
cat targets_Shared_HCT_hTERT.tsv | sed 's/\_.*\:chr/:chr/' > $outdir"/targets_Shared.tsv"



