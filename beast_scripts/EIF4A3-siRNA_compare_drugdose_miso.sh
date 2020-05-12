#!/bin/sh

# script to compare siRNA genes with high and low doses of the drugs

wd="/home/dyap/Projects/EIF4A3_paper/comparisons/compare-genes"

cd $wd 

path="/home/dyap/Projects/EIF4A3_paper/toTransfer/EIF4A3project/from_the_other_paper"

# comparing with siRNAs:
$path/plot_venn.py --inputFiles Hela_202.highdose Hela_202.lowdose miso_si48_genes --inputNames "MISO_HeLa_high_dose" "MISO_HeLa_low_dose" "MISO_siRNA_48h" --inputCols 1 --outputpref Miso_HeLa_high-low_dose-SiRNA48

