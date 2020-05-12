#!/bin/sh

# script to compare siRNA genes 

wd="/home/dyap/Projects/EIF4A3_paper/comparisons"

cd $wd 

awk '{if(sqrt($10*$10) >= 0.1 && $11 > 10) print $24, $3, $2, $10, $11}' TK_18_19_20_vs_TK_17_21_hg19_miso_bf.filtered.summary.annotated > 24h_genes

awk '{if(sqrt($10*$10) >= 0.1 && $11 > 10) print $24, $3, $2, $10, $11}' TK_22_23_24_vs_TK_17_21_hg19_miso_bf.filtered.summary.annotated > 48h_genes

awk '{print $1}' 48h_genes | sort | uniq -c | awk '{print "48h", $1, $3}' >> siRNA_new


cd /home/dyap/Projects/EIF4A3_paper/toTransfer/EIF4A3project/from_the_other_paper

# comparing with siRNAs:
./plot_venn.py --inputFiles Miso_EIF4A3.genes DiffSplice_EIF4A3.genes ../Miso_siRNA/48h_genenames --inputNames Miso_EIF4A3_others Diffsplice_EIF4A3_others $wd"/siRNA_new" --inputCols 1 --outputpref MisoEIF4A3Others_DiffSpliceOthers_siRNAOur48h

