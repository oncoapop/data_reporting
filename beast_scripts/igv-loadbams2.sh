#!/bin/sh

# Script to generate commands for loading bams into igv

# /share/lustre/abashash/EIF4A3/star_fusion_cufflinks_pipeline/OUTPUT/RUN/TK_22_23_24_star_cufflinks/outputs/results

input="/home/dyap/Projects/eIF4A3_NMD/RNAifiles"
output="/home/dyap/Projects/eIF4A3_NMD/igv-RNAi-input.txt"

echo "new" > $output
echo "genome hg19" >> $output
echo "snapshotDirectory /home/dyap/Projects/eIF4A3_NMD/Plots" >> $output

# get the control file
#line="HCT-116"
line="HeLa"
TKcon=`cat $input | awk -F"\t" '$2 == "0"' | grep $line | awk -F"\t" '{print $3}'`
echo $TKcon

	confile="/share/lustre/abashash/EIF4A3/star_fusion_cufflinks_pipeline/OUTPUT/RUN/"$TKcon"_star_cufflinks/outputs/results/*.bam"
	bam=`find $confile`
	echo "load "$bam >> $output
	echo $bam

TK=`cat $input | awk -F"\t" '$2 != "0"' | grep $line | awk -F"\t" '{print $3}'`
echo $TK

for i in `echo $TK`
	do
	file="/share/lustre/abashash/EIF4A3/star_fusion_cufflinks_pipeline/OUTPUT/RUN/"$i"_star_cufflinks/outputs/results/*.bam"
	bam=`find $file`
	echo "load "$bam >> $output
	echo $bam
	done

#gene="HRAS"
#gene="MAPK13"
#gene="G3BP1"
#gene="BCL2L1"
#gene="AURKB"
echo "goto "$gene >> $output
echo "sort" >> $output
echo "snapshot "$line"_"$drug"_"$gene >> $output




