#!/bin/sh

# Script to generate commands for loading bams into igv

# /share/lustre/projects/takeda_EIF4A3/samples/TK0nnn/illumina_wtss/TK0nnn/gsnap_aligned/

input="/home/dyap/Projects/eIF4A3_NMD/files"
output="/home/dyap/Projects/eIF4A3_NMD/igv-input.txt"

echo "new" > $output
echo "genome hg19" >> $output
echo "snapshotDirectory /home/dyap/Projects/eIF4A3_NMD/Plots" >> $output

# get the control file
line="HCT-116"
#line="HeLa"
TKcon=`cat $input | awk -F"\t" '$3 == "0"' | grep $line | awk -F"\t" '{print $4}'`
drug="T-3787595"
#drug="T-3787598"

confile="/share/lustre/projects/takeda_EIF4A3/samples/"$TKcon"/illumina_wtss/"$TKcon"/gsnap_aligned/*.bam"
	bam=`find $confile`
	echo "load "$bam >> $output
	echo $bam

TK=`cat $input | awk -F"\t" '$3 != "0"' | grep $line | grep $drug | awk -F"\t" '{print $4}'`
#echo $TK


for i in `echo $TK`
	do
	file="/share/lustre/projects/takeda_EIF4A3/samples/"$i"/illumina_wtss/"$i"/gsnap_aligned/*.bam"
	bam=`find $file`
	echo "load "$bam >> $output
	echo $bam
	done

gene="MAPK13"
#gene="G3BP1"
#gene="BCL2L1"
#gene="AURKB"
echo "goto "$gene >> $output
echo "sort" >> $output
echo "snapshot "$line"_"$drug"_"$gene >> $output




