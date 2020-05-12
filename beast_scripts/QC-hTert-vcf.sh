#!/bin/sh

# Script to QC failed posistions with the MiSeq vcf output

# Base directory
dir="/home/dyap/Projects/Single_Cell/184-hTert"
# Where the VCFs are lept (this is downloaded from BaseSpace)
subfol="/home/dyap/Projects/Single_Cell/184-hTert/3530550"

# summary file
outfile=$dir"/summary.txt"

echo "Query Position	No of samples with at least one variant read	Manifest Entry" > $outfile

for i in `cat $dir"/"failed-pos.txt | awk -F"_" '{print $2}'`
	do
	count=`grep $i $subfol/*.vcf | wc -l`
	id=`grep $i $subfol/*.vcf | awk -F"\t" '{print $1"_"$2}' | awk -F":" '{print $2}' | sort -u`
	manifest=`grep $i /share/lustre/archive/MiSeq/BaseSpace/3530550_184-hTert/A_184-hTert.AmpliconManifest | awk -F"\t" '{print $1}'`
	echo $id	$count	$manifest>> $outfile
	echo -ne "#"
	done

