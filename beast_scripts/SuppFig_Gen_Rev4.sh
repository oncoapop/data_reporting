#!/bin/sh

# Script to use the latest generated AmpliconManifest files
# at 

indir="/share/lustre/projects/breast_xeno_evolution/binomial_validation/amplicon_manifest_files"

# and SuppleFigfile.csv in 

# SA429 and SA501 have 2 sets
# SA429,493,494,495,496,499,500,501,530,531,532,533,534,535,536,542


list=`ls $indir | awk -F. '{print $1}' | sed 's/temp//'`

for i in $list
	do
	echo $i".AmpliconManifest"
	sample=`echo $i | awk -F- '{print $1}'`
	outfile="/home/dyap/Projects/PrimerDesign/SupplFig/"$i".csv"
	echo "Sample,SNV position,hg19 REF base,Status,Type,Annotation,Amplicon Chromosome,Left Primer Start,Left Primer End,Right Primer Start,Right Primer End" > $outfile

	pos=`cat $indir"/"$i".AmpliconManifest" | tail -n +7 | awk -F: '{print $2}'`
	echo $i

#	echo "==============================="
#	echo $pos
#	echo "==============================="

	for j in $pos
	do

	inref=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $3}'| awk -F"_" '{print $2}'`
	inst=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $5}'| awk -F"_" '{print $2}'`
	inan=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $6}'| sed 's/-/_/g' | awk -F"_" '{print $2}'`
	ingn=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $6}'| sed 's/-/_/g' | awk -F"_" '{print $3}' | awk -F" " '{print $1}' | sed 's/,/-/g'`


	inchr=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $2}'`
	insta=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $3}'`
	inend=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $4}'`
	leftlen=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $5}'`
	rightlen=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $6}'`

	inleft=`echo "$insta + $leftlen -1" | bc`
	inright=`echo "$inend - $rightlen +1" | bc`

		echo $sample","$j","$inref","$inst","$inan","$ingn","$inchr","$insta","$inleft","$inright","$inend >> $outfile

	done

done
