#!/bin/sh

# Script to use the latest generated AmpliconManifest files
# at 

indir="/share/lustre/projects/breast_xeno_evolution/binomial_validation/amplicon_manifest_files"

# and SuppleFigfile.csv in 

wd="/home/dyap/public_html"

outfile="/home/dyap/Projects/Tumour_Evol/Suppl_Table.csv"
failfile="/home/dyap/Projects/Tumour_Evol/Fail_Table.csv"
rm $outfile
rm $failfile

list=`ls $indir | awk -F. '{print $1}' | sed 's/temp//'`

for i in $list
	do
	echo $i".AmpliconManifest"
	sample=`echo $i | awk -F- '{print $1}'`
	pos=`cat $indir"/"$i".AmpliconManifest" | tail -n +7 | awk -F: '{print $2}'`
	
	echo "==============================="
	echo $pos
	echo $inref

	for j in $pos
	do

	echo $sample

	inref=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $3}'| awk -F"_" '{print $2}'`
	inst=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $5}'| awk -F"_" '{print $2}'`
	inan=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $6}'| awk -F"_" '{print $2}'`
	ingn=`grep $j $indir"/"$i".AmpliconManifest" | awk -F: '{print $6}'| awk -F"_" '{print $3}'`

	inchr=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $2}'`
	insta=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $3}'`
	inend=`grep $j $indir"/"$i".AmpliconManifest" | awk -F"\t" '{print $4}'`

	echo $sample","$j","$inref","$inan","$inst","$ingn","$inchr","$insta","$inleft","$inright","$inend >> $outfile

	done


done
