#!/bin/sh

# Script to use the latest generated AmpliconManifest files
# at 

indir="/share/lustre/projects/breast_xeno_evolution/binomial_validation/amplicon_manifest_files"

# and SuppleFigfile.csv in 

wd="/home/dyap/public_html"
sample="SA429"
ampsample="SA429-PrimerSet1"
# SA429 and SA501 have 2 sets
# SA429,493,494,495,496,499,500,501,530,531,532,533,534,535,536,542
#
#SA429,SA501,SA542
sufile="/home/dyap/public_html/Tumour_Xenograft_Rev-set2/SA429"
#sufile="/home/dyap/Projects/PrimerDesign/Tumour_Xenograft_Rev-test/primer3"

outfile="/home/dyap/Projects/PrimerDesign/SupplFig/"$sample"_Suppl_Table.csv"
failfile="/home/dyap/Projects/PrimerDesign/SupplFig/"$sample"_Fail_Table.csv"
rm -f $outfile
rm -f $failfile

list=`ls $indir | awk -F. '{print $1}' | sed 's/temp//'`

	echo $ampsample".AmpliconManifest"
	pos=`cat $indir"/"$ampsample".AmpliconManifest" | awk -F: '{print $2}'`
	
	echo "==============================="
	echo $pos
	echo "==============================="

	for j in $pos
	do

	echo $ampsample

	inref=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F: '{print $3}'| awk -F"_" '{print $2}'`
	inst=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F: '{print $5}'| awk -F"_" '{print $2}'`
	inan=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F: '{print $6}'| awk -F"_" '{print $2}'`
	ingn=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F: '{print $6}'| awk -F"_" '{print $3}' | awk -F" " '{print $1}'`

	inchr=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F"\t" '{print $2}'`
	insta=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F"\t" '{print $3}'`
	inend=`grep $j $indir"/"$ampsample".AmpliconManifest" | awk -F"\t" '{print $4}'`

	prchr=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $2}'`
	prsta=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $3}'`
	prend=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $4}'`
	prleft=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $5}'`
	prright=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $6}'`
	prlen=`grep $j $sufile/$sample"_SuppleFigFile.csv" | awk -F, '{print $7}'`
	
	echo $inchr  " = " $prchr
	echo $insta  " = " $prsta
	echo $inend  " = " $prend

		if [[ $inchr == $prchr ]] && [[ $insta == $prsta ]] && [[ $inend == $prend ]];
			then
			echo $sample","$j","$inref","$inan","$ingn","$insta","$inend","$prleft","$prright","$prlen >> $outfile
			else
			echo $sample","$j","$inref","$inan","$ingn","$insta","$inend >> $failfile

		fi


	done


