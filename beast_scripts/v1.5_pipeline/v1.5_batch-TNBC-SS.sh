#!/bin/sh

# Name of Project
Project="TNBC-SS"
# MiSeq reads
reads=200
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
posdir=$dir"/"$Project"/positions"

rm -f $posdir"/"*.tmp*
rm -f $posdir"/"*-hg18*
#rm -f $posdir"/"*-hg19
#rm -f $posdir"/"*-hg19.csv

# each sample name (space separated, no quotes, punctuations)
# These are TNBC hg18 positions!
build="hg18"
# Version of primer3 settings
ver=v3ss
#for sample in SA999 # testing

# for sample in SA425        # 14 Aug 2014 
# for sample in SA221        # 14 Aug 2014 
for sample in SA299   # 14 Aug 2014 

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi


	cat $posdir"/"$sample"_amplicon_nameField.txt" | sed 's/_/-/g' | tail -n+2 | uniq > $posdir"/"$sample".tmp1"

	for i in `cat $posdir"/"$sample".tmp1" | awk -F: '{print $2}' | awk -F- '{print $2}'`

		do
#		name=`grep $i $posdir"/"$sample".tmp1"`
		name=`grep $i $posdir"/"$sample".tmp1"  | awk -F: '{print $1}'`
		chr=`grep $i $posdir"/"$sample".tmp1"  | awk -F: '{print $2}' | awk -F- '{print $1}' | sed 's/chr//g'`
		pos=`grep $i $posdir"/"$sample".tmp1"  | awk -F: '{print $2}' | awk -F- '{print $2}'`

		# Format wanted
		# Sample-name,chromsome_#:start-end (if snv start=end)
		
		echo $name","$chr":"$pos"-"$pos >> $posdir"/"$sample"-"$build
		echo $name","$chr":"$pos"-"$pos 

		done

	echo "Sample "$sample" done."

	rm -f $posdir"/"*.tmp*

	# Run pipeline v1.5
	perl ~/Scripts/v1.5_pipeline/v1.5_primerdesign.pl Project=$Project Sample=$sample Ref=$build Miseq=$reads p3settings=$ver

	done

echo "Primer3 v1.5 Pipeline Completed for "$sample | mail $addr -s "Message from Primer3 pipeline"





exit;

