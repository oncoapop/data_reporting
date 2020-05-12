#!/bin/sh
# Script to generate PCR primers for supplemental figures
# Format 
# No SNV_ID , amplicon

# Set working directory and files
wdir="/home/dyap/Projects/Tumour_Evol/Summary"
cd $wdir

for i in `ls *isPCR_chk.csv`
	do
	sample=`echo $i | awk -F"_" '{print $1}' `
	outfile=$wdir"/"$sample"_SupplFig_allprimer.csv"
	
	for j in `cat $i`
		do
		chr=`echo $j | tr -d '"' | awk -F, '{print $2}' | awk -F"_" '{print $2}'`
		pos=`echo $j | tr -d '"' | awk -F, '{print $2}' | awk -F"_" '{print $3}'`
		start=`echo $j | tr -d '"' | awk -F, '{print $5}'`
		end=`echo $j | tr -d '"' | awk -F, '{print $6}'`
		len=`echo $j | tr -d '"' | awk -F, '{print $7}'`
		left=`echo $j | tr -d '"' | awk -F, '{print $8}' | tr -d " "`
		right=`echo $j | tr -d '"' | awk -F, '{print $9}' | tr -d " "`

		if [[ $pos != "" && $len != "NA" ]];
			then {
			      echo $sample,$chr"_"$pos,$start,$end,$left,$right >> $outfile
			     }
		fi
		echo "."	
		done

	done

