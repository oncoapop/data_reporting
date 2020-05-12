#!/bin/bash

# QC Script to grep primer name from read barcodes
# Runs with high number of unindexed reads will need this QC Step

# Run ID
run="170105_M02348_0134_000000000-APJJ6"
# Working Directory
dir="/share/lustre/archive/MiSeq/MiSeq_Output_Files/"$run"/Data/Intensities/BaseCalls"

# Source and Output files where Barcoded files stored
#sourcedir=$dir"/Alignment"
sourcedir=$dir"/Alignment2"
outdir="/home/dyap/Projects/MiSeq_Data/QC/170105_M02348_0134_000000000-APJJ6"

if [ ! -d "$outdir" ] 
	then
		mkdir $outdir
fi

inputfile=$sourcedir"/DemultiplexSummaryF1L1.txt"
samplesheet=$sourcedir"/SampleSheetUsed.csv"
#samplesheet=$outdir"/Corrected_Samplesheet.csv"

outfile=$outdir"/Matched_barcodes_new.csv"


rm -f $outfile

	for i in `cat $inputfile | tail -n +34 | awk -F"\t" '{print $1}'`
	do 	

	if [[ "$i" == "Index" || "$i" == "Index2" ]];
		then continue
		echo "skip"
	fi


	I7=`grep --color=always $i $inputfile | awk -F"\t" '{print $1}' `
	I5=`grep --color=always $i $inputfile | awk -F"\t" '{print $2}' `
	nr=`grep --color=always $i $inputfile | awk -F"\t" '{print $3}' `

	echo $i"="$nr 
	echo $i"="$nr >> $outfile

	match=`grep --color=always $i $samplesheet | awk -F, '{print $2","$5","$6","$7","$8}'`
	echo $match

	if [[ $match != "" ]]
		then
		echo $match"\t\t"$I7,$I5,$nr
		echo $match"\t\t"$I7,$I5,$nr >> $outfile
		echo "=============================" 


	fi

 	done;
	
echo "Press RETURN to cont..." ; read ans


exit;

