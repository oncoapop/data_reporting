#!/bin/sh

# Script to grep FLD primer name from read barcodes

# Working Directory
dir="/home/dyap/Projects/Single_Cell/ITH_case2_QC"

# Source and Output files where Barcoded files stored
sourcefile="/home/dyap/Projects/common/FLD_barcodes.csv"

inputfile=$dir"/barcodes-A4RND.txt"
samplesheet=$dir"/SampleSheet.csv"
outfile=$dir"/Matched_barcodes.csv"


cd $dir
rm -f $outfile

	for i in `cat $inputfile | awk -F"\t" '{print $1}'`
	do 	
	echo $i	
	{

	if [[ "$i" == "Index" ]];
		then continue
	fi

	echo -n `grep $i $sourcefile | awk -F"," '{print $1","$2","}' `;
	echo -n `grep $i $inputfile | awk -F"\t" '{print $3","} '`;
	grep $i $samplesheet | awk -F"," '{print $2} ';
	}  >> $outfile

 	done;
	
echo "Press RETURN to cont..." ; read ans


exit;

