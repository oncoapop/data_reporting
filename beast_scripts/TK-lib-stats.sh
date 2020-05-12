#!/bin/sh

# Script to automatically get the bam stats from flagstats
# source 

pwd="/home/dyap/Projects/eIF4A3_NMD"
pwd2="/home/dyap/Projects/eIF4A3_NMD/flagstats"
infile=$pwd"/IDs"
outfile=$pwd"/ID-stats"

echo -e "SA_ID"'\t'"GSC_lib"'\t'"Analysed_Lib" > $outfile

SAID=`cat $infile | awk -F"\t" '{print $2}' | tail -n +2`

echo "|| S/NO || Library || Total Reads || Aligned Reads || % Aligned Reads || STAR Total Reads || STAR Aligned reads || "

counter=1
for i in $SAID
	do

	# Test to see if the path and file exist and can be read

	path="/share/lustre/archive/"$i"/illumina_wtss/"$i"/gsnap_aligned"

	if [ -d "$path" ];
		then testfile=`ls $path/*.flagstat`
	fi 

	if [ -r "$testfile" ];
		then file=$testfile

			# Grep the no of reads and % aligned
			reads=`cat $file | grep "QC-passed" | awk -F" " '{print $1}'`
			align=`cat $file | grep -m1 "mapped" | awk -F" " '{print $1}'`
			alignpercent=`cat $file | grep -m1 "mapped" | awk -F"(" '{print $2}' | awk -F":" '{print $1}'`

			formatreads=`printf "%'d" $reads`
			formatalign=`printf "%'d" $align`

		else file="NONE"
	fi 

	testfile="/home/dyap/Projects/eIF4A3_NMD/flagstats/"$i".flagstat"

	if [ -r "$testfile" ];
		then file=$testfile

			# Grep the no of reads and % aligned
			starreads=`cat $file | grep "QC-passed" | awk -F" " '{print $1}'`
			staralign=`cat $file | grep -m1 "mapped" | awk -F" " '{print $1}'`
			sanpalignpercent=`cat $file | grep -m1 "mapped" | awk -F"(" '{print $2}' | awk -F":" '{print $1}'`

			formatstarreads=`printf "%'d" $starreads`
			formatstaralign=`printf "%'d" $staralign`

		else file="NONE"
	fi 

	# format numbers using printf "%'d" $number
	echo "| "$counter" |"$i" | "$formatreads" | "$formatalign" | "$alignpercent" | "$formatstarreads" | " $formatstaralign" |"


	# clear
	formatreads=""
	formatalign=""
	alignpercent=""
	formatstaralign=""
	formatstarreads=""
	counter=$((counter+1))

	done

exit;

