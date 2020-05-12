#!/bin/sh

# Script to automatically get the bam stats from flagstats
# source 

pwd="/home/dyap/Projects/eIF4A3_NMD"
pwd2="/home/dyap/Projects/eIF4A3_NMD/flagstats"
infile=$pwd2"/table.txt"
outfile=$pwd"/ID-stats"

SAID=`cat $infile | awk -F" " '{print $1}'`

counter=1
for i in $SAID
	do

	# Test to see if the path and file exist and can be read

	testfile="/home/dyap/Projects/eIF4A3_NMD/flagstats/"$i".flagstat"

	if [ -r "$testfile" ];
		then file=$testfile

			# Grep the no of reads and % aligned
			starreads=`cat $file | grep "QC-passed" | awk -F" " '{print $1}'`
			staralign=`cat $file | grep -m1 "mapped" | awk -F" " '{print $1}'`
			sanpalignpercent=`cat $file | grep -m1 "mapped" | awk -F"(" '{print $2}' | awk -F":" '{print $1}'`

			# format numbers using printf "%'d" $number
			formatstarreads=`printf "%'d" $starreads`
			formatstaralign=`printf "%'d" $staralign`

		else file="NONE"
	fi 

			line=`grep $i $infile | tr -d '\'`

			echo $line" & "$formatstarreads" \\\\"
	# clear
	formatstaralign=""
	formatstarreads=""
	counter=$((counter+1))

	done

exit;

