#!/bin/sh

# This Script was written by Damian Yap (May 2013)
# This only works of single plex primer amplification!

# When called within a script it uses the file /home/dyap/dyap_temp/Temp_cal-input
# Temp_cal-input is in one-line "name seq"  format

file="/home/dyap/dyap_temp/Temp_cal-input"
outfile="/home/dyap/dyap_temp/Temp_cal-output"
echo "" > $outfile

	for i in `cat $file | awk -F" " '{print $1}' | sort -u`
	do

	input=`grep $i $file | awk -F" " '{print $2}'` 
	amplicon=`echo $input | wc -c`
        A=`echo $input | grep "A" -o | wc -l`
        T=`echo $input | grep "T" -o | wc -l`
        C=`echo $input | grep "C" -o | wc -l`
        G=`echo $input | grep "G" -o | wc -l`
        tm=`echo "(( 64.9 + (41 * ($C+$G-16.4))/($A+$T+$C+$G) ))" | bc`

	id=`grep -m1 $i $file | awk -F" " '{print $3}'` 
	echo $i"_Tm="$tm"_Length="$amplicon"_ID="$id >> $outfile

	test=`grep -m1 $i $file | awk -F" " '{print $4}'` 
	if [ ! -z "$test" ]
		then 
		echo "These primer sets are duplicated:  "$test,$id	
		echo $i"_Tm="$tm"_Length="$amplicon"_ID="$test >> $outfile
	fi
	
	done

exit;

