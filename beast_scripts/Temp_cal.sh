#!/bin/sh

# This Script was written by Damian Yap (May 2013)

# When called within a script it uses the file /home/dyap/dyap_temp/Temp_cal-input
# Temp_cal-input is in one-line "name seq"  format

file="/home/dyap/dyap_temp/Temp_cal-input"
outfile="/home/dyap/dyap_temp/Temp_cal-output"
echo "" > $outfile

	for i in `cat $file | awk -F" " '{print $1}'`
	do

	input=`grep $i $file | awk -F" " '{print $2}'` 
	amplicon=`echo $input | wc -c`
        A=`echo $input | grep "A" -o | wc -l`
        T=`echo $input | grep "T" -o | wc -l`
        C=`echo $input | grep "C" -o | wc -l`
        G=`echo $input | grep "G" -o | wc -l`
        tm=`echo "(( 64.9 + (41 * ($C+$G-16.4))/($A+$T+$C+$G) ))" | bc`
	id=`grep $i $file | awk -F" " '{print $3}'` 

	echo $i"_Tm="$tm"_Length="$amplicon"_ID="$id >> $outfile
	done

exit;

