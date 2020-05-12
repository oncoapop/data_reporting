#!/bin/bash

wd="/home/dyap/ind231"
cd $wd

file="expanded_list_grouped_V5b.csv"
formatted="sift_input_select.csv"
audit="sift_input_audit2.csv"

rm -f $formatted
rm -fr $audit

# Remove the header
# select fro n < 5
#for i in `cat $file | tail -n +2 | awk -F, '{if ($15 < 5) print $0}'`
for i in `cat $file | tail -n +2 | awk -F, '{if ($2 == "BRCA1" || \
$2 == "BRCA2" || $2 == "PALB2" || $2 == "TP53" || $2 == "RAD51" ) print $0}'`
	do
	chr=`grep "$i" $file | awk -F"," '{print $1}'`
	gene=`grep "$i" $file | awk -F"," '{print $2}'`
	start=`grep "$i" $file | awk -F"," '{print $3}'`
	end=`grep "$i" $file | awk -F"," '{print $4}'`
	ref=`grep "$i" $file | awk -F"," '{print $5}'`
	alt=`grep "$i" $file | awk -F"," '{print $6}'`
	strand=`grep "$i" $file | awk -F"," '{print $7}'`


	# insertion
	if [ "$start" -eq "$end" ]
		then
		nt=`echo $alt | cut -b2-`
		echo $nt
		else
		nt="/"
	fi

	echo $chr","$start","$end","$strand"1,"$ref","$alt
	echo $chr","$start","$end","$strand"1,"$ref","$alt >> $audit
	echo $chr","$start","$end","$strand"1,"$nt >> $formatted
	echo $chr","$start","$end","$strand"1,"$nt
	echo $chr","$start","$end","$strand"1,"$nt >> $audit
	echo "======"
	echo "======" >> $audit

	done

split -d -l 10 $formatted submit

exit

