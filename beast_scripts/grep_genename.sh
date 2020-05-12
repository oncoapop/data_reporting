#!/bin/sh

# Script to grep gene name from R output from MS data

wd="/home/dyap/Projects/Takeda_T3/MS data/Expt1_ectopicCLK2"
cd "$wd"
pwd
for i in `ls ek_20141210_CLK*`
	do
	echo $i
	sample=`echo $i | awk -F_ '{print $3"_"$5}'`
	name=$sample"-genes"
	cat $i | awk -F"GN=" '{print $2}'| awk -F" " '{print $1}' | sort -u | sed '/^$/d' > $name
	done
