#!/bin/sh

# This script checks the presence of low freq SNPs from all samples

# Working directory & files
dir="/home/dyap/Projects/Single_Cell/variants"
source="lowfreq.csv"
tmp=$dir"/snp_pos.tmp"
dest="SNP_freq.csv"
snpfreq="SNP-allfreq.csv"

echo "" > $tmp

# This makes a list of all positions that there are low frequencies in the VBA0038 run1 & 3
cd $dir
rm $dest
rm $snpfreq

grep "chr" $source | awk -F"," '{print $2}' | awk -F":" '{print $2}' | tr -d '"' > $tmp

	for i in `cat $tmp | sort -u`
		do
		{
		count=`grep $i *.txt | grep "VBA" | wc -l`
		hets=`grep $i *.txt | grep "0/1" | wc -l`
		homs=`grep $i *.txt | grep "1/1" | wc -l`
		freq=`grep $i *.txt | awk -F, '{print $8}'`
		echo $i","$count","$hets","$homs >> $dest
		echo $freq >> $snpfreq
		}
	done

exit;

 
