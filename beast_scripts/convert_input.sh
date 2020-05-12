#!/bin/sh

# Script to convert input from chr:pos in filename - sample-raw
# To Sample_ID_chr_pos

Project="Tumour_Xenograft_Rev"
dir="/home/dyap/Projects/PrimerDesign/"$Project"/positions"
suffix="raw"
hg="hg19"

ls -al $dir

for i in `ls $dir | awk -F- '{print $1}'`
	do
	cd $dir
	awk -F" " '{print FILENAME","$1":"$2"-"$2}' $i"-"$suffix | sed 's/-raw//g' > $dir"/"$i"-"$hg
	cat $dir"/"$i"-"$hg
#	read ans
	done

exit;

