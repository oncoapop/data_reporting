#!/bin/sh

# Scipt to map qPCR plate with primer plate and also with the amplicon name, length as well as cal tm

# Enter the name of sample to run 
name="SA500"

platedir="/home/dyap/Projects/Tumour_Evol/"$name"/QC"

amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/"$name"/check-positions-QC/"$name"_design_space.csv"

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Tumour_Evol/"$name"/check-positions-QC/"

# Name processing - Do not change

# Source and Output directories where working files are stored
sourcedir=$dir"SNV/"

# contains the list of positions
infile=$sourcedir$name"_pos.txt"
# contains the output of this script
outfile=$dir$name"_Manifest"
wtpos=$dir$name"_WT_positions.csv"
manifest=$dir$name".AmpliconManifest"

readfile=$name"_anno.txt"

# This gets the WT position in the middle of a 11bp sequences for matching SNV
echo Processing files
cat $wtpos | awk -F, '{print $2"_"$3,$4}' | tr -d '"' | tr "c" "C" > $outfile"66.tmp"

	for i in `cat $wtpos | awk -F, '{print $3}' | tr -d '"'`
		do
			chr=`grep $i $wtpos | awk -F, '{print $2}' | tr -d '"'`
			grep $i $manifest 
		done
 


rm $outfile*

exit
