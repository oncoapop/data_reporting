#!/bin/sh

# Name of Project
Project="FL"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"

# Input the names of the samples, they are processed in series
# To do - rewrite to make them run in parallel
# 1 Apr 2014
inputsamples="FL1004-v1 FL1004-v2 FL1006-v1 FL1006-v2"

export samples=$inputsamples
~/Scripts/v4.1_pipeline/v4.1_parse2primer3.sh

# each sample name (space separated, no quotes, punctuations)
for sample in $samples
 

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

# copies the final version of the primer design process to the primer3 project directory
	cp $inpath"/"$sample"_primer3_output.txt.cat" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v4.1_pipeline/v4.1_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v4.1 ." >> ~/Projects/PrimerDesign/jobs.log

	done

