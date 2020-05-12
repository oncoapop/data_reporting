#!/bin/sh

# Name of Project
Project="Tumour_Xenograft"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath=$dir"/manual"

# each sample name (space separated, no quotes, punctuations)
for sample in SA535-BRCA

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log
	 echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"_p3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v3_pipeline/v3_display_pipeline.sh $Project $sample $reads

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log
	 echo $sample "completed." >> ~/Projects/PrimerDesign/jobs.log

	done
