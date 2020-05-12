#!/bin/sh

# Name of Project
Project="Tumour_Xenograft"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/share/lustre/asteif/primer_design/breast_xenograft"

# Example of full path and name of output from upstream pipeline
#/share/lustre/asteif/primer_design/breast_xenograft/SA533/SA533_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
for sample in SA534 SA535
# for sample in SA532 SA533 SA534 SA535 SA536

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

# for live version uncommment this line
#	cp $inpath"/"$sample"/"$sample"_primer3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

#	nohup ~/Scripts/v4.2_pipeline/v4.2_display_pipeline.sh $Project $sample $reads && echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log
	~/Scripts/v4.2_pipeline/v4.2_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.2." >> ~/Projects/PrimerDesign/jobs.log

	done

