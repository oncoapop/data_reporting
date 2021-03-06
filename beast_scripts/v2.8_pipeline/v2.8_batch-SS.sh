#!/bin/sh

# Name of Project
Project="Single-Cell-v2-3"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/share/lustre/asteif/pipelines/primer_design/breast_xenograft"

# Example of full path and name of output from upstream pipeline
#/share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA494/SA494_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
for sample in SA494 SA495

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"/set3/"$sample"_primer3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v2.8_pipeline/v2.8_display_pipeline.sh $Project $sample $reads
	echo $Project `date` >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"
