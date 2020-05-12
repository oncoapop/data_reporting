#!/bin/sh

# Name of Project
Project="TITAN-SS"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"

# Example of full path and name of output from upstream pipeline
# /home/dyap/Projects/PrimerDesign/TITAN-SS/primer3/DG1136g_primer3_output.txt 

# each sample name (space separated, no quotes, punctuations)
for sample in DG1136g

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"_p3_output" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v3.2_pipeline/v3.2_display_pipeline.sh $Project $sample $reads
	echo $Project `date` >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.2." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"



exit;

