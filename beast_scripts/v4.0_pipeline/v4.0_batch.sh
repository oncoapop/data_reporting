#!/bin/sh

# Name of Project
Project="FL"
# Name of Lab
lab="Steidl"
folder="2014Mar25"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$lab"/"$Project"/"$folder
inpath=$p3dir

# Example of full path and name of output from upstream pipeline
# /home/dyap/Projects/PrimerDesign/Steidl/FL/2014Mar25/
# FL1006T1_FL1006T2.validations.T2.primerInput.txt_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
for sample in FL1006T1_FL1006T2.validations.T2.primerInput.txt FL1006T1_FL1006T2.validations.T1.primerInput.txt FL1004T1_FL1004T2.validations.T2.primerInput.txt FL1004T1_FL1004T2.validations.T1.primerInput.txt


#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

#	if [ -d $dir"/"$Project ]; then
#			mkdir $p3dir
#		else
#			mkdir $dir"/"$Project
#			mkdir $p3dir
#	fi

# for live version uncommment this line
#	cp $inpath"/"$sample"/"$sample"_primer3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

#	nohup ~/Scripts/v4.0_pipeline/v4.0_display_pipeline.sh $Project $sample $reads && echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log
	~/Scripts/v4.0_pipeline/v4.0_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log

	done

