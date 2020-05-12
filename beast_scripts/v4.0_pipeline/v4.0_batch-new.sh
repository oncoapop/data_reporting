#!/bin/sh

# Name of Project
Project="FL"
# MiSeq reads
reads=150
# Project Directory

basedir="/home/dyap/Projects/PrimerDesign"
dir=$basedir"/"$Project
p3dir=$basedir"/"$Project"/primer3"
inpath=$dir

# Example of full path and name of output from upstream pipeline
# /home/dyap/Projects/PrimerDesign/FL/
# FL1006T1_FL1006T2.validations.T2.primerInput.txt_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
#for sample in FL1004_v1 FL1004_v2 FL1006_v1 FL1006_v2
for sample in FL1006_v1 

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log
         echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir
			mkdir $p3dir
	fi

# for live version uncommment this line
	export samples=$sample
	~/Scripts/v4.0_pipeline/v4.0_parse2primer3.sh 

#	nohup ~/Scripts/v4.0_pipeline/v4.0_display_pipeline.sh $Project $sample $reads && echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log
	~/Scripts/v4.0_pipeline/v4.0_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v4.0." >> ~/Projects/PrimerDesign/jobs.log

	done

