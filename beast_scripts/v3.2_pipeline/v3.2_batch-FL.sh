#!/bin/sh

# Name of Project
Project="FL"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"
inputsamples="FL1004_v1 FL1004_v2 FL1006_v1 FL1006_v2"

export samples=$inputsamples
#~/Scripts/v3.2_pipeline/v3.2_parse2primer3.sh

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

# for live version uncommment this line
	cp $inpath"/"$sample"_primer3_output.txt.cat" $p3dir"/"$sample"_primer3_output.txt"

#	nohup ~/Scripts/v3.2_pipeline/v3.2_display_pipeline.sh $Project $sample $reads && echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log
	~/Scripts/v3.2_pipeline/v3.2_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v3.2." >> ~/Projects/PrimerDesign/jobs.log

	done

