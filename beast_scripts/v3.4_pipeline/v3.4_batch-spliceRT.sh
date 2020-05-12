#!/bin/sh

# Name of Project
Project="eIF4A3"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/Projects/PrimerDesign/manual"

# Example of full path and name of output from upstream pipeline
# /home/dyap/Projects/PrimerDesign/manual/eIF4A3_NMD_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# Need MAY also to change line 160 of the primer_summaryRT2 to reflect the 
# the name change, if the output fails (do NOT include FS ="_")
 
# 21 Feb 2017 # eIF4A3_NMD
for sample in eIF4A3_NMD

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"_primer3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v3.4_pipeline/v3.4_display_pipelineRT.sh $Project $sample $reads
	echo $Project `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.4" >> ~/Projects/PrimerDesign/jobs.log

#	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"

