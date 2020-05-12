#!/bin/sh

# Name of Project
Project="Tumour_Xenograft_Rev"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/share/lustre/asteif/pipelines/primer_design/breast_xenograft"

# Example of full path and name of output from upstream pipeline
#/share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA533/SA533_primer3_output.txt
#/share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA530/SA530_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# 17 Oct 2013 run
#for sample in SA532 SA533 SA534 SA535 SA536
# 12 Nov 2013 run
for sample in SA542

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"/"$sample"_primer3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v2.8_pipeline/v2.8_display_pipeline.sh $Project $sample $reads
	echo $Project `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v2.8." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"
