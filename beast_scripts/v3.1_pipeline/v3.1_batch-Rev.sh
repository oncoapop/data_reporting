#!/bin/sh

# Name of Project
Project="Tumour_Xenograft_Rev-re"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/share/lustre/asteif/pipelines/primer_design/breast_xenograft"
# inpath="/share/lustre/asteif/pipelines/primer_design/single_cell"

# Example of full path and name of output from upstream pipeline
# /share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA533/SA533_primer3_output.txt
# /share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA530/SA530_primer3_output.txt
# /share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA499/SA499_primer3_output.txt
# /share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA499set2/SA499_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# 17 Oct 2013 run
# for sample in SA532 SA533 SA534 SA535 SA536
# 12 Nov 2013 run
# for sample in SA542
# 28 Nov 2013 run (re-design new set with old samples)
# for sample in SA429 SA501
# 19 Dec 2013 run 
# for sample in SA531
# 22 Jan 2014
# for sample in SA530
# 13 Jun 2014
# for sample in SA499
# 26 Jun 2014
set="set2"
for sample in SA499

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample$set"/"$sample"_primer3_output.txt" $p3dir"/"$sample$set"_primer3_output.txt"

	~/Scripts/v3.1_pipeline/v3.1_display_pipeline.sh $Project $sample$set $reads
	echo $Project `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"
