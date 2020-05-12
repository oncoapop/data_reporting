#!/bin/sh

# Name of Project
Project="Splice"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"
export database="/home/dyap/Projects/PrimerDesign/manual/splice3.2bit"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/Projects/PrimerDesign/manual"

# Example of full path and name of output from upstream pipeline
#/home/dyap/Projects/PrimerDesign/manual/Splice-sign1_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# 10 Feb 2015 run
for sample in Splice-sign2

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

	~/Scripts/v3.3_pipeline/v3.3_display_pipeline2.sh $Project $sample $reads
	echo $Project `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.3." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"
