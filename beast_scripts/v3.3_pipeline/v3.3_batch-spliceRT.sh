#!/bin/sh

# Name of Project
Project="Splice"
# MiSeq reads
reads=150
# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/Projects/PrimerDesign/manual"

# Example of full path and name of output from upstream pipeline
# /home/dyap/Projects/PrimerDesign/manual/hct116_htert_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# Need MAY also to change line 160 of the primer_summaryRT2 to reflect the 
# the name change, if the output fails (do NOT include FS ="_")
 
# 28 Apr 2014 run
# for sample in hct116_htert
# 19 Oct 2015 # original set
# for sample in Splice-RT
# 25 Oct 2015 # assoc factors
# for sample in Splice-RT2
# 16 Mar 2016 # new factors
for sample in Splice-RT3
# 6 Apr 2016 # CLK1,2
for sample in Splice-RT3

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

	~/Scripts/v3.3_pipeline/v3.3_display_pipelineRT.sh $Project $sample $reads
	echo $Project `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.3." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from Display pipeline"

	done

echo "Display Pipeline Completed." | mail $addr -s "Message from Display pipeline"

