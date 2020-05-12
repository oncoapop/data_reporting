#!/bin/sh

# Name of Project
Project="TNBC-SS"
# MiSeq reads
reads=200 # for latest single cell runs Aug 2014

# Project Directory
addr="damian@damianeva.org"

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"

# Example of full path and name of output from upstream pipeline
#/share/lustre/asteif/pipelines/primer_design/breast_xenograft/SA494/SA494_primer3_output.txt

# each sample name (space separated, no quotes, punctuations)
# for sample in DG1136g # run 26 Nov 2013
# for sample in DG1136g-hcm # run 27 Nov 2013
# for sample in DG1136g-dloh4 # run 27 Nov 2013
# Formerly used for TITAN runs
# for sample in SA221        # run 14 Aug 2014 
for sample in SA299        # run 14 Aug 2014
# for sample in SA425 # run 14 Aug 2014


#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

	cp $inpath"/"$sample"_p3_output.txt" $p3dir"/"$sample"_primer3_output.txt"

	# if run from old v1.5_pipeline splice
	posdir=$dir"/"$Project"/positions"
	cat $posdir"/"$sample"_positions.txt" | tr -d '"' | awk -F, '{print $4","$2":"$3"-"$3}' | sed 's/chr//2' | tail -n+2 > $posdir"/"$sample"_p3_positions.txt"

	~/Scripts/v3.1_pipeline/v3.1_display_pipeline.sh $Project $sample $reads
	echo $Project `date` >> ~/Projects/PrimerDesign/jobs.log; echo $Project $sample "completed using v3.1." >> ~/Projects/PrimerDesign/jobs.log

	echo $sample "Completed." | mail $addr -s "Message from v3.1 Display pipeline"

	done

echo "v3.1 Display Pipeline Completed." | mail $addr -s "Message from v3.1 Display pipeline"
