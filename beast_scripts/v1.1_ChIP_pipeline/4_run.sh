#!/bin/sh

# Name of Project
Project="CX5461"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis2"

# Example of full path and name of output from upstream pipeline
# /share/lustre/projects/chip-seq/APARICIO-218/IgG-CX54-10-6-R4_S8_bt2_alignment/outputs/results/bam

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 216 217 218 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \n >> $dir/jobs.log ; echo $sample "EXPT-"$expt" Started processing." >> $dir/jobs.log

		cd $dir"/"$Project"/EXPT-"$expt
		# using input as control
		Rscript /home/dyap/Scripts/v1.1_ChIP_pipeline/4.1_sushiPlot2.R
		# using IgG as control
		# Rscript /home/dyap/Scripts/v1.1_ChIP_pipeline/4_sushiPlot2.R

	 echo `date` \n >> $dir/jobs.log ; echo $sample "EXPT-"$expt" Finished processing." >> $dir/jobs.log


	done

/home/dyap/Scripts/v1.1_ChIP_pipeline/5_CopytoBeast.sh

exit
##############
