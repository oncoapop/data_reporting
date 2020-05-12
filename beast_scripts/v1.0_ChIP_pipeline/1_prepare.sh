#!/bin/sh

# Name of Project
Project="CX5461"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis"
inpath="/share/lustre/projects/chip-seq/APARICIO-"
suffix="_bt2_alignment/outputs/results/bam"

# Example of full path and name of output from upstream pipeline
# /share/lustre/projects/chip-seq/APARICIO-218/IgG-CX54-10-6-R4_S8_bt2_alignment/outputs/results/bam

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 216 217 218 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \n >> $dir/jobs.log ; echo $sample "started processing." >> $dir/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $dir"/"$Project"/EXPT-"$expt
		else
			mkdir $dir"/"$Project
			mkdir $dir"/"$Project"/EXPT-"$expt
	fi


	if [[ $expt == "216" ]];
		then
		macsdir=$inpath$expt"/ADHVY"
		else
		macsdir=$inpath$expt
	fi		

	cd $macsdir

	samples=`ls | sed 's/_bt2_alignment//' | sed '/yaml$/d'`

		for directory in $samples

		do
			
			fullpath=$macsdir/$directory$suffix
			cd $fullpath
			cp *.bam $dir"/"$Project"/EXPT-"$expt
			cp *.bai $dir"/"$Project"/EXPT-"$expt
		done


	done

exit
##############
