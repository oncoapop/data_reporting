#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis"

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo "Started summarizing for EXPT-"$expt >> $dir/macsjobs.log

	cd $dir"/"$Project"/EXPT-"$expt

	for f in `ls *_model.r`
		do
		/usr/bin/Rscript $f		
		done

	done

exit
##############

