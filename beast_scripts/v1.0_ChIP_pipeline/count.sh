#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/Projects/ChIPseqAnalysis"

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 218

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo "Started summarizing for EXPT-"$expt >> $dir/macsjobs.log

	cd $dir"/"$Project"/EXPT-"$expt

	for f in `ls *.bedgraph`
		do
		a=($(sort -n $f | sed -n '1s/^\([0-9]\+\).*$/\1/p;$s/^\([0-9]\+\).*$/\1/p'))
		done

	done

exit
##############

