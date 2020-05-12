#!/bin/sh

# Script to move data from Sam's account to my account for visualization
dir="218-RAD51"

# Prepare the directory in my space
cd "/home/dyap/Projects/ChIPseqAnalysis/Sam analysis"
mkdir $dir
cd $dir

# copy all the .bz2 files
cp /home/saparicio/APARICIO-"$dir"/*.bz2 .
bzip2 -d *.bz2

# rename the .bdg to .bedgraph
	for i in `ls *.bdg`
		do
		in=$i
		out=`echo $i | sed 's/bdg/bedgraph/'`
		mv $in $out
		done

exit;
