#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory
dir=/home/dyap/dyap_temp/ChIPseqAnalysis2
outpath=/share/Monco/Damian/ChIP-seq

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 216 217 218 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/jobs.log ; echo "Started processing "$expt >> $dir/jobs.log

        fullout="$outpath"/$Project/EXPT-$expt
	fullin=$dir"/"$Project"/EXPT-"$expt

	if [ -d $outpath"/"$Project ]; then
			mkdir -pv $fullout
                else
                        mkdir -pv $outpath"/"$Project
			mkdir -pv $fullout
        fi

	echo "Copying bedgraph files..."
	for file in `ls $fullin/*.bedgraph`
		do 
		echo "Copying "$file"..."
#		cp "$file" "$fullout"
		echo "=================="
		done

	echo "Copying bam files..."
	for file in `ls $fullin/*.bam`
		do 
		echo "Copying "$file"..."
		cp "$file" "$fullout"
		echo "=================="
		done

	echo "Copying bai files..."
	for file in `ls $fullin/*.bai`
		do 
		echo "Copying "$file"..."
		cp "$file" "$fullout"
		echo "=================="
		done

	echo "Copying PDF files..."
	for file in `ls $fullin/*.pdf`
		do 
		echo "Copying "$file"..."
		cp "$file" "$fullout"
		echo "=================="
		done


	 echo `date` \n >> $dir/jobs.log ; echo "Finished processing "$expt >> $dir/jobs.log

	done

exit
##############
