#!/bin/sh

wd="/share/scratch/dyap_temp/Pipeline/pipelines/primer_design_pipeline_v1.0"

cd $wd

samples="SA532 SA535"

museqsuffix="_mutSeqCombinedInfo.txt"

selected="_selectpos"

output="selected"

for s in $samples
	do

	rm -f $s"_"$output".tmp"

	head -n3 $s$museqsuffix > $s"_"$output".txt"
	echo "#CHROM	POS	RefBase	AltBase" >> $s"_"$output".txt"
	echo "These positions are not found in the mutseq file" > $s"_notfound.txt"
	echo "#CHROM	POS" >> $s"_notfound.txt"

	for p in `cat $s$selected | awk -F: '{print $2}'`
		do
		chr=`grep -m1 "$p" $s$museqsuffix | awk -F"\t" '{print $1}'` 
		ref=`grep -m1 "$p" $s$museqsuffix | awk -F"\t" '{print $5}'` 
		alt=`grep -m1 "$p" $s$museqsuffix | awk -F"\t" '{print $6}'`

		chrom=`grep -m1 "$p" $s$selected | awk -F: '{print $1}'`

		if [[ $chr == "" || $ref == "" || $alt == "" ]]
			then
			echo $chrom":"$p >> $s"_notfound.txt"
			echo "NOT FOUND : "$chrom":"$p 

			else
			echo $chr,$p,$ref,$alt 
			echo $chr,$p,$ref,$alt >> $s"_"$output".tmp"


		fi
		 
		done

	cat $s"_"$output".tmp" | tr "," "\t" >> $s"_"$output".txt"

	done

exit;


