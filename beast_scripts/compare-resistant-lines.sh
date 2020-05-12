#!/bin/sh

basedir="/share/lustre/projects/takeda_EIF4A3/pipeline_inbox/TAK-82/Data"

for sample in SA793 SA794 SA795
	do
	echo $sample

	wd=$basedir"/"$sample"/illumina_wgss"

	echo $wd
#	file=`ls $wd | grep COSMIC_v64.annotations.vcf`
	file=`ls $wd | grep dbSNP_v137.annotations.vcf`
#	file=`ls $wd | grep snvs.stats.genes.txt`

	sg="/home/dyap/Projects/EIF4A3_paper/comparisons/SG-genes.txt"

#	for i in `cat $sg`
#		do
#		echo $i
#	echo $file
#	cat $wd"/"$file | grep \|HRAS\| | wc -l
#	cat $wd"/"$file | grep \|HRAS\| | grep INDEL
#	cat $wd"/"$file | grep \|UPF | wc -l
#	cat $wd"/"$file | grep $i

	head -n2 $wd"/"$file
	cat $wd"/"$file | grep EIF4A3 | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'

#		done
	echo "###############"

	done

