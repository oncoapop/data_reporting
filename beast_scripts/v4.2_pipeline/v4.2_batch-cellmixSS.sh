#!/bin/sh

# Name of Project
Project="cellmix"
# MiSeq reads
reads=200
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"

# Input the names of the samples, they are processed in series
# To do - rewrite to make them run in parallel

# 28 Nov 2014 (initial run)
#inputsamples="SA036+40_L2-HCT-selected.vcf.txt"
#inputsamples="SA036+40_L2-HCT-selected-exon.vcf.txt"

# 3 Dec 2014 (after new significant vcf filter)
inputsamples="SA036+40_L2-HCT-selected.vcf.txt"

export samples=$inputsamples
~/Scripts/v4.2_pipeline/v4.2_parse2primer3.sh 230-250 220-260  

# each sample name (space separated, no quotes, punctuations)
for sample in $samples
 

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

# copies the final version of the primer design process to the primer3 project directory
	cp $inpath"/"$sample"_primer3_output.txt.cat" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v4.2_pipeline/v4.2_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v4.2 ." >> ~/Projects/PrimerDesign/jobs.log

	done

