#!/bin/sh

for i in `ls /share/lustre/archive/BOB/ctDNA/OUTPUT_tumor/RUN/F*_S*_targeted_sequencing_pipeline_single/outputs/results/realignment/*.bam`

	do
	echo $i
	cp "$i" /home/dyap/dyap_temp/bams/
	done


cd /home/dyap/dyap_temp/bams


for j in `ls *.bam`
	do

	echo "Indexing: "$j        

	samtools index "$j"

	done

