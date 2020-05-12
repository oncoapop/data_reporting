#!/bin/sh

# This script splits TNBC samples up by Sample number for submission to online primer design pipeline

base="/home/dyap/Projects"
Project="TNBC"
pos="positions"
p3="primer3"

p3dir=$base"/"$Project"/"$p3
posdir=$base"/"$Project"/"$pos

temp="/home/dyap/dyap_temp"

posfile=$posdir"/primerIn-TNBC-SNV-fix-List-AccountGermMut.txt"

cat $posfile | awk -F"\t" '{print $1}' | sort -u > $temp"/Samples.tmp"

for i in `cat $temp"/Samples.tmp"`
	do

	grep $i $posfile | awk -F"\t" '{print $5","$2":"$3"-"$3","$8}' > $posdir"/"$i
		
	done


exit;


