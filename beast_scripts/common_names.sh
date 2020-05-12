#!/bin/sh

wd="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation"
cd $wd

line1="HCT116"
line2="Hela"
drug1="T_595"
drug2="T_202"
design="highest"

L1=`ls | grep $line1 | grep $drug1 | grep $design | grep expression.res`
L2=`ls | grep $line1 | grep $drug2 | grep $design | grep expression.res`
file="Common_"$L1"-"$L2

outfile="/home/dyap/Projects/eIF4A3/"$file
rm -f $outfile

for i in `cat $L1 | awk -F" " '{print $1}' | tail -n +2`
	do
	echo $i
	grep "$i" $L2 | awk -F" " '{print $1}' | tail -n +2 >> $outfile
	done

