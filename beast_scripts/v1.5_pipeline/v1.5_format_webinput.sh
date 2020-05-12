#!/bin/sh

#################################################
##  Script to convert webinput to standard output
# for hg19 sample (hg18 formated by liftover.sh)
##       Dr Damian Yap , Research Scientist
##    oncoapop@sdf.org  Version 1.0 (Sep 2013)
##################################################

# Variables imported for cgi script
# $Project (project dir that is)
# $sample group together expts if you have to split them up
# $name = each submission has to have a unique number

basedir="/home/dyap"
projdir=$basedir"/Projects/PrimerDesign/"$Project
wd=$projdir"/positions"
p3dir=$projdir"/primer3"

cd $wd
# masked position file
# Format
# Sample_ID chr1:123456-123456 sequence

# input will always be hg19 (as that is only when this script is called)
input=$wd"/"$name"-hg19"

	echo "=============================================="
	cat $input

	echo "=============================================="
	echo "Please check the format of the data input:"
	echo "Chrom,Pos1,Pos2,Sample_ID,Seq"

for i in `cat $input | awk -F"," '{print $2}'`
	do
	sample=`grep $i $input | awk -F"," '{print $1}'`  
	chr=`grep $i $input | awk -F"," '{print $2}' | awk -F":" '{print $1}'`  
	pos1=`grep $i $input | awk -F"," '{print $2}' | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`  
	pos2=`grep $i $input | awk -F"," '{print $2}' | awk -F":" '{print $2}' | awk -F"-" '{print $2}'`

	seq=`grep $i $input | awk -F"," '{print $3}'`  

	# hg19
	pos2=`echo "$pos1 + 1" | bc`
	# Note the second position is overwritten by this for SNVs

	echo "chr"$chr,$pos1,$pos2,$sample,$seq >> $input".tmp" 

	echo "chr"$chr,$pos1,$pos2,$sample,$seq 

	done

cp $input".tmp" $input".csv"
rm -f *.tmp

exit;


