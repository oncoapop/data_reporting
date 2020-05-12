#!/bin/sh

#################################################
##  Script to convert -hg18 files to -hg19 files
##           using chain and liftOver 
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

liftover=$basedir"/bin/LiftOver/liftOver"
chain=$basedir"/bin/LiftOver/hg18ToHg19.over.chain"

cd $wd
# masked position file
# Format
# Sample_ID chr1:123456-123456 sequence

mpfile=$wd"/"$name"-hg18"

# Output will always be hg19
new=$wd"/"$name"-hg19"
old=$wd"/"$name"_hg18.tmp"

echo "Listing of "$mpfile
cat $mpfile
rm -f $old
rm -f $old".txt"

for i in `cat $mpfile | awk -F"," '{print $2}'`
	do
	sample=`grep $i $mpfile | awk -F"," '{print $1}'`  
	chr=`grep $i $mpfile | awk -F"," '{print $2}' | awk -F":" '{print $1}'`  
	pos=`grep $i $mpfile | awk -F"," '{print $2}' | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`  
	seq=`grep $i $mpfile | awk -F"," '{print $3}'`  

	# Do not add chr:pos info as they are hg18 and are confusing after conversion to 
	# hg19
	id=$sample
	pos2=`echo "$pos + 1" | bc`

	echo "chr"$chr" "$pos" "$pos2" "$id" "$seq >> $old

	done

echo "LIFTING OVER..."
unmap=$wd"/"$name"-unmap.txt"
$liftover $old $chain $new $unmap 

echo "==================================="

echo "These are the hg18 positions:"
cat $old

echo "These are the hg19 positions:"
cat $new

echo "==================================="

echo "These are the positions which cannot be mapped:"

cat $unmap

echo "==================================="

# This section standardizes the output to match the style:

cat $new | tr "\t" "," > $new".csv"

exit;


