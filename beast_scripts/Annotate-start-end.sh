#!/bin/sh

# This Script was written by Damian Yap (Jul 2013)
# This scipt is called by the R-script "Genome_Range_overlap.R"
# the output of the scripts are file1 and file2
# The positions therein are then automatically annotated by this script

path="/home/dyap/Projects/PPP2R2A"
samples="Con PPP2R2A EXP"
ver="hg18"

cd $path

for i in $samples
	do

	file1=$path"/"$i"-Annotate-start.csv"
	file2=$path"/"$i"-Annotate-end.csv"

	# Process the annotation file and runs ANNOVAR
	cat $file1  | sed 's/^.*"chr//'| tr -d '"' | tr "," " " > $file1".tmp"
	perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver $ver $file1.tmp ~/bin/ANNOVAR/annovar/humandb/ 
	cat $file1".tmp.variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $file1"_anno.txt" 

	cat $file2  | sed 's/^.*"chr//'| tr -d '"' | tr "," " " > $file2".tmp"                                    
	perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver $ver $file2.tmp ~/bin/ANNOVAR/annovar/humandb/
	cat $file2".tmp.variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $file2"_anno.txt"

	done

rm -f *.tmp
rm -f *.tmp.*


exit;
