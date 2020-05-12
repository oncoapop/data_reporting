#!/bin/sh

# Script to check the actual CG that matched factors to the CG that we designed for
# This will look at two levels.

# Firstly, check to see if there are exact matches (ie same breakpoints)
# Secondly, check to see if there are matches in both the conjoined gene (but maybe the breakpoints differ)

# Files
wd="/home/dyap/Projects/Takeda_T3/CG"
motifs="HCT116_hTERT_CG_high_motif_dens.tsv"
ordered="CG_primers_ordered.txt"
outfile="CG-matches.txt"
cd $wd


# Exact matches
echo "==============" > $outfile
echo "Exact matches:" >> $outfile
echo "==============" >> $outfile

for i in `cat $ordered`
	do
	grep "$i" $motifs >> $outfile
	echo $i
	grep "$i" $motifs 
	echo "======================"
	done

# Gene pair matches
echo "=================" >> $outfile
echo "Gene pair matches:" >> $outfile
echo "=================" >> $outfile

for j in `cat $ordered | sed 's/\@.*.*\:/:/' | sed 's/\@.*.$//'`
	do
	gene1=`echo $j | awk -F":" '{print $1}'`
	gene2=`echo $j | awk -F":" '{print $2}'`

	test=`grep "$gene1" $motifs | grep $gene2`
	if [[ $test != "" ]]
		then
		grep "$gene1" $ordered | grep $gene2 >> $outfile
		grep "$gene1" $motifs | grep $gene2 >> $outfile
	fi

	echo $j
	grep "$gene1" $motifs | grep $gene2
	echo "======================"
	done


