#!/bin/sh
# Script to match the SNV position with the primer sequences

##############
# Requirements
name="SA500all"
folder="/home/dyap/Projects/Tumour_Evol/SA500"

primerf="/home/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_order.txt"
clusf=$folder"/"$name"_clusters.tsv"

# preprocess the cluster file from dos to unix ( ie removing ^M char at end of line)
infile=$folder"/"$name"_cluster.txt"
# dos2unix -n $clusf $infile
cat $clusf | tr "\t" ":" > $infile

outfile=$folder"/"$name"_clusterprimer.csv"

# Write header to the output file	
echo "ID,Chr,Pos,Cluster,Annotation,Left Primer,Right Primer" > $outfile

echo "Processing..."
#############
for i in `grep ":chr" $infile`

	do 
	# Get the info from the processed cluster file ($infile)
	anno=`echo $i | awk -F: '{print $1}'`
	chr=`echo $i | awk -F: '{print $2}'`
	pos=`echo $i | awk -F: '{print $3}'`
	clus=`echo $i | awk -F: '{print $4}'| tr -d "\r\n"`

	# echo $i
	# echo $clus
	# Format the match
	matchr=`echo $chr | tr "c" "C"`
	match=$matchr"_"$pos
	# echo $match

	# Match the primers from primer file (primerf)
	left=`grep -A2 $match $primerf | grep "LEFT_0_" | sed 's/^.*=//'`
	right=`grep -A2 $match $primerf | grep "RIGHT_0_" | sed 's/^.*=//'`
	# echo LEFT $left
	# echo RIGHT $right

	# Format the output file
	echo $match","$chr","$pos","$clus","$anno","$left","$right >> $outfile
	done
echo "done."
###########

exit
