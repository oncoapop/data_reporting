#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# To add specific adaptors to primer3 generated primers for fluidigm barcoding system of MiSeq
# takes any tsv file name<tab>left_primer<tab>right_primer
# and adds adaptors to each end respectively

dir="/home/dyap/Projects/PrimerDesign/manual"
clear 
cd $dir

echo "Adding adaptors to sequences..."

sourcefile=$dir"/splice-signature-test.txt"

orderleft="Splice_test_left.csv"
orderright="Splice_test_right.csv"

rm -f $orderleft
rm -f $orderright

# Illlumina Adaptors (5'->3')
fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

# Forward adaptor for Fluidigm
# fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
# ra="TACGGTAGCAGAGACTTGGTCT"


ID=`cat $sourcefile | awk -F"\t" '{print $1}'`

test1=`cat $sourcefile | awk -F"\t" '{print $1}' | wc -l`
test2=`cat $sourcefile | awk -F"\t" '{print $1}' | sort -u | wc -l`

if [[ $test1 == $test2 ]];

	then echo "Cluster_IDs are unique"

	else    {
		echo "Cluster IDs are NOT unique...";
		exit;
		}
fi
		for i in $ID

		do
		label=`grep $i $sourcefile | awk -F"\t" '{print $1}'`
		left0=`grep $i $sourcefile | awk -F"\t" '{print $2}'`
		right0=`grep $i $sourcefile | awk -F"\t" '{print $3}'`

#		echo $left0
#		echo $right0

			if [[ $left0 != "" ]];
			
				then
					echo $label","$fa$left0 >> $orderleft	
			fi

			if [[ $right0 != "" ]];
			
				then
					echo $label","$ra$right0 >> $orderright
			fi

		
		done

echo "done."
exit;
