#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# To add specific adaptors to primer3 generated primers for fluidigm barcoding system of MiSeq

exit;
#dir="/home/dyap/Projects/Tumour_Evol/Rearrangements/"
clear 
cd $dir

echo "Adding adaptors to sequences..."

sourcefile="destruct_valbk_primer_order.txt"

orderleft="Destruct_left.csv"
orderright="Destruct_right.csv"

rm $orderleft
rm $orderright

# Illlumina Adaptors (5'->3')
# fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
# ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

# Forward adaptor for Fluidigm
fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
ra="TACGGTAGCAGAGACTTGGTCT"


grep "chimeric" $sourcefile | awk -F"_" '{print $1}' > pattern

test1=`grep "chimeric" $sourcefile | awk -F"_" '{print $1}' | wc -l`
test2=`cat pattern | wc -l`

if [[ $test1 == $test2 ]];

	then echo "Cluster_IDs are unique"

	else    {
		echo "Cluster IDs are NOT unique...";
		exit;
		}
fi
		for i in `cat pattern`

		do
		label=`grep $i $sourcefile | awk -F"_" '{print $1"_"$2"-"$4}'`
		left0=`grep -A2 $i $sourcefile | grep "PRIMER_LEFT_0_SEQUENCE=" | sed 's/PRIMER_LEFT_0_SEQUENCE=//'` 
		right0=`grep -A2 $i $sourcefile | grep "PRIMER_RIGHT_0_SEQUENCE=" | sed 's/PRIMER_RIGHT_0_SEQUENCE=//'` 

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
