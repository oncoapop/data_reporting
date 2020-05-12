#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# To take the primers for each Multiplex reaction and pool into fasta file
# submit to http://www.thermoscientificbio.com/webtools/multipleprimer/

dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/SNV/"
clear 
cd $dir

# EXPORT NAME = SAnnn
		
infile=$name"_p3_order.txt"

sourcefile=$name"_pos.txt"

outfile=$name"_muptest.txt"
outfile2=$name"_w-adaptors.txt"

echo "File for submission to http://www.thermoscientificbio.com/webtools/multipleprimer/" > $outfile
echo "File for ordering with adaptors added - please check orientation!!" > $outfile2
# Forward adaptor for Fluidigm (SS-project)
fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (SS-project) (5'->3')
ra="TACGGTAGCAGAGACTTGGTCT"

	for i in `cat $sourcefile | awk -F- '{print $1}' | tr ":" "_" | tr "c" "C"`
		do

		left0=`grep -A4 $i $infile | grep "PRIMER_LEFT_0_SEQUENCE=" | sed 's/PRIMER_LEFT_0_SEQUENCE=//'` 
		right0=`grep -A4 $i $infile | grep "PRIMER_RIGHT_0_SEQUENCE=" | sed 's/PRIMER_RIGHT_0_SEQUENCE=//'` 
		left1=`grep -A4 $i $infile | grep "PRIMER_LEFT_1_SEQUENCE=" | sed 's/PRIMER_LEFT_1_SEQUENCE=//'` 
		right1=`grep -A4 $i $infile | grep "PRIMER_RIGHT_1_SEQUENCE=" | sed 's/PRIMER_RIGHT_1_SEQUENCE=//'` 

		echo $left0
		echo $right0
		echo $left1
		echo $right1

			if [[ $left0 != "" ]];
			
				then
					echo ">"$i"_left0" >> $outfile
					echo ">"$i"_left0+adaptor" >> $outfile2
					echo $left0 >> $outfile	
					echo $fa$left0 >> $outfile2	
			fi

			if [[ $right0 != "" ]];
			
				then
					echo ">"$i"_right0" >> $outfile
					echo ">"$i"_right0+adaptor" >> $outfile2
					echo $right0 >> $outfile
					echo $ra$right0 >> $outfile2
			fi

			if [[ $left1 != "" ]];
			
				then
					echo ">"$i"_left1" >> $outfile
					echo ">"$i"_left1+adaptor" >> $outfile2
					echo $left1 >> $outfile
					echo $fa$left1 >> $outfile2
			fi

			if [[ $right1 != "" ]];
			
				then
					echo ">"$i"_right1" >> $outfile
					echo ">"$i"_right1+adaptor" >> $outfile2
					echo $right1 >> $outfile
					echo $ra$right1 >> $outfile2
			fi


		done

exit;
