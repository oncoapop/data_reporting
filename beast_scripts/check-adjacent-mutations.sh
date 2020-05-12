#!/bin/sh

# This script is written to check to see if there are any adjacent mutations
# in the first instance there are quite a few mutations most somatic

# need to filter for germline ones

# output path

path="/home/dyap/Projects/ctDNA/AT094/validation/"
outpath="$path/mutations/"

# The summary file
input=$outpath"/summary.txt"

# The position file
pos=$path"/expt1.csv"

# The output file
outfile=$outpath"/adj_mutations.txt"

rm -f $outfile

# Get the sample name & positions from position file
patient=`cat $pos | awk -F"," '{print $1","$2","$3}' | sed 's/BOB0//'`

for i in `echo $patient`
	do

	echo "------------------------------------------------" >> $outfile

	echo $i
	sam=`echo $i | awk -F, '{print $1}'`
	chr=`echo $i | awk -F, '{print $2}'`
	mut=`echo $i | awk -F, '{print $3}'`

	lowlim=`echo "$mut-100" | bc`
	uplim=`echo "$mut+100" | bc`

	id=`cat $pos | grep $sam | grep $chr | grep $mut`
	echo $id >> $outfile
#	echo $sam 
#	echo $chr
#	echo $mut 
#	echo $lowlim 
#	echo $uplim 

	check=`grep "$sam" $input | awk -F" " '{print $2}' | sort -u`


	norm=`grep "G$sam" $input | awk -F" " '{print $2}' | sort -u`
	ffpe=`grep "F$sam" $input | awk -F" " '{print $2}' | sort -u`

				statnorm="norm"				
				statgerm="germ"				

	echo "========"

	echo $check

	for j in `echo $check`
		do

		# Check sample name and see if there are any mutations within 100 bp of the mutation

		statnorm=`echo $norm | grep $j`
		statffpe=`echo $ffpe | grep $j`

		if [[ $j > $lowlim && $j < $uplim && $j != $mut ]]
			then
				echo "Adjacent mutation : "$j >> $outfile


			if [[ $statnorm == $statffpe ]]
				then
					echo "Germline mutation : "$j >> $outfile
					statnorm="norm"				
					statgerm="germ"				
			fi

		fi
				
		done

	echo "=====" >> $outfile

	done



