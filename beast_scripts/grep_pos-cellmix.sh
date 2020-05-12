#!/bin/sh

# Script to grep the unique and shared positions primers from file
# Need to redesign new primers to fill in the drop out

wd="/home/dyap/Projects/Tumour_Evol/Cell_line_Mixing"

# This is the list of primers that we need to keep the primer sequence
# These primers must be checked on the Bam files to see if they overlap 
# or are adjacent to any SNP

#anno=$wd"/Annotated.csv"
anno=$wd"/Redesign.csv"

rm -f $anno
#file=$wd"/keep.txt"
file=$wd"/redesign.txt"
source=$wd"/CellMix.csv"

for i in `cat $file`
	do
	for type in hct116 htert shared
		do
		tab=`echo $i | sed s'/Chr//' | tr "_" "\t" | tr -d " "`

		checkfile=$wd"/"$type"_pos.txt"
#		echo $checkfile
		lab=`grep "$tab" $checkfile`
		label=`echo $lab | awk -F" " '{print "Chr"$1"_"$2}'`
 
		if [[ $label != "Chr_" ]];
			then 
			left=`grep $i $source | awk -F"," '{print $5}' | tr -d '"' | tr -d " "`	
			right=`grep $i $source | awk -F"," '{print $6}' | tr -d '"' | tr -d " "`
#			echo $label","$left","$right","$type >> $anno	
			echo $label","$type >> $anno	

		left=""
		right=""
		type=""
		fi
		
		done
	done

	for type in hct116 htert shared
		do
		
		count=`grep $type $anno | wc -l` 

		echo $type" = "$count

		done

