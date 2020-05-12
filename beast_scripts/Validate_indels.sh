#!/bin/sh

# script to read information from VCFs 
# and design positions from the indels


# written by Damian Yap (Jun 2013)
# Ubuntu 12.04

# Source directory
indir="/home/dyap/Projects/Tumour_Evol/Indel" 

# Destination directory
outdir="/home/dyap/Projects/Tumour_Evol/Indel" 

cd $indir
for i in `ls *.tsv`
	do
	{
	fname=`echo $i | sed 's/variantCalls-rmNI_L_Mere-orgCol-ALLsamp-ID-//'`
	outfile=$outdir"/"$fname"_sum.txt"
	echo $outfile
	echo "Sample,Chr_pos,ref,alt,ref-length,alt-length" > $outfile
	# The only unique value is the variant position 
	cat $i | awk -F"\t" '{print $2}' > $i".pattern"

	for j in `cat $i".pattern"`
		 do
		{
		chr=`grep $j $i | awk -F"\t" '{print $1}'`
		pos=`grep $j $i | awk -F"\t" '{print $2}'`
		snpid=`grep $j $i | awk -F"\t" '{print $3}'`
		ref=`grep $j $i | awk -F"\t" '{print $4}'`
		alt=`grep $j $i | awk -F"\t" '{print $5}'`
		reflen=`grep $j $i | awk -F"\t" '{print $4}' | wc -c`
		altlen=`grep $j $i | awk -F"\t" '{print $5}' | wc -c`
		sample=`grep $j $i | awk -F"\t" '{print $11}' | sed 's/-.*//'`

		echo $sample","$chr":"$pos","$ref","$alt","$reflen","$altlen >> $outfile
 		}
		done
	
	rm -f $i".pattern"
	echo done.

	}

	done

exit;
