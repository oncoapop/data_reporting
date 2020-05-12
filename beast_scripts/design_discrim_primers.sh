#!/bin/sh

# This script takes the fasta format two transcripts (one consensus and one NMD)
# Designs discrimininatory PCR primers to be able to identify

# As well as common primers which amplify and then use sequencing to distinguish between them
# provided that they length difference is not too large

dir="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation"
outdir="/home/dyap/Projects/eIF4A3"

isofile=$outdir"/isoforms.fa"
cp $dir"/isoforms.fa" $isofile

file=$dir"/HCT116_T_595_normalized_highest_expression.res"

for i in `cat $file | awk -F" " '{print $3}' | tail -n +2 | sort -u`
	do
	gene_id=`grep "$i" $file | awk -F" " '{print $2}'`
	gene=`grep "$i" $file | awk -F" " '{print $1}'`
	NMD=`grep "$i" $file | awk -F" " '{print $3}'`
	CON=`grep "$i" $file | awk -F" " '{print $4}'`

	N=$NMD"_coding"
#	C=$CON"_coding"
	C=$CON
	
	NMD_Seq=`grep -A+1 "$N" $isofile | tail -n +2` 
	CON_Seq=`grep -A+1 '$C' $isofile | tail -n +2` 

#	echo $gene,$gene_id,$NMD,$NMD_Seq,$CON,$CON_Seq
#	echo $NMD_Seq
	len=`echo $NMD_Seq | wc -c`

	echo $gene,$gene_id,$NMD,$len
	done 
