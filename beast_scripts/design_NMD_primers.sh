#!/bin/sh

# This script takes the fasta format of the NMD transcript
# A list of selected transcripts as
# /home/dyap/Projects/eIF4A3/Plots/selected_transcripts.csv
# Generated manually using the R-script
# Plot_NMD-compare.R
# transcripts which increase in T-595 but not T-598

dir="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation"
outdir="/home/dyap/Projects/eIF4A3"

isofile=$outdir"/isoforms.fa"
cp $dir"/isoforms.fa" $isofile

file=$outdir"/Plots/selected_transcripts.csv"
outfile=$outdir"/design.csv"

for i in `cat $file | awk -F"," '{print $3}' | tr -d '"' | tail -n +2`
	do
	gene=`grep "$i" $file | awk -F"," '{print $2}' | tr -d '"'`
	trans=$i
	pre=`grep -A1 "$i" $isofile | grep -A1 "premRNA" | tail -n +2`
	cod=`grep -A1 "$i" $isofile | grep -A1 "coding" | tail -n +2`
	len=`echo $cod | wc -c`

	echo $gene,$trans,$cod,$len
	echo $gene"_"$trans,$cod,$len >> $outfile
	echo "+++++++++++++++++++"

	done


