#!/bin/sh

# This Script was written by Damian Yap (Aug 08)
# to parse MCB validation results (Aug 08)
# modified Mar  13 to format to FASTA
# the output of indexed barcodes which is all on one line

# Working Directory
dir="/home/dyap/Projects/MiSeq_Data/130624_M00897_0044_000000000-A4D1T"

outdir="/home/dyap/Projects/MiSeq_Data/fasta"

cd $dir

#	for i in `ls *R1*.fastq`
	for i in `ls SA495-NSG_S12_L001_R1_001.fastq`
	do
# Convert fastq into fasta
	outfile=`echo $i | sed 's/_S12//' | sed 's/_L001_R1_001.fastq/.fa/'`

	echo $i into $outfile

	awk '/^@M/{gsub(/^@/,">",$1);print;getline;print}' $i > $outdir/$outfile

	echo _________
	done;

echo ----------
# cd $dir/$dir2

#	for i in `ls *R2*.fastq`
	for i in `ls SA495-NSG_S12_L001_R2_001.fastq`
	do
# Convert fastq into fasta
	outfile=`echo $i | sed 's/_S12//' | sed 's/_L001_R2_001.fastq/.fa/'`

	echo $i into $outfile

	awk '/^@M/{gsub(/^@/,">",$1);print;getline;print}' $i >> $outdir/$outfile

	echo _________
	done;
echo ----------!
exit;
