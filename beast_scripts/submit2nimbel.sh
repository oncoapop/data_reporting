#!/bin/sh

#################################################
##  Script to convert output to exon bed format
##       Dr Damian Yap , Research Scientist
##    oncoapop@sdf.org  Version 1.0 (Jul 2015)
##################################################

# Format of input
# CHROMOSOME    START   STOP    NAME
#chr1    15949181        15949360        -_15_CASP9_SE_exon1
#chr1    15949181        15949360        -_14_CASP9_SE_exon1

# Format of output
#chrom start stop gene score strand type
#1 chr15 73062668 73062770 ENST00000362710 0.000000 -1 exon

pwd="/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15/Roche_150725_HG19_Splice-Sign2_RNA"
cd $pwd
input="Splice-Sign2_regions.txt"
outfile="exon_regions.bed"

echo -e "chrom\tstart\tstop\tgene\tscore\tstrand\ttype" > $outfile

for i in `cat $input | awk -F"_" '{print $2"_"$3"_"$4"_"$5}'`
	do
	echo $i
	chrom=`grep -m1 $i $input | awk -F"\t" '{print $1}'`
	start=`grep -m1 $i $input | awk -F"\t" '{print $2}'`
	end=`grep -m1 $i $input | awk -F"\t" '{print $3}'`
	gene=`grep -m1 $i $input | awk -F"\t" '{print $4}' | awk  -F"_" '{print $2"_"$3"_"$4"_"$5}' | sed 's/\r$//'`
	strand=`grep -m1 $i $input | awk -F"\t" '{print $4}' | awk  -F"_" '{print $1}'`
	
	echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\texon" 
	echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\texon" >> $outfile

	done

exit;

