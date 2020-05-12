#!/bin/sh

# script to read information from VCFs

# written by Damian Yap (Jun 2013)
# Ubuntu 12.04

# Source directory
indir="/home/dyap/Projects/Single_Cell/VBA0038-run3-ASROCK_MiSeq" 

# Destination directory
outdir="/home/dyap/Projects/Single_Cell/variants" 

cd $indir
for i in `ls VBA*.*vcf`
	do
	{
	outfile=$outdir"/"$i"_run3sum.txt"
	echo $outfile
	echo "Chr_pos,ref,alt,dpsnpid,fliter,depth,genotype,vfreq,gene_id,descp" > $outfile
	# The only unique value is the variant position 
	grep ^chr $i | awk -F"\t" '{print $2}' > $i".pattern"

	for j in `cat $i".pattern"`
		 do
		{
	chr=`grep $j $i | awk -F"\t" '{print $1}'`
	pos=`grep $j $i | awk -F"\t" '{print $2}'`
	snpid=`grep $j $i | awk -F"\t" '{print $3}'`
	ref=`grep $j $i | awk -F"\t" '{print $4}'`
	alt=`grep $j $i | awk -F"\t" '{print $5}'`
	qual=`grep $j $i | awk -F"\t" '{print $6}'`
	filt=`grep $j $i | awk -F"\t" '{print $7}'`
	ac=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $1}'`
	af=`grep $j $i  | awk -F"\t" '{print $8}' | awk -F";" '{print $2}'`
	an=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $3}'`
	dp=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $4}'`
	dq=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $5}'`
	ti=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $6}'`
	gi=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $7}' | awk -F"," '{print $1}'`
	fc=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $8}' | awk -F"," '{print $1}'`
	anno=`grep $j $i | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' `
	geno=`grep $j $i | awk -F"\t" '{print $10}' | awk -F":" '{print $1}'`
	wtdp=`grep $j $i | awk -F"\t" '{print $10}' | awk -F":" '{print $2}' | awk -F"," '{print $1}'`
	mtdp=`grep $j $i | awk -F"\t" '{print $10}' | awk -F":" '{print $2}' | awk -F"," '{print $2}'`

	depth=`grep $j $i | awk -F"\t" '{print $10}' | awk -F":" '{print $3}'`
	vfreq=`grep $j $i | awk -F"\t" '{print $10}' | awk -F":" '{print $6}'`

	echo $chr"_"$pos","$ref","$alt","$snpid","$filt","$dp","$geno","$vfreq","$gi","$anno >> $outfile
 		}
		done
	rm -f $i".pattern"
	echo done.

	}

	done

exit;
