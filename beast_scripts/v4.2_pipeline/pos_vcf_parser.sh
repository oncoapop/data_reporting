#!/bin/sh

# Script to take positions and grep only those vcfs of interest 
# ie those that we want to validate by targeted deep sequencing

# Directory
	dir="/home/dyap/dyap_temp/vcf"

##################################
##################################
###   Change these variables    ##
# vcf
#	vcf="/share/lustre/projects/epigenetic/museq_SingleSample/SA039_A24723_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
#	vcf="/share/lustre/projects/epigenetic/museq_SingleSample/SA040_A24724_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
#	vcf="/home/dyap/dyap_temp/vcf/HCT116-masked.vcf"
	vcf="/home/dyap/dyap_temp/vcf/hTertL2-masked.vcf"
#	vcf="/home/dyap/dyap_temp/vcf/L2-H-shared-masked.vcf"

# output file 
#	file="SA036_HCT116.vcf"
	file="SA040_HTERTL2.vcf"
#	file="SA036,SA040_HCT116-HTERTL2-shared.vcf"

# IF you change this name change the SP later

# positions
# Format
# Chr1_12345678<tab>shared

#	input="/home/dyap/Projects/PrimerDesign/cellmix/hTert"
#	input="/home/dyap/Projects/PrimerDesign/cellmix/Shared"
#	input="/home/dyap/Projects/PrimerDesign/cellmix/HCT116"
##################################
##################################

# Working Directories
	wd=$dir"/newPipeline"
	out=$wd"/"$file
	tmp=$wd"/temp"

	cp $vcf $out

#===========================================
# IF you want to do selection uncomment this
# Grepping the headers of the vcf file 
#	grep "#" $vcf > $out

# Chr,pos
#	chrpos=`cat $input | awk -F"\t" '{print $1}'`
#	echo $chrpos

#for match in $chrpos
#	do
#	echo $match
#	pos=`echo $match | awk -F"_" '{print $2}'`
#	chr=`echo $match | awk -F"_" '{print $1}' | sed 's/Chr//'`
#	echo $chr,$pos

#	grep $pos $vcf > $tmp

#	cat $tmp

#	grep $chr $tmp >> $out
#	done

# cat $out

#echo -e "SP_:chr_pos:REF_:ALT_:ST_cell-line"'\t'"Prob"'\t'"Geno"'\t'"Maskseq" > $out.txt 
#vcf-query $out -f 'chr%CHROM\_%POS:REF\_%REF:ALT\_%ALT:ST\_cell-line\t%INFO/PR\t%INFO/GT\t%INFO/FS\n' >> $out.txt
#
##################################
###   Change these variables    ##
# Change the SAXXX in this line if you change the name above
echo -e "SA_chr_pos"'\t'"REF"'\t'"ALT"'\t'"Prob"'\t'"Geno"'\t'"Maskseq" > $out.tmp 
vcf-query $out -f 'SA040\_chr%CHROM\_%POS\t%REF\t%ALT\t%INFO/PR\t%INFO/GT\t%INFO/FS\n' >> $out.tmp
##################################

# This are the passing conditions, please change as necessary
# PR (ie $4) must be > 0.9
# GT (ie $5) must be 0/1 ie het

# This picks out those that PASS the conditions
awk -F"\t" '{ if ($4 > 0.9 && $5=="0/1") print $0}' $out.tmp | tail -n +1 > $out.txt

#if [ $( cat check.txt | wc -l ) -gt "1" ];
#then
#	head check.txt
#	num=`cat check.txt | wc -l`
#    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++"
#    echo "Failed - all conditions not met for "$num" positions"
#    echo "****************************************************"

#else

#	echo "Showing first lines of file: "$out.txt
#	head $out.txt | column -t

#    echo "***************************"
#    echo "All conditions met - passed"
#    echo "***************************"

#fi

bname=`basename $out.txt`
	cp $out.txt "/home/dyap/Projects/PrimerDesign/cellmix/"$bname


