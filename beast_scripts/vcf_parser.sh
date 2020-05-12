#!/bin/sh

# Script to take positions from vcfs of passing criteria 

# Directory
	dir="/home/dyap/dyap_temp/vcf"

##################################
##################################
###   Change these variables    ##
# vcf
#	vcf="/share/lustre/projects/epigenetic/museq_SingleSample/SA039_A24723_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
#	vcf="/share/lustre/projects/epigenetic/museq_SingleSample/SA040_A24724_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
	vcf="/home/dyap/dyap_temp/vcf/HCT116-masked.vcf"
#	vcf="/home/dyap/dyap_temp/vcf/hTertL2-masked.vcf"
#	vcf="/home/dyap/dyap_temp/vcf/L2-H-shared-masked.vcf"

# output file 
	file="SA036_HCT116.vcf"
#	file="SA040_HTERTL2.vcf"
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

##################################
###   Change these variables    ##
# Change the SAXXX in this line if you change the name above
echo -e "SA_chr_pos"'\t'"REF"'\t'"ALT"'\t'"Prob"'\t'"Geno"'\t'"Maskseq" > $out.tmp 
vcf-query $out -f 'SA036\_chr%CHROM\_%POS\t%REF\t%ALT\t%INFO/PR\t%INFO/GT\t%INFO/FS\n' >> $out.tmp
##################################

# This are the passing conditions, please change as necessary
# PR (ie $4) must be > 0.9
# GT (ie $5) must be 0/1 ie het

# This picks out those that PASS the conditions
awk -F"\t" '{ if ($4 > 0.9 && $5=="0/1") print $0}' $out.tmp | tail -n +1 > $out.txt

bname=`basename $out.txt`
	cp $out.txt "/home/dyap/Projects/PrimerDesign/cellmix/"$bname


