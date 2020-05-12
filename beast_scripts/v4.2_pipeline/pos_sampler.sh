#!/bin/sh

# Script to take positions and grep only those vcfs of interest
# samples EQUALLY SPACED positions over the length of chromsomes

dir="/home/dyap/dyap_temp/vcf"

##################################
##   Change these variables     ##
##				##
# Use HTert-L2 (SA040) vcf calls from
#	vcf=$dir"/SA040-HTERTL2-masked.vcf"
#	file="SA040_HTERTL2-selected.vcf"

# Use HCT116 (SA036) vcf calls from
 	vcf=$dir"/SA036_HCT116WT-masked.vcf"
 	file="SA036_HCT116-selected.vcf"

# Use vcf-isec intersected file from
# 	vcf=$dir"/L2-H-shared-masked.vcf"
# 	file="SA036+40_L2-HCT-selected.vcf"


# No of positions per chromosome required
no=20
##				##
##				##
##################################

# Working Directories
        wd=$dir"/newPipeline"
        out=$wd"/"$file
        tmp=$wd"/temp"

grep "#" $vcf > $out
date > $out.log
echo "===============" >> $out.log
echo $vcf >> $out.log
echo $no " positions samples per chromsome in "$out >> $out.log

# Get the chromosome range
	chr=`cat $vcf | cut -f 1 | sort -n | uniq | sed '/^#.*/d'`

	j=$wd"/tmpfile1"
	k=$wd"/tmpfile2"

	for i in $chr
		do

# Gets all the vcfs of the chromosome $i
		grep ^$i $vcf > $k

# Lists all the positions
		cat $k | cut -f 2 | sort -n | sed '/^#.*/d' > $j
# List last position
		last=`cat $j | tail -n 1`
# Counts all the positions
		nvar=`cat $j | wc -l`

                sam=`echo "$last/$no" | bc`
                sam2=`echo "$nvar/$no" | bc`

		echo "==============================================================="
		echo "Chr="$i," Last (bp)="$last," Bin size (bp)= "$sam
		echo "Chr="$i," #var="$nvar," Sample every nth pos= "$sam2

		echo "===============================================================" >> $out.log
		echo "Chr="$i," Last (bp)="$last," Bin size (bp)= "$sam >> $out.log
		echo "Chr="$i," #var="$nvar," Sampled every nth pos= "$sam2 >> $out.log

# Sampling a position from file every $sam2 positions
		cat $k | awk -v sam="$sam2" 'NR%sam==0' >> $out

		done

rm -f $wd"/tmp*"

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

rm -f *.tmp

