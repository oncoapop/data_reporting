#!/bin/sh

# Script to copy vcf into directory to be read
# input directory
source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment"

# output directory
dest="/home/dyap/Projects/ctDNA/AT094"

#cp $source/*[0-9].vcf $dest/
 
# Find somatic mutations
# vcf present on the FFPE but not in the Genome (G -WT, buffy coat)

cd $dest

# Get all the sample IDs

for i in `ls *[0-9].vcf  | awk -F. '{print $1}' | grep F | sed 's/F//'`
	do
	echo $i
	pos=`grep "^chr" G$i.vcf | awk -F"\t" '{print $2}'`	

	#echo $pos

	# germline mutations
	grep "$pos" F$i.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_germline"

	# somatic mutations
	grep -v "$pos" F$i.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_somatic"

	echo "###################################################"
	done
	


