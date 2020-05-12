#!/bin/sh

# Script to take positions and grep only those vcfs of interest 
# ie those that we ALREADY HAVE primers we want to keep (just check they do not overlap SNP)
# We already only filtered based on interesting mutations (see feeder script vcf_zipper_index.sh)

source ~dyap/.bashrc

# Directory
	dir="/home/dyap/dyap_temp/vcf"

cd $dir

# Get the positions from /home/dyap/Projects/PrimerDesign/cellmix/Check_Shared
chk="/home/dyap/Projects/PrimerDesign/cellmix/Check_Shared"

	# Get the VCF headers!
	grep "#" HCT116_SA036.vcf > HCT-shared-check.vcf
	grep "#" hTertL2_SA040.vcf > hTertL2-shared-check.vcf
	grep "#" hTertL2-HCT116-shared.vcf > HCT116_hTertL2-shared-check.vcf

for i in `cat $chk`
	do
	chr=`echo $i | awk -F"_" '{print $1}'`
	pos=`echo $i | awk -F"_" '{print $2}'`
	echo "==========================="
	echo $chr ":" $pos
	HCT=`awk -F"\t" -v pos="$pos" chr="$chr"'{ if ($1==chr && $2==pos) print $0 }' HCT116_SA036.vcf`
	echo "HCT116 WT :"$HCT
	echo $HCT >> HCT-shared-check.vcf
	L2=`awk -F"\t" -v pos="$pos" chr="$chr"'{ if ($1==chr && $2==pos) print $0 }' hTertL2_SA040.vcf`
	echo "hTert-L2 :" $L2
	echo $L2  >> hTertL2-shared-check.vcf

	refH=`echo $HCT | awk -F" " '{print $4}'`
	altH=`echo $HCT | awk -F" " '{print $5}'`

	refL=`echo $L2 | awk -F" " '{print $4}'`
	altL=`echo $L2 | awk -F" " '{print $5}'`

	echo $refH,$refL
	echo $altH,$altL

		if [[ $refH=$refL ]] && [[ $altH==$altL ]]; then
			echo $L2 >> HCT116_hTertL2-shared-check.vcf
		else
			echo "No match."
		fi

	echo "Done."

	done

echo "Positions in file: " 
cat $chk | wc -l

echo "Positions matched in shared vcf file: "
grep -c "PR="  hTertL2-shared-check.vcf 
grep -c "PR="  HCT-shared-check.vcf
grep -c "PR="  HCT116_hTertL2-shared-check.vcf

cat HCT-shared-check.vcf | bgzip -c > HCT-shared-check.vcf.gz
cat  hTertL2-shared-check.vcf | bgzip -c >  hTertL2-shared-check.vcf.gz
sed '/^$/d'  HCT116_hTertL2-shared-check.vcf | tr " " "\t" | bgzip -c >  HCT116_hTertL2-shared-check.vcf.gz

# Fill in with SNV masked RefSeq for Design Space (FS=)
echo "Filling in primer design space...."
fill-fs -r /share/scratch/pipeline_temp/GenomeAnalysisPipeline/reference/GRCh37-lite.fa -c 300 -l 300 HCT116_hTertL2-shared-check.vcf.gz > L2-H-shared-masked_check.vcf
