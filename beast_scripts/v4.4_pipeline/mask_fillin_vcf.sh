#!/bin/sh

# Script to take positions and grep only those vcfs of interest 
# ie those that we want to validate by targeted deep sequencing
# We already only filtered based on interesting mutations (see feeder script vcf_zipper_index.sh)

source ~dyap/.bashrc

# Directory
	dir="/home/dyap/dyap_temp/vcf"

# vcfs

	TL9="hTertL9_SA039_selected.vcf" ## gains at slower rate

	TL2="hTertL2_SA040_selected.vcf"	## used for ss Seq (Chr20 gainer)
#	TL2="TertL2m_SA040.vcf" 		## remove one position <20 to prevent error

	H36="HCT116_SA036_selected.vcf" ## WT (4 lanes)
	H37="HCT116_SA037_selected.vcf" ## BRCA2-/- 18 (5 lanes)
	H38="HCT116_SA038_selected.vcf" ## BRCA2 -/- 46 (5 lanes)

cd $dir


# compare shared (common for HCT116 WT and hTERTL2)
echo "Checking for common positions...."
echo "================================="
vcf-isec -f -n +2 -a $TL2.gz $H36.gz | bgzip -c > hTertL2-HCT116-shared.vcf.gz
tabix -p vcf hTertL2-HCT116-shared.vcf.gz

# present in only HTERTL2
echo "Checking for hTERTL2 unique positions"
echo "====================================="
vcf-isec -f -c $TL2.gz $H36.gz > Unique-$TL2
vcf-isec -f -c $TL2.gz $H36.gz | bgzip -c > Unique-$TL2.gz
tabix -p vcf Unique-$TL2.gz

# present in only HCT116 WT
echo "Checking for HCT116 WT unique positions"
echo "======================================="
vcf-isec -f -c $H36.gz $TL2.gz > Unique-$H36
vcf-isec -f -c $H36.gz $TL2.gz | bgzip -c > Unique-$H36.gz
tabix -p vcf Unique-$H36.gz

# Fill in with SNV masked RefSeq for Design Space (FS=)
echo "Filling in primer design space...."
fill-fs -r /share/scratch/pipeline_temp/GenomeAnalysisPipeline/reference/GRCh37-lite.fa -c 300 -l 300 hTertL2-HCT116-shared.vcf.gz > L2-H-shared-masked.vcf 
fill-fs -r /share/scratch/pipeline_temp/GenomeAnalysisPipeline/reference/GRCh37-lite.fa -c 300 -l 300 Unique-$H36.gz > SA036_HCT116WT-masked.vcf 
fill-fs -r /share/scratch/pipeline_temp/GenomeAnalysisPipeline/reference/GRCh37-lite.fa -c 300 -l 300 Unique-$TL2.gz > SA040-HTERTL2-masked.vcf 
