#!/bin/sh

# Name of Project
Project="cellmix3"
# MiSeq reads
reads=150
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp"

# Input the names of the samples, they are processed in series
# To do - rewrite to make them run in parallel

# 28 Nov 2014 (initial run)
#inputsamples="SA036+40_L2-HCT-selected.vcf.txt"
#inputsamples="SA036+40_L2-HCT-selected-exon.vcf.txt"

# 1 Dec 2014 (after new significant vcf filter)
# inputsamples="SA036_HCT116-selected.vcf.txt SA040_HTERTL2-selected.vcf.txt SA036+40_L2-HCT-selected.vcf.txt"

# 2 Mar 2015 (get positions preselected against for adjacent SNPs)
#/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix2/DAH54/targets_54.tsv
#/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix2/DAH55/targets_55.tsv
#/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix2/Shared/targets_Shared.tsv

# 2 Mar 2015 (get positions preselected against for adjacent SNPs)
#/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix2/HCT/targets_HCT.tsv
#/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix2/hTERT/targets_hTERT.tsv

# 2 Mar 2015 (for DAH preselected for no adjacent SNPs)
# inputsamples="DAH54 DAH55 Shared"

# 2 Mar 2015 (for preselected for no adjacent SNPs)
# inputsamples="HCT hTERT"

########## LATEST RUN #################

# 16 Mar 2015 (get positions preselected against for adjacent SNPs)
/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix3/DAH54/targets_54.tsv
/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix3/DAH55/targets_55.tsv
/home/dyap/R-patched/bin/Rscript ~/Scripts/v4.3_pipeline/v4.3_GetSeq.R --no-save --no-restore --args cellmix3/Shared/targets_Shared.tsv

inputsamples="DAH54 DAH55 Shared"


export samples=$inputsamples
~/Scripts/v4.3_pipeline/v4.3_parse2primer3.sh 150-160 140-190

# each sample name (space separated, no quotes, punctuations)
for sample in $samples
 

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

# copies the final version of the primer design process to the primer3 project directory
	cp $inpath"/"$sample"_primer3_output.txt.cat" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v4.3_pipeline/v4.3_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v4.3 ." >> ~/Projects/PrimerDesign/jobs.log

	done

