#!/bin/sh

# Script to take positions and grep only those vcfs of interest 
# This was modified in Jun 2015 to check context of single region chr:start-end

############################
# ENTER CHR: START-END HERE
#
chr=19
start=55115752
end=55115771
#
###########################


source ~dyap/.bashrc

# Directory
	dir="/home/dyap/dyap_temp/vcf"

# vcfs
	vcf_TertL9="/share/lustre/projects/epigenetic/museq_SingleSample/SA039_A24723_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
	vcf_TertL2="/share/lustre/projects/epigenetic/museq_SingleSample/SA040_A24724_hg19_RR0911_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
	vcf_HCT_36="/share/lustre/projects/cell_lines/museq/hct116_single_sample/SA036_A46777_hg19_RR1027_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
	vcf_HCT_37="/share/lustre/projects/cell_lines/museq/hct116_single_sample/SA037_A45342_hg19_RR1027_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"
	vcf_HCT_38="/share/lustre/projects/cell_lines/museq/hct116_single_sample/SA038_A45343_hg19_RR1027_MutationSeqSS.annotSnpEff.flagDBsnp.flag1000gen.vcf"


# output file 
	outL9=$dir"/hTertL9_SA039"
	outL2=$dir"/hTertL2_SA040"
	outWT=$dir"/HCT116_SA036"
	outBR18=$dir"/HCT116_SA037"
	outBR46=$dir"/HCT116_SA038"

# Zip and index
# At this step we want to filter with significant positions in the VCFs
# grep NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED

	cd $dir

	echo "Selecting significant variants from "$vcf_TertL9
#	cp $vcf_TertL9 $outL9.vcf
#	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outL9".vcf" > $outL9"_selected.vcf"
	cat hTertL9_SA039.vcf | awk -v chr="$chr" -v start="$start" -v end="$end" -F"\t" '{ if ($1 == chr && $2 > start && $2 < end )  print $0 }' > $outL9"_gRNA_selected.vcf"

#	echo "Selecting significant variants from "$vcf_TertL2
#	cp $vcf_TertL2 $outL2.vcf
#	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outL2".vcf" > $outL2"_selected.vcf"

#	echo "Selecting significant variants from "$vcf_HCT_36
#	cp $vcf_HCT_36 $outWT.vcf
#	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outWT".vcf" > $outWT"_selected.vcf"

#	echo "Selecting significant variants from "$vcf_HCT_37
#	cp $vcf_HCT_37 $outBR18.vcf
#	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outBR18".vcf" > $outBR18"_selected.vcf"

#	echo "Selecting significant variants from "$vcf_HCT_38
#	cp $vcf_HCT_38 $outBR46.vcf
#	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outBR46".vcf" > $outBR46"_selected.vcf"

cd $dir
# We already only filtered based on interesting mutations (see feeder script vcf_zipper_index.sh)

source ~dyap/.bashrc

# Directory
	dir="/home/dyap/dyap_temp/vcf"

for i in `ls *_gRNA_selected.vcf`
	do
	echo "Zipping "$i
	bgzip -c $i > $i.gz
	echo "Indexing "$i
	tabix -p vcf $i.gz
	done

# exit;
# If there are possible entries, then comment the exit and run the following script

# Script to take positions and grep only those vcfs of interest 
# ie those that we want to validate by targeted deep sequencing

# vcfs

#	TL9="hTertL9_SA039_selected.vcf" 	## gains at slower rate
	TL9r="hTertL9_SA039_gRNA_selected.vcf" 	## test gRNA against sequenced genome 

#	TL2="hTertL2_SA040_selected.vcf"	## used for ss Seq (Chr20 gainer)
#	TL2="TertL2m_SA040.vcf" 		## remove one position <20 to prevent error

#	H36="HCT116_SA036_selected.vcf" ## WT (4 lanes)
#	H37="HCT116_SA037_selected.vcf" ## BRCA2-/- 18 (5 lanes)
#	H38="HCT116_SA038_selected.vcf" ## BRCA2 -/- 46 (5 lanes)

cd $dir


# Fill in with SNV masked RefSeq for Design Space
echo "Filling in sequence region with any non-reference sequences...."
# fill-fs -r /share/scratch/pipeline_temp/GenomeAnalysisPipeline/reference/GRCh37-lite.fa -c 300 -l 300 $TL9r".gz" > L9_gRNA-masked.vcf 
