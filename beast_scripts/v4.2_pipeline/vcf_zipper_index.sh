#!/bin/sh

# Script to take positions and grep only those vcfs of interest 
# ie those that we want to validate by targeted deep sequencing

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
	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outL9".vcf" > $outL9"_selected.vcf"

	echo "Selecting significant variants from "$vcf_TertL2
#	cp $vcf_TertL2 $outL2.vcf
	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outL2".vcf" > $outL2"_selected.vcf"

	echo "Selecting significant variants from "$vcf_HCT_36
#	cp $vcf_HCT_36 $outWT.vcf
	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outWT".vcf" > $outWT"_selected.vcf"

	echo "Selecting significant variants from "$vcf_HCT_37
#	cp $vcf_HCT_37 $outBR18.vcf
	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outBR18".vcf" > $outBR18"_selected.vcf"

	echo "Selecting significant variants from "$vcf_HCT_38
#	cp $vcf_HCT_38 $outBR46.vcf
	grep "#\|NON_SYNONYMOUS_CODING\|NON_SYNONYMOUS_START\|STOP_GAINED\|STOP_LOST\|CODON_DELETION\|CODON_INSERTION\|EXON_DELETED" $outBR46".vcf" > $outBR46"_selected.vcf"

cd $dir

for i in `ls *_selected.vcf`
	do
	echo "Zipping "$i
	bgzip -c $i > $i.gz
	echo "Indexing "$i
	tabix -p vcf $i.gz
	done

exit;

# Links to mask_fillin_vcf.sh

