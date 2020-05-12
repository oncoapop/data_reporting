#!/bin/sh

path="/home/dyap/bin/bedtools2/bin"
cmd="fastaFromBed"

# location of custom genome in fasta format
inpath="/home/dyap/dyap_temp/genomes"
#inpath="/share/lustre/projects/chip-seq/APARICIO-211/reference/"
ref="hg19_rDNA.fa"
# bzip2 -cd  $inpath"/"$ref".bz2" | tar xvf -

for expt in 216 217 218 219
	do

	# bedpath
	bedpath="/home/dyap/Projects/ChIPseq/APARICIO-"$expt"_Galaxy"
	file1="drug-treat_IgG"
	file2="drug-treat_input"

	cd $bedpath
	$path"/"$cmd -fi $inpath"/"$ref -bed $file1.bed -fo $file1.fa
	$path"/"$cmd -fi $inpath"/"$ref -bed $file2.bed -fo $file2.fa

	echo $inpath"/"$ref
	echo $file1.bed
	echo $file1.fa

	done

echo "Done!" 
exit;

# Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta> 

# Options: 
#  	      -fi     Input FASTA file
#  	      -bed    BED/GFF/VCF file of ranges to extract from -fi
#  	      -fo     Output file (can be FASTA or TAB-delimited)
#  	      -name   Use the name field for the FASTA header
#  	      -split  given BED12 fmt., extract and concatenate the sequencesfrom the BED "blocks" (e.g., exons)
#  	      -tab    Write output in TAB delimited format.
#                - Default is FASTA format.

#  	      -s      Force strandedness. If the feature occupies the antisense,
#                strand, the sequence will be reverse complemented.
#                - By default, strand information is ignored.

#  	      -fullHeader     Use full fasta header.
#                - By default, only the word before the first space or tab is used.
