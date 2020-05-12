#!/bin/sh

path="/home/dyap/bin/bedtools2/bin"
cmd="fastaFromBed"

# location of custom genome in fasta format
inpath="/home/dyap/dyap_temp/genomes"
#inpath="/share/lustre/projects/chip-seq/APARICIO-211/reference/"
ref="hg19_rDNA.fa"
# bzip2 -cd  $inpath"/"$ref".bz2" | tar xvf -

for expt in 217
	do

	echo "======================================================"
	# bedpath
	echo $expt
	echo "---------"
	bedpath="/home/dyap/Projects/ChIPseq/APARICIO-"$expt"_Galaxy"
	cd $bedpath

	file1="drug-treat_IgG"
	file2="drug-treat_input"

#	file1="U2OS-gH2AXIPvsIgG-CX5461_10-7"
#	file2="U2OS-gH2AXIPvsIgG-nodrug"

# Modify depending on drug conc
# If there is a header to the file, it has to be removed.
#	cat diff-IgG.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 && $9 == 1 )  print $1"\t"$2"\t"$3 }' > $file1.bed
#	cat diff-input.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 && $9 == 1 )  print $1"\t"$2"\t"$3 }' > $file2.bed
	cat diff-IgG.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 )  print $1"\t"$2"\t"$3 }' > $file1.bed
	cat diff-input.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 )  print $1"\t"$2"\t"$3 }' > $file2.bed
#       cat diff-IgG.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 )  print $1"\t"$2"\t"$3 }' > $file1.bed
#       cat diff-input.bed | tail -n +2 | awk -F"\t" '{ if ($6 == 0 && $7 == 1 )  print $1"\t"$2"\t"$3 }' > $file2.bed

	echo "Converting from Bed to Fasta..."
	$path"/"$cmd -fi $inpath"/"$ref -bed $file1.bed -fo $file1.fas
	$path"/"$cmd -fi $inpath"/"$ref -bed $file2.bed -fo $file2.fas

	echo $inpath"/"$ref
	echo $file1.bed
	echo $file1.fas

	echo "Motif Discovery by MEME..."

# Need to filter fasta seq <8 bp
	perl /home/dyap/Scripts/remove_small.pl 8 $file1.fas > $file1-fil.fa
	perl /home/dyap/Scripts/remove_small.pl 8 $file2.fas > $file2-fil.fa

	meme -maxsize 10000 -maxsites 5 -dna -nmotifs 10 -nsites 5 $file1-fil.fa -oc $file1-MEME-IgG-drug
	meme -maxsize 10000 -maxsites 5 -dna -nmotifs 10 -nsites 5 $file2-fil.fa -oc $file2-MEME2-IgG-nodrug
	echo "======================================================"


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
