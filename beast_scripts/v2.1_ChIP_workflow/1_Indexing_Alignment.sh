#!/bin/sh

# This script sets up the files and parameters for the ChIP-Seq workflow

# Working directories
fq="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/150306_M02348_0082_000000000-ADHR5/Data/Intensities/BaseCalls"
bm="/home/dyap/Projects/PPP2R2A/test"

# reference genome
ref="/share/lustre/reference/genomes/ucsc.hg19.fasta"
#ref="/home/dyap/bin/bwa-0.7.17/bwakit/hs38.fa"

# command - ALIGNER (BWA)
cmd="/home/dyap/bin/bwa-0.7.17/bwa"

# Temp storage
tmp=$bm"/tmp"
rm -f $tmp

# Finding the files in the Fastq directory
cd $fq
fn=`ls *fastq.gz`

# Generating the respective indexes for the fastqs
for i in $fn
	do
	echo $i >> $tmp
	$cmd aln $ref $i > $bm"/"$i"_aln.sai"
	done

# Performing the alignment using aligner
names=`cat $tmp | awk -F"_" '{print $1}' | sort | uniq`

for j in $names
	do
	fastqs=`ls $fq"/"$j*fastq.gz`
	sais=`ls $bm"/"$j*sai`	

# checking
#	echo "+++++++++++++++"
#	echo $fastqs
#	echo $sais
#	echo "________"
#	echo $names
#	echo "++++++++++++++++++++++++++++++++++"
#	echo $j

	$cmd sampe $ref $sais $fastqs > $bm/$j"_aln.sam"

	done

