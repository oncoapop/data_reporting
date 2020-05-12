#!/bin/sh

## script used to extract fasta sequence from list of identifiers
##      in the refGene table

## given an accession list of identifiers known to be in the refGene table
## fetch the transcript coordinates from the refGene table
## score is zero for all items

cmd="twoBitToFa /home/dyap/dyap_temp/genomes/hg19.2bit"

cd /home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15/Roche_150729_custom
cat exon_regions.bed | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > accessionList.bed

## separate out the reverse and forward strand items, reformat to
## match twoBitToFa format chrN:start-stop, and with name included
awk -F'\t' '/\+1$/ {printf "%s:%d-%d\t%s\n", $1,$2,$3,$4}' accessionList.bed \
        > forwardStrand.bed
awk -F'\t' '/-1$/ {printf "%s:%d-%d\t%s\n", $1,$2,$3,$4}' accessionList.bed \
        > reverseStrand.bed

## fetch each sequence from the 2bit file, insert name in fasta output
## forward strand items
cat forwardStrand.bed | while read N
do
    POSITION=`echo "${N}" | cut -f1`
    NAME=`echo "${N}" | cut -f2`
    $cmd:${POSITION} stdout | sed -e "s/^>/>${NAME} /"
done > forward.fa

## reverse strand items
cat reverseStrand.bed | while read N
do
    POSITION=`echo "${N}" | cut -f1`
    NAME=`echo "${N}" | cut -f2`
    $cmd:${POSITION} stdout | sed -e "s/^>/>${NAME} /"
done > reverse.fa

## and finally, reverse compliment the reverse strand sequence
faRc -keepName -keepCase reverse.fa rev_compliment.fa
