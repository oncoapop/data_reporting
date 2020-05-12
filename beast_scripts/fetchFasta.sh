#!/bin/sh

## script used to extract fasta sequence from list of identifiers
##      in the refGene table

## given an accession list of identifiers known to be in the refGene table
## fetch the transcript coordinates from the refGene table
## score is zero for all items
cat accessionList.txt | while read A
do
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
        -e 'select chrom,txStart,txEnd,name,0,strand
                from refGene where name ="'${A}'";' hg18
done > accessionList.bed 2> sqlErrors.txt
## This does not capture any sqlErrors when the name is not found ?

## separate out the reverse and forward strand items, reformat to
## match twoBitToFa format chrN:start-stop, and with name included
awk -F'\t' '/\+$/ {printf "%s:%d-%d\t%s\n", $1,$2,$3,$4}' accessionList.bed \
        > forwardStrand.bed
awk -F'\t' '/-$/ {printf "%s:%d-%d\t%s\n", $1,$2,$3,$4}' accessionList.bed \
        > reverseStrand.bed

## fetch each sequence from the 2bit file, insert name in fasta output
## forward strand items
cat forwardStrand.bed | while read N
do
    POSITION=`echo "${N}" | cut -f1`
    NAME=`echo "${N}" | cut -f2`
    twoBitToFa hg18.2bit:${POSITION} stdout | sed -e "s/^>/>${NAME} /"
done > forward.fa

## reverse strand items
cat reverseStrand.bed | while read N
do
    POSITION=`echo "${N}" | cut -f1`
    NAME=`echo "${N}" | cut -f2`
    twoBitToFa hg18.2bit:${POSITION} stdout | sed -e "s/^>/>${NAME} /"
done > reverse.fa

## and finally, reverse compliment the reverse strand sequence
faRc -keepName -keepCase reverse.fa rev_compliment.fa
