#!/bin/sh

# Script to map MISO events to the gene name

###############################

misodir="/home/dyap/Projects/EIF4A3_paper/comparisons/MISO"
pattern="/home/dyap/Projects/EIF4A3_paper/comparisons/*.annotated"

output1=$misodir"/miso_si48_genes"
tempfile1="/home/dyap/temp/tempfile_miso48"
rm -fr $tempfile1

for i in `cat $misodir/miso_48h_0.1_10 | sort -u | awk -F" " '{print $1}'`
        do
        grep -m1 -w -R "$i" $pattern | head -1 | awk -F"\t" '{print $21}' >> $tempfile1   
        grep -m1 -w -R "$i" $pattern | head -1 | awk -F"\t" '{print $21}' 

        done

cat $tempfile1 | sort -u > $output1

tempfile2="/home/dyap/temp/tempfile_miso_highT202_hela"
rm -fr $tempfile2

output2=$misodir"/miso_Hela_202high_genes"
for i in `cat $misodir/miso_Hela_202.highdose | sort -u | awk -F" " '{print $1}'`
        do
        grep -m1 -w -R "$i" $pattern | head -1 | awk -F"\t" '{print $21}' >> $tempfile2   
        grep -m1 -w -R "$i" $pattern | head -1 | awk -F"\t" '{print $21}' 
 
        done

cat $tempfile2 | sort -u > $output2

