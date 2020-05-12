#!/bin/sh

# Script to map VAST events to the gene name

vastdir="/home/dyap/Projects/EIF4A3_paper/comparisons/VAST"

output1=$vastdir"/vast_si48_genes"
tempfile1="/home/dyap/temp/tempfile_si48_vast"
echo "Saving to "$output1" ..."
rm -fr $tempfile1

for i in `cat $vastdir/vast_48h_0.1_10 | sort -u | awk -F" " '{print $1}'`
	do
	grep -m1 -w -R "$i" /share/lustre/abashash/EIF4A3/vast_tools_pipeline | head -1 | awk -F":" '{print $2}' | awk -F"\t" '{print $1}' >> $tempfile1
	grep -m1 -w -R "$i" /share/lustre/abashash/EIF4A3/vast_tools_pipeline | head -1 | awk -F":" '{print $2}' | awk -F"\t" '{print $1}'
	done

cat $tempfile1 | sort -u > $output1

echo "Saving to "$output1" ..."

tempfile2="/home/dyap/temp/tempfile_HeLaT202high_vast"
rm -fr $tempfile2

output2=$vastdir"/vast_Hela_202high_genes"
for i in `cat $vastdir/vast_Hela_202.highdose | sort -u | awk -F" " '{print $1}'`
        do
        grep -m1 -w -R "$i" /share/lustre/abashash/EIF4A3/vast_tools_pipeline | head -1 | awk -F":" '{print $2}' | awk -F"\t" '{print $1}' >> $tempfile2   
        grep -m1 -w -R "$i" /share/lustre/abashash/EIF4A3/vast_tools_pipeline | head -1 | awk -F":" '{print $2}' | awk -F"\t" '{print $1}'
        done

cat $tempfile2 | sort -u > $output2

echo "Saving to "$output2" ..."

echo "VAST...completd"

