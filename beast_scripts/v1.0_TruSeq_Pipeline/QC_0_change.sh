#!/bin/sh

# This script changes the samplesheet with the MSR allocated sample IDs
# which is in the different order as the Old Sample sheet

wdir="/home/dyap/Projects/ctDNA/AT094/variants"
cd $wdir

miseqdir="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094"

samplesheet=$miseqdir"/SampleSheet.csv"
sampleIDs="SampleSheetNew.txt"
output="NewSampleSheet"

rm -f $output
rm -f $output.tmp

# Each sample is represented by TWO libraries (A and B)
# if the structure is different DO NOT use this script

# Find the line number of the header and get only the samplelines from the samplesheet
line=`grep -n "Sample_ID" $samplesheet | awk -F":" '{print $1}'`
start=`echo "$line+1" | bc`  
echo $samplesheet,$start
cat $samplesheet | tail -n +$start > "SampleLines.csv"

# Get the real sample IDs from the demultiplexed run
# always removing S0 which is undetermined
ls $miseqdir"/Data/Intensities/BaseCalls" | grep fastq | awk -F"_" '{print $1"_"$2}' | sort -tS -k2n -u | tail -n +2 > $sampleIDs

echo "##########################################################"

# Get the unique SampleNames from the New Samplesheet
for i in `cat $sampleIDs | awk -F"_" '{print $1}'`
	do
	echo $i
	old1=`grep -m1 "$i" $samplesheet | awk '{ sub("\r$", ""); print }'`
	new1=`grep -m1 "$i" $sampleIDs`
	old2=`grep -m2 "$i" $samplesheet  | awk '{ sub("\r$", ""); print }' | tail -n +2 `
	new2=`grep -m2 "$i" $sampleIDs | tail -n +2`

	echo $old1","$new1 
	echo $old2","$new2 


	echo $old1","$new1 >> $output.tmp
	echo $old2"," $new2 >> $output.tmp

	cat $output.tmp

	done

echo "##########################################################"
