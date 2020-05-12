#!/bin/sh

# New script to compare the original VCFs
# 
# Script to copy vcf into directory to be read

export LC_ALL=C

# input directory
#vcfsource="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/"
source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/VariantCallingLogs"

# output directory
dest="/home/dyap/Projects/ctDNA/AT094/variants"

# Find somatic mutations
# vcf present on the FFPE but not in the Genome (G -WT, buffy coat)
cd $dest

# make the samplesheet as MiSeq Reporter has rearranged the samples IDs
cp $source/*[0-9].vcf $dest/

ls *[0-9].vcf | awk -F"." '{print $1}' > SampleSheetNew.txt

rm -f Buffy*
rm -f FFPE*
rm -f TUMOURonly_*
rm -f WTonly_*
rm -f germline_*
rm -f F
rm -f G
rm -f *.tmp
rm -f *.key

# Get the list of unique sample IDs in the samplesheet (unique and sorted)
cat SampleSheetNew.txt | awk -F"_" '{print $1}' | sed s'/[A-Z]*//' | sed '/^$/d' | sort -u > sampleIDs

for file in `cat sampleIDs`
	do
#	echo $file
	F1=`ls *.vcf | grep -m1 $file"_"`
	F2=`ls *.vcf | grep -m2 $file"_" | tail -n +2`
	F3=`ls *.vcf | grep -m3 $file"_" | tail -n +3`
	F4=`ls *.vcf | grep -m4 $file"_" | tail -n +4`
	echo $F1,$F2,$F3,$F4

# Get the list of all the variant position per sample
# Need to sort files in the same way for comm to work properly
	grep -v "#" $F1 | awk -F"\t" '{print $1","$2","$3","$4","$5}' | sort -u > $F1.key
	grep -v "#" $F2 | awk -F"\t" '{print $1","$2","$3","$4","$5}' | sort -u > $F2.key
	grep -v "#" $F3 | awk -F"\t" '{print $1","$2","$3","$4","$5}' | sort -u > $F3.key
	grep -v "#" $F4 | awk -F"\t" '{print $1","$2","$3","$4","$5}' | sort -u > $F4.key

# File the common keys between FFPE and Buffy coat samples

N1=`echo $F1 | awk -F"." '{print $1}'`
N2=`echo $F2 | awk -F"." '{print $1}'`
N3=`echo $F3 | awk -F"." '{print $1}'`
N4=`echo $F4 | awk -F"." '{print $1}'`

echo $N1, $N2, $N3, $N4 

# Comm fileA fileB
# col1=unique to fileA
# col2=unique to fileB
# col3=common to both fileA & fileB
# -12 (suppress col 1 & 2) ie only print out common names in both files

# Comparing and output the position on both strands / libraries which are common
comm -12 $F1.key $F2.key > "FFPE_"$N1"_"$N2
comm -12 $F3.key $F4.key > "Buffy_"$N3"_"$N4

N5=`echo $N1 | awk -F"_" '{print $1}'`
N6=`echo $N3 | awk -F"_" '{print $1}'`

# Comparing and output the positions present in both samples from same patient (ie germline)
comm -12 "FFPE_"$N1"_"$N2 "Buffy_"$N3"_"$N4 > "germline_"$N5"_"$N6

# Comparing and output the positions present in only Buffy sample from same patient (ie WT only)
comm -13 "FFPE_"$N1"_"$N2 "Buffy_"$N3"_"$N4 > "WTonly_"$N5"_"$N6

# Comparing and output the positions present in only FFPE sample from same patient (ie tumour only)
comm -23 "FFPE_"$N1"_"$N2 "Buffy_"$N3"_"$N4 > "TUMOURonly_"$N5"_"$N6

	done

