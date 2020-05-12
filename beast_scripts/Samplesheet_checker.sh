#!/bin/sh
# Script to convert lookup barcodes sequences from FLD list
# Also checks the Data section of MiSeq samplesheets 
# Version 1.1
# by Damian Yap, Molecular Oncology, BCCRC
# on beast.cluster.bccrc.ca
# 3 Jul 2013

# This is the path and name of the file that we want to identify the barcodes
# and check sample number and plate mappings
############################
qpath="/home/dyap/dyap_temp"
input=$qpath/"Samplesheet_input.csv"
###########################

# DO NOT CHANGE THIS PATH AND FILE NAME
# Path and name of file that contains all the barcode mappings
inpath="/home/dyap/Projects/MiSeq_Data/MiSeq_QC"
source=$inpath/"FLD_barcodes.txt"

# Output file
##################################
out=$qpath/"Checked_Samplesheet.csv"
rm -f out
#################################


query=$qpath"/Samplesheet.tmp"
rm $qpath/*.tmp

# Read sample sheet and check it
# Convert into from MAC to UNIX format
awk '{ gsub("\r", "\n"); print $0;}' $input > $query.tmp 

clear

grep -B40 "Sample_ID" $query.tmp > $input".header"
diff $input".header" $query".tmp" | grep ">" | sed 's/^>\ //' > $query
 
echo "Checking for duplicate sample IDs - they must be unique"
echo "If any Sample_ID appears below, it is duplicated"
cat $query | awk -F, '{print $1}' | sort | uniq -d
cat $query | awk -F, '{print $1}' | sort -u > $qpath"/SID_temp.tmp"

read ans

rm $out
touch $out
echo "Processing Samplesheet..."

# Read the barcodes from one file and match them with source mapping file
for i in `cat $query | awk -F, '{print $1}'`
	do
		echo $i
		SampleID=`grep -w ^$i $query | awk -F, '{print $1}'`
		Name=`grep -w ^$i $query | awk -F, '{print $2}'`
		PLATE=`grep -w ^$i $query | awk -F, '{print $3}'`
		WELL=`grep -w ^$i $query | awk -F, '{print $4}'`
		Project=`grep -w ^$i $query | awk -F, '{print $5}'`
		Manifest=`grep -w ^$i $query | awk -F, '{print $6}'`
		Folder=`grep -w ^$i $query | awk -F, '{print $7}'`
		FLD_ID=`grep -w ^$i $query | awk -F, '{print $8}'`
		Des=`grep -w ^$i $query | awk -F, '{print $9}'`

		 if [[ $FLD_ID != "" ]];
			then {
				BC=`grep $FLD_ID $source | awk -F"\t" '{print $2}'`
				PLATE=`grep $FLD_ID $source | awk -F"\t" '{print $4}'`
				WELL=`grep $FLD_ID $source | awk -F"\t" '{print $3}'`
				}
			else
			
			 	if [[ $PLATE != "" && $WELL != "" ]];
					then {
					BC=`egrep ^.{22}.*$PLATE $source | grep -w $WELL | awk -F"\t" '{print $2}'`
					FLD_ID=`egrep ^.{22}.*$PLATE $source | grep -w $WELL | awk -F"\t" '{print $1}'`
					}
				
				fi
			fi
		echo $SampleID","$Name","$PLATE","$WELL","$Project","$Manifest","$Folder","$FLD_ID","$BC","$Des >> $out".tmp"
	done

echo "Checking for duplicate Barcodes - they must be unique"
echo "If any Barcode appears below, it is duplicated"
cat $out".tmp" | awk -F, '{print $9}' | sort | uniq -d

read ans

echo "Generating new Samplesheet..."
cat $input".header" > $out
echo "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Project,Manifest,GenomeFolder,Index_ID,Index,Description"  >> $out
cat $out".tmp" >> $out

echo $out" has been generated in "$qpath

rm $qpath/*.tmp

exit ;

