#!/bin/sh

# Script to counter check SA IDs are correctly mapped to GSC IDs that are returned
# source ls /share/lustre/archive/<SA ID>/illumina_wgss

# In the email - you cut and paste the paths to a file and then change the filename in infile
# /projects/analysis/analysis30/A08453/HGH7FCCXY_1/A08453_HGH7FCCXY_1_gsc_library.summary

# This script will return all the matches in the file hierarchy - ie previously analysed
# Hence this is only useful for resequencing matches

pwd="/home/dyap/Projects/TNBC_Database"
infile=$pwd"/Reseq_Jan2018"
outfile=$pwd"/SAID-GSC-match"
tempfile=$pwd"/temp"
	 ls /share/lustre/archive/SA*/illumina_wgss > $tempfile

echo -e "SA_ID"'\t'"GSC_lib"'\t'"Analysed_Lib" > $outfile

GSCID=`cat $infile | awk -F"/" '{print $5}'`
n=0

for i in $GSCID
	do
	n=`echo "$n + 1" | bc`
	match=`grep -B1 "$i" $tempfile`
	SAID=`echo $match | awk -F"/" '{print $5}'`
	SAID2=`echo $match | awk -F"/" '{print $10}'`

	echo $n" : "$i" = "$SAID" ,  "$SAID2
	done
