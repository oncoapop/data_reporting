#!/bin/sh

# Script to counter check SA IDs are correctly mapped to GSC IDs
# source ls /share/lustre/archive/<SA ID>/illumina_wgss

pwd="/home/dyap/Projects/TNBC_Database"
infile=$pwd"/ID-lib4"
outfile=$pwd"/ID-lib4-nomatch"
outfile2=$pwd"/ID-lib4-all"

echo -e "SA_ID"'\t'"GSC_lib"'\t'"Analysed_Lib" > $outfile
echo -e "SA_ID"'\t'"GSC_lib"'\t'"Analysed_Lib" > $outfile2

SAID=`cat $infile | awk -F"\t" '{print $1}'`

for i in $SAID
	do
	echo $i

	path="/share/lustre/archive/"$i"/"
	if [ -d "$path" ];
		then testlib=`ls $path`
		else testlib="NONE"
	fi 

		echo  -e $i'\t'$testlib'\t'$stat >> $outfile2

	

	done
