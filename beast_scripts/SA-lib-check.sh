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
	GSClib=`grep -w $i $infile | awk -F"\t" '{print $18}'`

	path="/share/lustre/archive/"$i"/illumina_wgss"
	if [ -d "$path" ];
		then testlib=`ls $path`
		else testlib="NONE"
	fi 

	 if [[ "$testlib" == "$GSClib" ]]
		then echo "MATCH!"
			stat=""
	 fi

	 if [[ "$testlib" != "$GSClib" ]]
		then echo "NO MATCH!"
			stat="******"
		echo $i,$GSClib,$testlib
		echo $i,$GSClib,$testlib >> $outfile
	 fi

		echo  -e $i'\t'$GSClib'\t'$testlib'\t'$stat >> $outfile2

	

	done
