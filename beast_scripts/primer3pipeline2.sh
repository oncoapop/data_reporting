#!/bin/sh

# This Script was written by Damian Yap (Aug 2013)
# WSOP2013-001 version 4.0

# Primer3 is installed on beast
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Project Directory
Project="TNBC"
dir="/home/dyap/Projects/"$Project
tmp="/home/dyap/dyap_temp/TNBC"

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"


cd $dir
ls

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here 
type="indel"
name=$Project"-"$type

# Source and Output directories where Barcoded files stored
sourcedir=$posdir
outdir=$p3dir
settings=$p3dir"/primer3_settings.txt"
#raw=$posdir"/primerIn-TNBC-"$type"-fix-List-AccountGermMut.txt"
raw=$posdir"/primerIn-TNBC-"$type"-fix-List-altRefSeq-AccountGermMut.txt"


# Part of the pipeline, use default output of fasta2primer3.sh
filename=$name"_primer3_input.txt"
outfile=$tmp"/"$name"_p3_output"
infile=$tmp"/"$filename

# Final output
primerlist=$tmp"/"$name"_primerlist.txt"

# Remove old file
if [ -f $infile ]
	then
	echo "Overwritting "$infile". Press ENTER to continue or Ctrl-C to exit."
	read ans
	rm -f $infile
fi

if [ -f $primerlist ]
	then
	echo "Overwritting "$primerlist". Press ENTER to continue or Ctrl-C to exit."
	read ans
	rm -f $primerlist
fi

if [ -f  $raw ]
	then echo $raw
	else echo $raw is not a valid file
	     exit;
fi

echo "Generating primer3 input file..."

for i in `cat $raw | awk -F"\t" '{print $3}'`
	do
 	sa=`grep $i $raw | awk -F"\t" '{print $1}'`
 	chr=`grep $i $raw | awk -F"\t" '{print $2}'`
 	start=`grep $i $raw | awk -F"\t" '{print $3}'`
 	end=`grep $i $raw | awk -F"\t" '{print $4}'`
	id=$sa"_"$chr"_"$start
 	ref=`grep $i $raw | awk -F"\t" '{print $6}'`
 	alt=`grep $i $raw | awk -F"\t" '{print $7}'`
	# echo $id

	# indel $ref-$alt. If +ve, then DELETION, if -ve then INSERTION
	reflen=`echo $ref | wc -c`
	altlen=`echo $alt | wc -c`
	indel=`echo "$reflen - $altlen" | bc`


	# Note that the sequence already contains the LONGEST 
	# insert whether ref or alt


		# for INSERTIONS
		# Alt is longest, therefore add (-$len) to width
		if [[ $indel < 0 ]]; 
			then 	{
				arr="ins-"$altlen
				# echo "insertion"
				len=`echo "$altlen + 25" | bc`
				} else {
					# for deletions
					# ref is longest
					# echo "deletion"
					arr="del-"$altlen
					len=`echo "0 + 25" | bc`
					}
		fi

 	indelid=`grep $i $raw | awk -F"\t" '{print $5}'`
	myid=$id"-"$arr
 	seq=`grep $i $raw | awk -F"\t" '{print $10}'`

		if [[ $snvid != "ID" ]]; 
			then 	{
			echo "SEQUENCE_ID="$indelid >> $infile
			echo "SEQUENCE_TEMPLATE="$seq >> $infile
			echo "P3_COMMENT="$myid >> $infile
			echo "SEQUENCE_TARGET=290,"$len >> $infile

			echo "=" >> $infile	
				} else {
					continue;
					}
		fi

	done 


echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

echo "#########################################################"
echo "Running primer3...."

cd $outdir
echo "primer3_core -output="$outfile" -p3_settings_file="$settings" < "$infile
primer3_core -output=$outfile -p3_settings_file=$settings < $infile

echo "Done. Press ENTER to continue..."

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


echo "Press Enter to show..."
read ans

grep -B10 "PAIR_NUM_RETURNED=0" $outfile 
grep -B10 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

export Project=$Project
export type=$type

~/Scripts/p3out2isPCRin.sh

export Project=$Project
export type=$type

~/Scripts/primer-check-isPCR.sh

export Project=$Project
export type=$type
export raw=$raw

~/Scripts/primer_order.sh

exit;
