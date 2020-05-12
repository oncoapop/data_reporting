#!/bin/bash

# Primer Design Phase I
# Taking input from researchers and parsing into primer3 format

###############################################################
# Admin settings

# Project
project="FL"

# DO NOT CHANGE AUTOMATED DATE SETTINGS
# Output saved by date in folders
# Format YYYYMMMDD
date=`date +%Y%b%d`

wd="/home/dyap/Projects/PrimerDesign/"$project
outdir=$wd"/primer3"

	if [ -d "$outdir" ]
		then echo "Directory "$outdir;
		else mkdir $wd"/primer3"
	fi
clear

# DO NOT CHANGE THESE COMMON SETTINGS
p3set="v4"
settings="/home/dyap/Projects/PrimerDesign/settings/"$p3set"_pipeline_p3_settings.txt"

cd $wd
echo $samples

# This is the length of each MiSeq run (300 cycle kit ie 2x151 bp)
Miseq=150

for s in $samples
        do
	outfile=$outdir"/"$s"_primer3_input.txt"
	p3outfile=$outdir"/"$s"_primer3_output.txt"
	infile=$wd"/"$s

	rm -f $outfile
	rm -f $p3outfile

	for i in `cat $infile | awk -F"\t" '{print $1}'`
		do

		ID=`grep -m1 $i $infile | sed 's/N_/N_chr/g' | awk -F"\t" '{print $1}'`
		Ref=`grep -m1 $i $infile | awk -F"\t" '{print $2}'`
		Alt=`grep -m1 $i $infile | awk -F"\t" '{print $3}'`
		Type=`grep -m1 $i $infile | awk -F"\t" '{print $4}'`
		Pos=`grep -m1 $i $infile | awk -F"\t" '{print $5}'`
		Size=`grep -m1 $i $infile | awk -F"\t" '{print $6}'`
		# clean up to make unrecognised characters and spaces N in the Seq
		Seq=`grep -m1 $i $infile | awk -F"\t" '{print $7}' | tr -c 'ATCGatcg' 'N' `

		llim=`echo "$Pos - 10" | bc`
		echo $llim, $Pos
		range=`echo $llim"-"$Pos`
		echo $range

		if [[ "$Pos" -eq "target" ]];
			then continue
			else
		# This is the context matching for viewing only
		cxtSeq=`grep -m1 $i $infile | awk -F"\t" '{print $7}' | cut -c$range`
		fi

		seqLength=`echo $Seq | wc -c`

	if [[ "$seqLength" -lt "$Miseq" ]];
		# Not enough design space, next.
		then continue;
		else
	{
	echo "SEQUENCE_ID="$ID;
	echo "SEQUENCE_TEMPLATE="$Seq;
	echo "P3_COMMENT="$cxtSeq":"REF_"$Ref":"VB_"$Alt":"$Type;

# <start>,<length>
        echo "SEQUENCE_TARGET="$llim",20";

# Specifies the optimal primer sizes
#        echo "PRIMER_OPT_SIZE=22";
#        echo "PRIMER_MIN_SIZE=18";
#        echo "PRIMER_MAX_SIZE=26";

# Specifies the number of unknown bases in any primer
        echo "PRIMER_NUM_NS_ACCEPTED=0";

# Specifies the product range (or ranges)
        echo "PRIMER_PRODUCT_SIZE_RANGE=160-175 145-190";

        echo "PRIMER_GC_CLAMP=2";
        echo "PRIMER_NUM_RETURN=5";
        echo "PRIMER_EXPLAIN_FLAG=1";
        echo "=";
        } >> $outfile
	fi
		done


echo "Running Primer3..."

cd $wd

primer3="/opt/apps/primer3/primer3-2.3.5/src/primer3_core"
$primer3 -output=$p3outfile -p3_settings_file=$settings < $outfile

echo $s" Done."
done

