#!/bin/sh

# This Script was written by Damian Yap (Aug 2013)
# Modified to generate rerun of failed primers

# Primer3 2.3.5 is installed on pleiades.myseqdept.org
# and needs a specific input
# This script runs primer3
# on .txt files in the format (record sep "=")

# Need to change this for each file
# The name is exported from fasta2primer3.sh script
# if not, uncomment & input name here
# Project="TOV"
# type="SNV2"
# name=$Project"-"$type

echo "Name of project: "$name
echo "If this is correct, press ENTER or else Ctrl-C to exit..."

# Project Directory
dir="/home/dyap/Projects/PrimerDesign/"$Project

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

tmp="/home/dyap/dyap_temp"

if [ -e $tmp ];
        then
        echo "tmp directory is ok"
        else
        mkdir $tmp
        chmod og-wrx $tmp
fi

# Source and Output directories where Barcoded files stored
sourcedir=$posdir
outdir=$tmp
settings="/home/dyap/Projects/PrimerDesign/settings/"$p3set"_pipeline_p3_settings.txt"
raw=$posdir"/"$name"-hg19.csv"
uniq=$posdir"/"$name"_uniqpos.txt"

# Part of the pipeline, use default output of R-script GetRanges_TOV_SNV.R
in1=$p3dir"/"$name"_p3_design.txt"
primerlist=$outdir"/"$name"_primerlist.txt"
outfile=$outdir"/"$name"_p3_output"
infile=$outdir"/"$name"_primer3_input.txt"


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
if [ -f  $in1 ]
        then echo $in1
        else echo $in1 is not a valid file
             exit;
fi

echo "Generating primer3 input file..."

# tail -n +2 removes the header
cat $in1 | tail -n +2 | awk -F"," '{print $4}' | tr -d '"' | sort -u > $uniq
for i in `cat $uniq`
	do
 	id=`grep -m1 $i $in1 | awk -F"," '{print $2}' | tr -d '"'`
 	chr=`grep -m1 $i $in1 | awk -F"," '{print $3}' | tr -d '"'`
 	pos=`grep -m1 $i $in1 | awk -F"," '{print $4}'`
 	seq=`grep -m1 $i $in1 | awk -F"," '{print $8}' | tr -d '"' |  tr  -c 'ATCGatcg' 'N'`

	len=`echo $seq | wc -c`
	width=20
	sta=`echo "($len-1)/2 - ($width/2)" | bc`

	echo "SEQUENCE_ID="$id >> $infile
	echo "SEQUENCE_TEMPLATE="$seq >> $infile
	echo "P3_COMMENT="$chr"_"$pos >> $infile
	echo "SEQUENCE_TARGET="$sta","$width >> $infile
#        echo "PRIMER_PRODUCT_OPT_SIZE=1" >> $infile
        echo "PRIMER_GC_CLAMP=2" >> $infile
        echo "PRIMER_MAX_SIZE=25" >> $infile
	echo "=" >> $infile

	done 

echo "File to output " $outfile
echo "Input file (with full path): "$infile
echo Output to this directory $outdir

cd $outdir
primer3="/usr/local/bin/primer3_core"
$primer3 -output=$outfile -p3_settings_file=$settings < $infile

echo "Done. Press ENTER to continue..."

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


echo "Press Enter to show..."

grep -B13 "PAIR_NUM_RETURNED=0" $outfile 
grep -B13 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

echo "####################################################"

###################################
# reiteration module goes in here #
###################################

# Part of the redesign pipeline
rename=$tmp"/"$name"_p3_redesign.txt"

#cat $outfile.failed | sed 's/PRIMER_PRODUCT_SIZE_RANGE=260-280\ 250-300/PRIMER_PRODUCT_SIZE_RANGE=240-340\ 
#220-380\ 200-420/'| sed 's/PRIMER_MAX_SIZE=25/PRIMER_MAX_SIZE=32/' | sed 
#'s/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | grep -A19 "PRIMER_SEQUENCE_ID=" | sed 's/--/=/' > $rename

cat $outfile".failed" | sed 's/PRIMER_MAX_SIZE=25/PRIMER_MAX_SIZE=34/' | sed 's/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | grep -A7 "SEQUENCE_ID=" | sed 's/--/=/' > $rename

echo "=" >> $rename

reoutfile=$tmp"/"$name"_p3_output2"
$primer3 -output=$reoutfile -p3_settings_file=$settings  < $rename

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $rename

echo "Number of outputs in " $reoutfile
grep "^=" -c $reoutfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $reoutfile

grep -B13 "PAIR_NUM_RETURNED=0" $reoutfile > $reoutfile.failed

######################################

cat $outfile > $outfile".txt"
cat $reoutfile >> $outfile".txt"


export dir=$dir
export tmp=$tmp
export Project=$Project
export type=$type
export name=$name

echo "####################################################"
echo "Generating in-silico PCR input file"

~/Scripts/v1.5_pipeline/v1.5_p3out2isPCRin.sh

echo "####################################################"
echo "Checking amplicons/primers using in-silico PCR"

~/Scripts/v1.5_pipeline/v1.5_primer-check-isPCR.sh

export dir=$dir
export tmp=$tmp
export Project=$Project
export sample=$sample
export raw=$raw
export name=$name
export reads=$Miseq

echo "####################################################"
echo "Generating order file..."

~/Scripts/v1.5_pipeline/v1.5_primer_order.sh

echo "Files transferred to html path for viewing..."

exit;

