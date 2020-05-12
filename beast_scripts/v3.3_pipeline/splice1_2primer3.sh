#!/bin/sh

# This is the script to read the files from the splicing project
# input format http://wiki.mo.bccrc.ca/x/RIpC

# Sample,splice_id,seq1,seq1_winlen
# SA123,gene1.ex1@cord-gene2.int3@cord,=========[XY]-------------------,2,

temp="/home/dyap/dyap_temp"
indir="/home/dyap/Projects/PrimerDesign/manual"
outdir="/home/dyap/Projects/PrimerDesign/manual"

outfile=$outdir"/hct116_htert_primer3_output.txt"
infile=$indir"/hct116_htert_primer3_input.txt"

# This module selects only a subset of the cluster ID to work with.
# If you want to use the all the results comment out this section

	tfas=$temp"/Splice1_primer3_input.fa"

rm -f $tfas

#<<<-------SELECTED SVs--------->>>
# This is the file has we have used to short-list the ones that we want to work on
sel="/home/dyap/Projects/PrimerDesign/manual/hct116_htert_primer_pipeline_input.csv"
sltd="/home/dyap/Projects/PrimerDesign/manual/hct116_htert_primer_pipeline_input.txt"


# This is in MAC format so we have to translate that into unix format
# awk '{ gsub("\r", "\n"); print $0;}' $sel > $sltd
# echo "converted"

# To extract only the primer_id (since sample primer set for all samples used)
cat $sel | awk -F, '{print "hct116-htert,"$2","$3","$4}' | sort -u > $sltd

for i in `cat $sltd |  awk -F, '{print $2}'` 
	do

	if [[ $i =~ "splice_id" ]];
		then 	{
		echo "Skipping header..."
		continue
			 }
	fi

	echo $i


# String (position) Field
# $1	Sample
# $2	cluster_id (unique)
# $3	sequence
# $4 	window length

#	sa=`grep $i $sltd | awk -F, '{print $1}'` 
# all primers to be run on all samples 

	id=`grep $i $sltd | awk -F, '{print $2}'` 

	seq=`grep $i $sltd | awk -F, '{print $3}'` 
	win=`grep $i $sltd | awk -F, '{print $4}'` 


# Version 1.1 (with manual reiterations)
# The breakpoint is enclosed by the "[ ]"; windows length=$4 
# This counts the first instance for any other character other than A-Z which happens to be `[`
# Hence this leftdist represents the target start
	leftdist=`echo $seq | grep -P -o '(?<='$match')[A-Z]*(?='$match')' | wc -c`
	right=`echo "$leftdist + $win" | bc`

	cxt=`echo $seq | awk -F"[" '{print $2}' | awk -F"]" '{print $1}'`

	{
	echo "PRIMER_SEQUENCE_ID="$id"_"$leftdist"-"$right; 
	echo "SEQUENCE_TEMPLATE="$seq ;
	echo "P3_COMMENT="$cxt ;
	echo "SEQUENCE_TARGET="$leftdist","$win;
	echo "PRIMER_OPT_SIZE=21";
        echo "PRIMER_MIN_SIZE=18";
        echo "PRIMER_MAX_SIZE=26";
        echo "PRIMER_NUM_NS_ACCEPTED=0";
        echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190"
        echo "PRIMER_FILE_FLAG=0";
        echo "PRIMER_PICK_INTERNAL_OLIGO=0";
        echo "PRIMER_MIN_TM=58.0";
        echo "PRIMER_OPT_TM=60.0";
        echo "PRIMER_MAX_TM=62.0";
        echo "PRIMER_GC_CLAMP=1";
        echo "PRIMER_NUM_RETURN=5";
	echo "PRIMER_MAX_POLY_X=4";
        echo "PRIMER_EXPLAIN_FLAG=1";
        echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
        echo "=";
        } >> $tfas

	done

echo "Completed" 

#<<<-------SELECTED SVs--------->>>

# clean up to make unrecognised characters and spaces

cat $tfas  | tr -d "[]" | sed 's/=\ /=/' > $infile

echo "Running Primer3..."

cd $outdir

#primer3_core -format_output -output=primer3_view.txt < $infile
primer3_core -output=$outfile < $infile

echo "done."

echo "Press enter to continue..."
read ans

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


#echo "Press Enter to show..."

#read ans
#grep -B25 "PAIR_NUM_RETURNED=0" $outfile | more
grep -B25 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

echo Press Return to continue...
read ans

echo "Total failed"
grep -c "PAIR_NUM_RETURNED=0" $outfile
total=`grep -c "PAIR_NUM_RETURNED=0" $outfile`

exit;
	
