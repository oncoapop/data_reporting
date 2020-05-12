#!/bin/bash

# Primer Design Phase I
# Taking input from researchers and parsing into primer3 format
# This improved version takes the input for the amplicon length

if [ "$1" = "" ] || [ "$2" = "" ] || [[ $1 != *[0-9]"-"[0-9]* ]] || [[ $2 != *[0-9]"-"[0-9]* ]] ; then
        echo "usage:  "$0" <optimal size range> <secondary size range>, where range is <123-456>"
        exit 0
fi

# Optimal size range (start-end)
size1=$1

# Seconday size  range (start-end)
size2=$2


###############################################################
# Admin settings

# Project
project="cellmix"

wd="/home/dyap/Projects/PrimerDesign/"$project
dir=$wd
outdir=$wd"/primer3"
temp="/home/dyap/dyap_temp"

	if [ -d "$outdir" ]
		then echo "Directory "$outdir;
		else mkdir $wd"/primer3"
	fi
clear

# DO NOT CHANGE THESE COMMON SETTINGS
p3set="v4"
settings="/home/dyap/Projects/PrimerDesign/settings/"$p3set"_pipeline_p3_settings.txt"

cd $wd


# This is the length of each MiSeq run (300 cycle kit ie 2x151 bp)
# exported from v4.2_batch
# reads=150

for s in $samples
        do
		echo "Generating Primer3 input file for "$s" ..."

	outfile=$temp"/"$s"_primer3_input.txt"
	p3outfile=$temp"/"$s"_primer3_output.txt"
	infile=$wd"/"$s

	rm -f $outfile
	rm -f $p3outfile
	rm -f $p3outfile.cat

	for i in `cat $infile | awk -F"\t" '{print $1}' | tail -n +2`
		do

	####################################################
	#   This parser must be changed for varying input  #
	####################################################


		ID=`grep -m1 $i $infile | awk -F"\t" '{print $1}'`
		Ref=`grep -m1 $i $infile | awk -F"\t" '{print $2}'`
		Alt=`grep -m1 $i $infile | awk -F"\t" '{print $3}'`
		Pos=301 ## by definition


		####################################################################
		# clean up to make unrecognised characters and spaces N in the Seq #
		####################################################################

#		Seq=`grep -m1 $i $infile | awk -F"\t" '{print $6}' | tr -d '[' | tr -d ']' | tr -c 'ATCGatcg' 'N' `
		Seq=`grep -m1 $i $infile | awk -F"\t" '{print $6}' | tr -d '['  | sed 's/\/.*\]//g'`


		################################################################################
		# For matching of context, we get the 5' 10 bp of the SNV or nreakpoint or end #
		################################################################################

		llim=`echo "$Pos - 5" | bc`
		rlim=`echo "$Pos + 5" | bc`
#		echo $llim, $rlim
		range=`echo $llim"-"$rlim`
#		echo $range

		if [[ "$Pos" -eq "target" ]];
			then continue
			else

		#################################################
		# This is the context matching for viewing only #
		#################################################

		cxtSeq=`echo $Seq | cut -c$range`
		fi

		seqLength=`echo $Seq | wc -c`

	if [[ "$seqLength" -lt "$Miseq" ]];
		# Not enough design space, next.
		then continue;
		else
	{
	echo "SEQUENCE_ID="$ID;
	echo "SEQUENCE_TEMPLATE="$Seq;

	###################################################################
	# The Context, ALT, REF alleles are carried through as P3_COMMENT #
	###################################################################

	echo "P3_COMMENT="$cxtSeq"@"REF_"$Ref":"VB_"$Alt":"


	#############################################
	# Primer3 Parameters specific to the expt   #
	# Note see ~/Projects/PrimerDesign/settings #
	#     for general settings for Primer3      #
	#############################################

	# <start>,<length>
	plim=`echo "$Pos - 10" | bc`
        echo "SEQUENCE_TARGET="$plim",20";

	# Specifies the optimal primer sizes
	echo "PRIMER_OPT_SIZE=22";
	echo "PRIMER_MIN_SIZE=18";
	echo "PRIMER_MAX_SIZE=26";

	# Specifies the number of unknown bases in any primer
	# For N Masking to prevent primer overlap
	# This must be uncommented 
        echo "PRIMER_NUM_NS_ACCEPTED=0";

	# Specifies the product range (or ranges)
	# Must be within range of $reads
	# dlim1=`echo "0.9 * $reads" | bc` # For automated calculation etc
        echo "PRIMER_PRODUCT_SIZE_RANGE="$size1" "$size2;

        echo "PRIMER_GC_CLAMP=2";
        echo "PRIMER_NUM_RETURN=5";
        echo "PRIMER_EXPLAIN_FLAG=1";
        echo "=";
        } >> $outfile
	fi
		done



cd $wd

################################################################################
# Primer3 on beast (note different versions on beast may yield varying output) #
################################################################################

primer3="/opt/apps/primer3/primer3-2.3.5/src/primer3_core"

echo "Running Primer3..."

$primer3 -output=$p3outfile -p3_settings_file=$settings < $outfile

echo $s" Done."

echo "########################"
echo "#   First Pass check   #"
echo "########################"

echo " "
echo "Number of sequences in input file"
cat $infile | wc -l

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep "PAIR_NUM_RETURNED=0" -c $p3outfile

cat $p3outfile |  grep -B18 "PAIR_NUM_RETURNED=0" > $p3outfile".failed"


###################################
# Reiteration module goes in here #
###################################

# Part of the redesign pipeline
rename=$temp"/"$s"_p3_redesign.txt"


# Reiteration - basic GC clamp 2->0 only
cat $p3outfile.failed | sed 's/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | grep -A9 "SEQUENCE_ID=" | sed 's/--/=/' > $rename

echo "=" >> $rename

echo "Re-running primer3..."
reoutfile=$temp"/"$s"_p3_output2"
$primer3 -output=$reoutfile -p3_settings_file=$settings  < $rename

echo "########################"
echo "#  Reiteration check   #"
echo "########################"

echo "Number of sequences in input file"
grep -c "SEQUENCE_TEMPLATE=" $rename

echo "Number of outputs in " $reoutfile
grep "^=" -c $reoutfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $reoutfile

grep -B18 "PAIR_NUM_RETURNED=0" $reoutfile > $reoutfile.failed

##############################################
# combine primer3 output from both processes #
##############################################

cat $p3outfile >> $p3outfile.cat
cat $reoutfile >> $p3outfile.cat

echo $s" completed"
echo " "
echo "############"


done

