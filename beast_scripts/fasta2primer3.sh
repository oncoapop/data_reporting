#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# Primer3 is installed on beast
# and needs a specific input
# This script takes the output of GetSeq.R (ie csv)
# and parses it into primer3 input format

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/"

# Source and Output directories where Barcoded files stored
sourcedir=$dir
outdir="SNV/"
outfile=primer3_input
file=$dir$outdir$outfile

echo Output file = $file

# This files needs $name which is exported from the calling scipt
# IF this is run independently, please uncomment and specify string here
# name="SA029"

# get the name of source file
source=$name"_positions.csv"


cd $sourcedir
rm -f $dir$outdir/$outfile
rm -f $dir$outdir/$outfile.txt
rm -f $dir$outdir/*.tmp

# Preprocessing that file to make it all in one line (separated by ",")
cat $source | sed 's/.*chr/ID=/' | sed s'/",/_/' | tr -d '"' | sed 's/,/@SEQUENCE_TEMPLATE=/' > $file.tmp

	
	for i in `grep "ID=" $file.tmp` 
	do
	{ 
	echo $i;

# ID inputted into comment field as well
	echo "P3_COMMENT="$i| sed 's/@.*$//'| sed 's/ID=/Chr/';

# The SNV is in the middle, get position from GetSeq.R it is position 151
# <start>,<length>
	echo "SEQUENCE_TARGET=125,50";

# Specifies the optimal primer sizes
	echo "PRIMER_OPT_SIZE=22";
	echo "PRIMER_MIN_SIZE=18";
	echo "PRIMER_MAX_SIZE=27";

# Specifies the number of unknown bases in any primer
	echo "PRIMER_NUM_NS_ACCEPTED=0";

# Specifies the product range (or ranges)
	echo "PRIMER_PRODUCT_SIZE_RANGE=175-225";

# Specifies the min overlap of multiple primers
#	echo "PRIMER_MIN_THREE_PRIME_DISTANCE=0";

# Specifies where the primer pairs can be 
# left_start>,<left_length>,<right_start>,<right_length>
	echo "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,125,175,125";

# DEFAULTS FROM WEB INTERFACE
# http://primer3.wi.mit.edu

# Flag to retain 
	echo "PRIMER_FILE_FLAG=0";
	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
	echo "PRIMER_MIN_#=58.0";
	echo "PRIMER_OPT_#=60.0";
	echo "PRIMER_MAX_#=62.0";
	echo "PRIMER_GC_CLAMP=1";
	echo "PRIMER_NUM_RETURN=2";
	echo "PRIMER_EXPLAIN_FLAG=1";
	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
	echo "=";
	} >> $dir$outdir/$outfile
	done

cat $dir$outdir/$outfile | sed 's/ID=/PRIMER_SEQUENCE_ID=Chr/' | tr "@" "\n" | tr "#" "TM" > $dir$outdir/$outfile.txt

# Clean up comment for testing
rm -f $dir$outdir/$outfile
rm -f $dir$outdir/*.tmp


exit;
