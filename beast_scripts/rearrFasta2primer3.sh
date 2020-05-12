#!/bin/sh

# This Script was written by Damian Yap (Jun 2013)

# Primer3 is installed on beast
# and needs a specific input
# This script takes the output of python / unix scripts 
# as part of the STRUCTURAL REARRANGEMENT WORKFLOW
# and parses it into primer3 input format

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/"

# Source and Output directories where Barcoded files stored
sourcedir=$dir
outdir="destruct_results/"
outfile=primer3_input
file=$dir$outdir$outfile

echo Output file = $file

# This files needs $name which is exported from the calling scipt
# IF this is run independently, please uncomment and specify string here
# name="SA029"

name="SVs"
# get the name of source file
source=$name"_selected.txt"


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

# Version 1.0 (Experimentally validated 96% by qPCR melt temp)
# The breakpoint is denoted by the "[" and "]" encoded by the python script

	echo "PRIMER_OPT_SIZE=20";
	echo "PRIMER_MIN_SIZE=17";
	echo "PRIMER_MAX_SIZE=27";
	echo "PRIMER_NUM_NS_ACCEPTED=0";
	echo "PRIMER_PRODUCT_SIZE_RANGE=140-170"
	echo "PRIMER_FILE_FLAG=0";
	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
	echo "PRIMER_MIN_TM=58.0";
	echo "PRIMER_OPT_TM=60.0";
	echo "PRIMER_MAX_TM=62.0";	
	echo "PRIMER_GC_CLAMP=2";
	echo "PRIMER_PAIR_NUM_RETURNED=1";
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
