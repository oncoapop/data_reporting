#!/bin/bash

# v3.3 using modified v3.1 
# for splciing signature primer design

# primer3output in temp directory
# p3outfile=$tmp"/"$name"_primer3_output.txt"

# Script to get the top hits of isPCR
# Generate an order file for those primers - done
# Generate Miseq Manifest - done
# Generate HTML view file - done
# Output FAIL Positions and info - done
# scp that to godel.cluster.bccrc.ca - not done


if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ]; then
 	echo "usage:  "$0" <Project> <Sample> <miseq>"
     	exit 0
fi

# Project Name (Directory)
Project=$1

# This is the sample name
name=$2

# Reads of a run on MiSeq
reads=$3

echo "##############################################################################################"
echo "Primer3_output must be saved in "$Project"/primer3/"$name"_primer3_output.txt"
echo "##############################################################################################"

#read ans

# Project Directory
dir="/home/dyap/Projects/PrimerDesign/"$Project

echo "-------------------------------------------------" >> $dir"/primer_design.log"
date >> $dir"/primer_design.log"

if [ -d "$dir" ]; then

	# positions
	posdir=$dir"/positions"
	if [ -d "$posdir" ]; then
		echo "Position directory ok"
	else
		mkdir $posdir
	fi

	# Annotations
	anndir=$dir"/Annotate"
	if [ -d "$anndir" ]; then
		echo "Annotation directory ok"
	else
		mkdir $anndir
	fi

	# primer3 output
	p3dir=$dir"/primer3"
	if [ -d "$p3dir" ]; then
		echo "Primer3 directory ok"
	else
		mkdir $p3dir
		echo "Script requirements not met"
		echo "Exiting..."
		exit;
	fi

fi
# Tmp file
tmp="/home/dyap/dyap_temp"

# Module to generate positions from primer3_output
p3outfile=$p3dir"/"$name"_primer3_output.txt"

# Module for annotation
annofile=$anndir"/"$name"_Annotate.csv"
annofile2=$anndir"/"$name"_anno.txt"

# De novo generation of positions
file=$name"_p3_positions.txt"
posfile=$posdir"/"$file

# Skips the GetSeq.R script
pos2file=$posdir"/"$name"_positions.txt"
designfile=$p3dir"/"$name"_p3_design.txt"

rm -f $posfile
rm -f $pos2file

echo "Generation of pos file"
for i in `cat $p3outfile | grep "PRIMER_SEQUENCE_ID="`
	do
	sample="splice1"
	chr=`echo $i  | awk -F= '{print $2}'`
	sta="99"
	end="999"
	mseq=`grep -A1 $i $p3outfile | grep SEQUENCE_TEMPLATE= | awk -F= '{print $2}'`
	cxt=`grep -A2 $i $p3outfile | grep P3_COMMENT= | awk -F= '{print $2}'`
	cpos1=`grep -A3 $i $p3outfile | grep P3_COMMENT2= | awk -F= '{print $2}' | awk -F"-" '{print $1}'`
	cpos2=`grep -A3 $i $p3outfile | grep P3_COMMENT2= | awk -F= '{print $2}' | awk -F"-" '{print $2}'`

	echo $chr","$chr":%"$cpos1"%"$cpos2"%,"$mseq >> $posfile
	echo "NN,"$chr","$sta","$chr"_%"$cpos1"%"$cops2"%" >> $pos2file
	echo $chr"_"$cpos1"-"$cpos2","$chr","$cpos1","$cpos2",WT,"$cxt","$mseq >> $designfile
	echo -ne '#'
	done

echo "positions regenerated for p3 output file."

# Generate "_isPCR-input"
# Perform in-silico PCR (using command-line isPCR)

export dir=$dir
export name=$name
export tmp=$tmp

date >> $dir"/primer_design.log"

~/Scripts/v3.3_pipeline/v3.3_inSilicoPcr2.sh

date >> $dir"/primer_design.log"

export Project=$Project
export name=$name
export file=$file
export reads=$reads
export dir=$dir

~/Scripts/v3.3_pipeline/v3.3_primer_summary2.sh

date >> $dir"/primer_design.log"
echo "===========================================================" >> $dir"/primer_design.log"

exit;


