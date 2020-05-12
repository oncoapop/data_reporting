#!/bin/bash

# v3 using Adi's primer3 pipeline

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
Atmp="/home/dyap/dyap_temp"

# Module to generate positions from primer3_output
p3outfile=$p3dir"/"$name"_primer3_output.txt"

# Module for annotation
annofile=$anndir"/"$name"_Annotate.csv"
annofile2=$anndir"/"$name"_anno.txt"

# De novo generation of positions
file=$name"_p3_positions.txt"
posfile=$posdir"/"$file

rm -f $posfile

echo "Generation of pos file"
for i in `cat $p3outfile | grep "PRIMER_SEQUENCE_ID="`
	do
	sample=`echo $i | awk -F= '{print $2}' | awk -F_ '{print $1}'`
	chr=`echo $i  | awk -F= '{print $2}' | awk -F_ '{print $2}' | tr -d "chr"`
	snv=`echo $i  | awk -F= '{print $2}' | awk -F_ '{print $3}'`
	mseq=`grep -A1 $i $p3outfile | grep SEQUENCE_TEMPLATE= | awk -F= '{print $2}'`

	echo $sample"_"$chr"_"$snv","$chr":"$snv"-"$snv","$mseq >> $posfile
	echo -ne '#'
	done

echo "positions regenerated for p3 output file."

# Generate "_isPCR-input"
# Perform in-silico PCR (using command-line isPCR)

export dir=$dir
export name=$name
export tmp=$tmp

date >> $dir"/primer_design.log"

~/Scripts/v3.2_pipeline/v3.2_inSilicoPcr.sh


# R script
export Project=$Project
export name=$name
export file=$file

date >> $dir"/primer_design.log"

Rscript ~/Scripts/v3.2_pipeline/v3.2_GetSeq.R --no-save --no-restore --args $Project/$name/$file

# ANNOVAR
if [ -f $annofile ]
        then
        echo "Reading from file: "$annofile
        tail -n +2 $annofile | awk -F, '{print $2" "$3" "$4" "$5" "$6" "$7}' | tr -d '"' | sed 's/chr//' > $annofile2
        # Annotation of all SNV positions
        annotate="/home/dyap/INSTALL/annovar/annotate_variation.pl"
        build="hg19"
        # dir of databases (trailing slash needed)
        dbdir="/home/dyap/INSTALL/annovar/humandb/"

                # Annotate the all the SNV in hg19
                cd $annodir

		while :;do echo -n .;sleep 1;done &
                perl $annotate -buildver $build -outfile $anndir"/"$name $annofile2 $dbdir
		kill $!; trap 'kill $!' SIGTERM

                echo "Annotation completed."
        else
        echo $annofile": File not found. Exiting...."
        exit 1;
fi

export Project=$Project
export name=$name
export file=$file
export reads=$reads
export dir=$dir

~/Scripts/v3.2_pipeline/v3.2_primer_summary.sh

date >> $dir"/primer_design.log"
echo "===========================================================" >> $dir"/primer_design.log"

exit;


