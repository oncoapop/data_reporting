#!/bin/bash

# v3.3 using modified v3.1 
# for splciing run on primer design

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
	sample=$name
	id=`echo $i  | awk -F= '{print $2}'`
	mseq=`grep -A1 -m1 $i $p3outfile | grep SEQUENCE_TEMPLATE= | awk -F= '{print $2}'`
	cxt=`grep -A5 -m1 $i $p3outfile | grep P3_COMMENT= | awk -F= '{print $2}'`

	echo $sample"_"$id","$mseq >> $posfile
	echo $sample"_"$id","$mseq >> $pos2file
	echo $id","$cxt","$mseq >> $designfile
	echo -ne '#'
	done

echo "positions regenerated for p3 output file."

# Generate "_isPCR-input"
# Perform in-silico PCR (using command-line isPCR)

export dir=$dir
export name=$name
export tmp=$tmp

date >> $dir"/primer_design.log"

~/Scripts/v3.3_pipeline/v3.3_inSilicoPcrRT.sh


# R script
export Project=$Project
export name=$name
export file=$file
export dir=$dir
export name=$name
export tmp=$tmp

date >> $dir"/primer_design.log"


###########################################################################################
# Annotation module does not work since we are using RNA transcriptome rather than genome
###########################################################################################

# ANNOVAR
isout=$dir"/primer3/"$name"_isPCR-output.fa"
anno=$tmp"/"$name"_anno.tmp"
anno1=$tmp"/"$name"_anno1.tmp"
anno2=$tmp"/"$name"_anno2.tmp"
#rm -f $annofile

#	grep ">" $isout | awk -F" " '{print $1" "$2}'| tr -d ">" | sed 's/+/-/' > $anno
#	grep ">" $isout | awk -F" " '{print $1}'| tr -d ">" | sed 's/[+-]/:/' | awk -F":" '{print $1":"$2"-"$2+1}' > $anno1
#	grep ">" $isout | awk -F" " '{print $1}'| tr -d ">" | sed 's/[+-]/:/' | awk -F":" '{print $1":"$3"-"$3+1}' > $anno2

#        twoBitToFa /home/dyap/dyap_temp/genomes/hg19.2bit $anno1.tmp -seqList=$anno1 
#        twoBitToFa /home/dyap/dyap_temp/genomes/hg19.2bit $anno2.tmp -seqList=$anno2 

#	cat $anno1.tmp | tr "\n" " " | tr ">" "\n" | sort -u | sed 's/chr//' | tr "\-\+\:" " " | awk -F" " '{print $1,$2,$2,$4,"C"}'> $anno1.tmp2
#	cat $anno2.tmp | tr "\n" " " | tr ">" "\n" | sort -u | sed 's/chr//' | tr "\-\+\:" " " | awk -F" " '{print $1,$2,$2,$4,"G"}'> $anno2.tmp2

#	cat $anno1.tmp2 | tail -n +2 > $annofile
#	cat $anno2.tmp2 | tail -n +2 >> $annofile


#if [ -f $annofile ]
#        then
#        echo "Reading from file: "$annofile

#        cat $annofile | sed 's/chr//' > $annofile2
#        # Annotation of all SNV positions
#        annotate="/home/dyap/INSTALL/annovar/annotate_variation.pl"
#        build="hg19"
#        # dir of databases (trailing slash needed)
#        dbdir="/home/dyap/INSTALL/annovar/humandb/"

#                # Annotate the all the SNV in hg19
#                cd $annodir

#               while :;do echo -n .;sleep 1;done &
#                perl $annotate -buildver $build -outfile $anndir"/"$name $annofile2 $dbdir
#               kill $!; trap 'kill $!' SIGTERM

#                echo "Annotation completed."
#        else
#        echo $annofile": File not found. Exiting...."
#        exit 1;
# fi


export Project=$Project
export name=$name
export file=$file
export reads=$reads
export dir=$dir

#~/Scripts/v3.3_pipeline/v3.3_primer_summaryRT.sh
~/Scripts/v3.3_pipeline/v3.3_primer_summaryRT2.sh

date >> $dir"/primer_design.log"
echo "===========================================================" >> $dir"/primer_design.log"

exit;


