#!/bin/sh

# This Script was written by Damian Yap (May 2013) v1.0
# Updated Jul 2013 for v2.0 for ITH project
# Updated to use patient specific SNP/SNVs
# Version 2.5 for ITH Project Sep 2013

# The input is an all-in-one file with different numbers of SNV/indel 
# per sample of ITH 
# PatID  Chr:SNV-SNV	  SNV/SNP-masked Sequence +/- 400bp
# PAT1   16:89858425-89858425    GTATTCNNCAGAAGNGAGACTGTGCACACCCAAACA

# Version 2.0 Updates
# input ranges (<50bp) on same chromosome rather than just single positions
# outputs an all in one file

# Version 2.5 Updates
# Accept input from SNP masked genome as above

# Common ID 
# Sample_chrn_startpos
# ie PAT12_chr1_1234455

############
name="ITH2"
Project="ITH"
type="SNV2"
###########


echo "Please confirm the name of the sample that this script will process: " $name
echo "If incorrect, please exit by pressing ctrl-c"
read ans

export name=$name
export Project=$Project
export type=$type
 
# calls the script which takes the output of GetRanges.R 
# and parses it into primer3 input format
# input file is $name"_design.txt"
# This script is specific for ITH conditions (V2.0)

# Project specific data processing
# Run once only
# cat /share/lustre/projects/ith/germline_mutations/targets_with_masked_sequences/PAT* | tr "atcg" "ATCG" > ~/Projects/ITH/primer3/ITH2_design.txt

export name=$name

echo "Running fasta2primer3ITH.sh..."
~/Scripts/fasta2primer3ITH2.sh

export name=$name

# calls the primer3 script for single plex conditions (ITH)
echo "Running primer3ITH.sh..."
~/Scripts/primer3ITH.sh

export name=$name

# Process the annotation file and runs ANNOVAR
anodir="/home/dyap/Projects/ITH/Annotate"
echo "Annotating ITH positions..."
cat $anodir"/"$name"_Annotate.csv"  | awk -F, '{print $2" "$3" "$4" "$5" "$6}' | sed 's/^.*"chr//'| tr -d '"' > $anodir"/"$name
perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $anodir"/"$name ~/bin/ANNOVAR/annovar/humandb/ 
cat $anodir"/"$name".variant_function"  | awk -F" " '{print "Chr"$3"_"$4,$6,$1,$2}' > $anodir"/"$name"_anno.txt" 

export name=$name

echo $name_anno.txt

export name=$name

# new script to do in-silico PCR with the primers and chose the best one.
echo "Preparing format for isPCR submission and checking..."
~/Scripts/p3out2isPCRin.sh

echo "Performing in-silico PCR to check primers and amplicons..."
~/Scripts/primer-check-isPCR.sh

echo "Checking primers and amplicons for MiSeq run and generating manifest..."
~/Scripts/primer_check_MiSeq.sh

echo "Formatting for ordering of primers..."
~/Scripts/primer_orderITH.sh


# export to Godel (needs password)
#export godelpath="/var/www/html/workflow/primer3check/"
#export exportname="/home/dyap/Projects/ITH/primer3/"$name"_p3_check.html"

#echo "Transferring to Godel...."
#scp $exportname dyap@godel.cluster.bccrc.ca:$godelpath


echo "Done."
exit;
