#!/bin/sh

# This Script was written by Damian Yap (May 2013)

# Cut and paste the positions in format 
# chrn:xxx-xxx
# into a text file created in nano (to prevent MAC format problems)
# naming convention

# use alphabets for subsequent set NOT SA999-1 (not valid identifier!)
name="SA494"

# $name_pos.txt
# ( This file name is used for later scripts, so keep format and store in positions/SNV directory)
# change input for GetSeq.R
# output
# $name_positions.csv

echo "Please confirm the name of the sample that this script will process: " $name
echo "If incorrect, please exit by pressing ctrl-c"
read ans

export name=$name


# calls the script which takes the output of GetSeq.R 
# and parses it into primer3 input format
# input file is $name"_positions.csv"
# This script is specific for Peter's conditions (V1.0, GC=2)
~/Scripts/fasta2primer3PE.sh


# calls the primer3 script for single plex conditions (Tumour Het - Peter's project)
~/Scripts/primer3default.sh

echo $name
export $name

# Process the annotation file and runs ANNOVAR
cat /home/dyap/Projects/Tumour_Evol/SNV/$name-Anno_SNVs.csv  | sed 's/^.*"chr//'| tr -d '"' | tr "," " " > /home/dyap/Projects/Tumour_Evol/Annotate/$name
perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 /home/dyap/Projects/Tumour_Evol/Annotate/$name ~/bin/ANNOVAR/annovar/humandb/ 
cat /home/dyap/Projects/Tumour_Evol/Annotate/$name.variant_function  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > /home/dyap/Projects/Tumour_Evol/positions/SNV/$name"_anno.txt" 

echo $name_anno.txt

# calls the script to generate list of primers as well as html output to do quick check
~/Scripts/primer_summaryPE.sh

# Add forward adaptor to forward primers : ACACTGACGACATGGTTCTACA
# Add reverse adaptor to reverse primers : TACGGTAGCAGAGACTTGGTCT
~/Scripts/add-adapt.sh

# export to Godel (needs password)
export godelpath="/var/www/html/workflow/primer3check/"
export exportname="/home/dyap/Projects/Tumour_Evol/positions/SNV/ITH_p3_check.html"

~/Scripts/put-on-Godel.sh
exit;
