#!/bin/sh

# This Script was written by Damian Yap (May 2013) v1.0
# Updated Jul 2013 for v2.0 for ITH project

# The input is an all-in-one file with different numbers of SNV/indel 
# per sample of ITH

# GetRanges.R now can process short indels (on the same chromosome <40bp)
# https://github.com/AparicioLab/Ovarian_ITH_Primer_design/blob/master/GetRange-ITH-SNV.R
# as well as SNVs
# Output of the file is in ITH/primer3 folder $name_design.txt

# Version 2.0 Updates
# input ranges (<50bp) on same chromosome rather than just single positions
# outputs an all in one file
# "S/No","ID","Chr","Start","End","Indel" or "SNV","Context","Design"

# Common ID 
# Sample_chrn_startpos
# ie SA457_chr1_1234455

############
name="ITH"
###########

echo "Please confirm the name of the sample that this script will process: " $name
echo "If incorrect, please exit by pressing ctrl-c"
read ans

export name=$name

# calls the script which takes the output of GetRanges.R 
# and parses it into primer3 input format
# input file is $name"_design.txt"
# This script is specific for ITH conditions (V2.0)
echo "Running fasta2primer3ITH.sh..."
#~/Scripts/fasta2primer3ITH.sh

#export name=$name

# calls the primer3 script for single plex conditions (ITH)
echo "Running primer3ITH.sh..."
#~/Scripts/primer3ITH.sh

#export name=$name

# Process the annotation file and runs ANNOVAR
anodir="/home/dyap/Projects/ITH/Annotate"
echo "Annotating ITH positions..."
cat $anodir"/"$name"_Annotate.csv"  | awk -F, '{print $2" "$3" "$4" "$5" "$6}' | sed 's/^.*"chr//'| tr -d '"' > $anodir"/"$name
perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $anodir"/"$name ~/bin/ANNOVAR/annovar/humandb/ 
cat $anodir"/"$name".variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $anodir"/"$name"_anno.txt" 

export name=$name

echo $name_anno.txt

# calls the script to generate list of primers as well as html output to do quick check
echo "Generating summary..."
~/Scripts/primer_summaryITH.sh

# Add forward adaptor to forward primers : ACACTGACGACATGGTTCTACA
# Add reverse adaptor to reverse primers : TACGGTAGCAGAGACTTGGTCT
echo "Adding FLD adaptors...."
~/Scripts/add-adapt.sh

# export to Godel (needs password)
export godelpath="/var/www/html/workflow/primer3check/"
export exportname="/home/dyap/Projects/ITH/primer3/ITH_p3_check.html"

echo "Transferring to Godel...."
~/Scripts/put-on-Godel.sh

echo "Done."
exit;
