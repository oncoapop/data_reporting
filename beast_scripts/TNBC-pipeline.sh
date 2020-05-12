#!/bin/sh

# This Script was written by Damian Yap (May 2013) v1.0
# Updated Jul 2013 for v2.0 for TNBC project
# Updated Aug 2013 V3.0 for SNP/SNV masking of germline positions

# The input is an all-in-one file with different numbers of SNV/indel 
# per sample of TNBC

# GetRanges.R now can process short indels (on the same chromosome <40bp)
# https://github.com/AparicioLab/TNBC-pipeline/blob/master/GetSeq_SNV.R
# as well as SNVs
# Output of the file is in TNBC/primer3 folder $name_design.txt

# Version 3.0 Updates
# input ranges (<50bp) on same chromosome rather than just single positions
# Already masks SNVs and SNPs Germline with N to avoid primers in those regions

# outputs an all in one file
# "Sample","Chr","Start","End","ID","startPosRegion","endPosRegion","seqRegion"

# Common ID 
# Sample_chrn_startpos
# ie SA457_chr1_1234455

############
name="indel"
###########

# Send notification of progress by email to:
addr="dyap-on-beast@damianeva.org"

echo "Please confirm the name of the sample that this script will process: " $name
echo "If incorrect, please exit by pressing ctrl-c"
read ans

export name=$name

# calls the script which takes the output of GetRanges.R 
# and parses it into primer3 input format
# input file is $name"_design.txt"
# This script is specific for TNBC/SNV/InDel conditions (V3.0)
 message="Running fasta3primer3TNBC.sh..."
 echo $message
 echo $message | mail $addr -s "Message from Script #1"

~/Scripts/fasta2primer3TNBC.sh


# calls the primer3 script for single plex conditions (TNBC)
 message="Running primer3TNBC.sh..."
 echo $message
 echo $message | mail $addr -s "Message from Script #2"

~/Scripts/primer3TNBC.sh

echo $name
export $name

# Process the annotation file and runs ANNOVAR
anodir="/home/dyap/Projects/TNBC/Annotate"
message="Annotating TNBC positions..."
echo $message
echo $message | mail $addr -s "Message from Script #3"

cat $anodir"/"$name"_Annotate.csv"  | awk -F, '{print $2" "$3" "$4" "$5" "$6}' | sed 's/^.*"chr//'| tr -d '"' > $anodir"/"$name
perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $anodir"/"$name ~/bin/ANNOVAR/annovar/humandb/ 
cat $anodir"/"$name".variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $anodir"/"$name"_anno.txt" 

export name=$name

echo $name_anno.txt

# calls the script to generate list of primers as well as html output to do quick check
message="Generating summary..."
echo $message
echo $message | mail $addr -s "Message from Script #4"

~/Scripts/primer_summaryTNBC.sh

# Add forward adaptor to forward primers : ACACTGACGACATGGTTCTACA
# Add reverse adaptor to reverse primers : TACGGTAGCAGAGACTTGGTCT
echo "Adding FLD adaptors...."
~/Scripts/add-adapt.sh

# export to Godel (needs password)
export godelpath="/var/www/html/workflow/primer3check/"
export exportname="/home/dyap/Projects/TNBC/primer3/SNV_p3_check.html"

echo "Transferring to Godel...."
#~/Scripts/put-on-Godel.sh

message="Done."
echo $message
echo $message | mail $addr -s "Message from Script #5"
echo "TNBC pipeline done" | mail $addr -s "Message from TNBC pipeline"

exit;
