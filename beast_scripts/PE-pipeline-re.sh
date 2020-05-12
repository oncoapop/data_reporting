#!/bin/sh

# This Script was written by Damian Yap (May 2013)

# This is the second round after modification of conditions
# edit /share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/primer3_input.txt
# based on the /share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/SA500_p3_output.failed


name="SA493"

echo "Please confirm the name of the sample to be reiterated: " $name
read ans

export name=$name

# calls the primer3 script for single plex conditions (Tumour Het - Peter's project)
~/Scripts/primer3default.sh

echo $name
export $name

echo $name_anno.txt

# calls the script to generate list of primers as well as html output to do quick check
~/Scripts/primer_summaryPE.sh


# things to do
# Add forward adaptor to forward primers : ACACTGACGACATGGTTCTACA
# Add reverse adaptor to reverse primers : TACGGTAGCAGAGACTTGGTCT

exit;
