#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# Cut and paste the positions in format 
# chrn:xxx-xxx
# into a text file created in nano (to prevent MAC format problems)
# naming convention

name="SA029"

# $name_pos.txt
# ( This file name is used for later scripts, so keep format and store in positions/SNV directory)
# change input for GetSeqSS2.R
# output
# $name_positions.csv

export $name


# calls the script which takes the output of GetSeqSS.R 
# and parses it into primer3 input format
# input file is $name"_positions.csv"
# ~/Scripts/fasta2primer3.sh

# calls the primer3 script for multiplex conditions (single cell - Jas' project)
# ~/Scripts/Mplex-primer3.sh

# calls the script to generate list of primers as well as html output to do quick check
~/Scripts/primer_summary.sh

# things to do
# Add forward adaptor to forward primers : ACACTGACGACATGGTTCTACA
# Add reverse adaptor to reverse primers : TACGGTAGCAGAGACTTGGTCT

exit;
