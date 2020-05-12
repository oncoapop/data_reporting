#!/bin/sh

# Scipt to generate the Manifest file for subset of positions for
# Xenographs samples or any single cell project given primers
# This uses the file that is used to generate the PCR primers

# Enter the name of sample to run
name="SA029"

platedir="/home/dyap/Projects/Single_Cell/positions/SNV"

# amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_amplicons.fa"
amplicon="/home/dyap/Projects/Single_Cell/positions/"$name"_positions.csv"

# Working Directory
dir="/home/dyap/Projects/Single_Cell/positions/SNV"


