#!/bin/sh

# This Script was written by Damian Yap (Jul 2013)

# This is the script for reiteration of primer3 conditions
# edit /home/dyap/Projects/ITH/primer3/primer3_input.txt
# based on the /home/dyap/Projects/ITH/primer3/ITH_p3_output.failed

dir="/home/dyap/Projects/ITH/primer3"
# Source and Output directories where files stored
sourcedir=$dir
outdir=$dir

# Part of the redesign pipeline
rename=$name"_p3_redesign.txt"

echo "Please confirm the name of the sample to be reiterated: " $name
read ans

export name=$name

grep "PRIMER_SEQUENCE_ID=" /home/dyap/Projects/ITH/primer3/ITH_p3_output /home/dyap/Projects/ITH/primer3/ITH_p3_passed
cat /home/dyap/Projects/ITH/primer3/$name_p3_output.failed | sed 's/SEQUENCE_TARGET=170,60/SEQUENCE_TARGET=190,20/'| sed 's/PRIMER_PRODUCT_SIZE_RANGE=160-200\ 140-220/PRIMER_PRODUCT_SIZE_RANGE=240-260\ 270-300\ 300-420/'| grep -A19 PRIMER_SEQUENCE_ID= | sed 's/--/=/' > $rename
echo "=" >> $rename



exit;
