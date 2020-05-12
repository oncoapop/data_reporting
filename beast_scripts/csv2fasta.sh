#!/bin/sh

# This Script was written by Damian Yap (Mar 2013)

# Use to convert a single line multi col file into FASTA format
# using only 1 col

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/"

name="VBA0038"

# Source and Output directories where files are stored
csvfile=$name"_positions.csv"
unixfile=$name"-unix.txt"
fafile=$name"_amplicons"

# In the CSV file col 2 Chr, col 3=pos, col 4=seq

cd $dir
rm -f $fafile

tr '\r' '\n' <$csvfile >$unixfile

        for i in `cat $unixfile  | awk -F, '{print $2]'
`
        do
        chr=`grep -m 1 $i $unixfile | awk -F\, '{print $2}'
        pos=`grep -m 1 $i $unixfile | awk -F\, '{print $3}'

        echo ">"$chr":"$pos >> $fafile

# Change $2 for forward direction and $3 for rev complement

        grep -m 1 $i $unixfile | awk -F\, '{print $4}' >> $fafile

        done;


exit;


