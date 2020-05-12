#!/bin/sh

# This Script was written by Damian Yap (Mar 2013)

# Download and install latest blast package from NCBI
# Need to make database with command
# makeblastdb -in database.fsa -parse_seqids -dbtype nucl

# take the output of autoblast in .fa format

# Working Directory
dir="/home/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq"

# Source and Output directories where Barcoded files stored
sourcedir=""
blastout=$dir/Blastfasta

echo Output to this directory $blastout

# get the blast matches

cd $dir/$sourcedir

#	for j in `ls *-*`
#	do
j=Nucleus-001

	out=$blastout/$j


# This gets all the Blast query headers
# incl those with no hits in database
	grep ^Query $j > $out.tmp

	echo $out.tmp

# This gets only those with alignment lines and the line before
# which from the temp file can only be another alignment or
# the query header (this drops any queries with no hits!)
	grep "Query " -B1 $out.tmp > $out.aln

	echo $out.aln

	grep "Query=" $out.aln > $out.aln.tmp

	echo $out.aln.tmp


# To generate template hits

	grep ">" $j  | sed 's/^.*_//' | sort > $out.sort

#	done;


echo Done!

echo Output in the following directories:
echo $blastout

rm -f *.tmp

exit;
