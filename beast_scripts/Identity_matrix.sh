#/bin/sh

# script to grep identities from blast aligned files

# written by Damian Yap (Apr 2013)
# written on Ubantu 12, now on Redhat Cluster

# Source directory
indir='/home/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq' 
infile=Nucleus-014

# Destination directory
outdir='/home/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq/Blastfasta' 

cd $outdir
rm *-stats
rm *-validreads
rm *-reads.csv
rm *-IDs

cd $indir

	for i in `ls Nucleus*`

	do

	echo $i
	infile=$i

	grep "Identities" $infile | awk -F" " '{print $3}' | sed 's/\/.*//' > $outdir/$infile-stats

	echo "Total number of reads in "$infile" :"
	grep "Identities" $infile | awk -F" " '{print $3}' | sed 's/\/.*//' |  awk '{ SUM += $1} END { print SUM }'

	# grep only reads from 140 - 159 (max 151 for MiSeq) ie from Blast aligned reads
	cat $outdir/$infile-stats | grep "\(\([1][4-5][0-9]\)\|\([1][0-5][0-9]\+\)\)" > $outdir/$infile-validreads

	echo "Valid number of reads: "
	grep -c "\(\([1][4-5][0-9]\)\|\([1][0-5][0-9]\+\)\)" $outdir/$infile-stats

	echo "Total number of reads: "
	grep -c "\(\([0-9]\)\)" $outdir/$infile-stats

	# Get all the Amplicon IDs from reads 140 - 159 
	grep -a4 "Identities = \(\([1][4-5][0-9]\)\|\([1][0-5][0-9]\+\)\)" $infile | grep ">" | sed 's/^.*_//' > $outdir/$infile-IDs

	# Get all the sequences that are from valid aligned reads (starts from 1-49) into one line per read

	grep -a11 "Identities = \(\([1][4-5][0-9]\)\|\([1][0-5][0-9]\+\)\)" $infile | grep "^Query" | sed 's/Query\ \ \(\([0-9]\)\|\([1-4][0-9]\)\)\ \ \ /x x @/' | awk -F" " '{print $3}' | tr -d "\n" | tr "@" "\n" > $outdir/$infile-read.txt

	# Distribution of reads across amplicons
	sort -n $outdir/$infile-IDs | uniq -c > $outdir/$infile-distrib.csv

	echo ++++++++++++++++++++++++++++++++++++++++++++++++++

	done


echo All done.

exit;

