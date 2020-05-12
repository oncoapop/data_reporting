#!/bin/sh

# Script that will check the fastq files for PCR primers 
# to see if these will be expressed.
# These must be RNA-seq files.

# These are the untreated controls
path="/share/lustre/archive/SA464/illumina_wtss/A14559/sequence"
cd $path

primerpath="/home/dyap/Projects/PrimerDesign/Splice/primer3"
primerfile=$primerpath"/hct116_htert_primer_order.txt"

name=`echo $path | awk -F"/" '{print $5}'`
outfile="/home/dyap/Projects/Takeda_T3/primerQC/"$name".txt"
rm -f $outfile

for file in `ls *2.fastq` 
	do
	echo "File query= "$file
	
	for j in `cat $primerfile | awk -F"," '{print $1}' | tail -n +2`
		do
		echo $j
		fwd=`grep $j $primerfile | awk -F"," '{print $6}'`
		fc=`grep $fwd $file | wc -l`
		fwdrc=`echo $fwd | awk  'BEGIN {
                              j = n = split("A C G T", t)
                              for (i = 0; ++i <= n;)
                              map[t[i]] = t[j--]
                                      }
                      {
                              if (/LEFT/) print
                      else {
                              for (i = length; i; i--)
                              printf "%s", map[substr($0, i, 1)]
                              print x
                              }
                      }'`

		frcc=`grep $fwdrc $file | wc -l`

		rev=`grep $j $primerfile | awk -F"," '{print $8}'`
		rc=`grep $rev $file | wc -l`
		revrc=`echo $rev | awk  'BEGIN {
                              j = n = split("A C G T", t)
                              for (i = 0; ++i <= n;)
                              map[t[i]] = t[j--]
                                      }
                      {
                              if (/LEFT/) print
                      else {
                              for (i = length; i; i--)
                              printf "%s", map[substr($0, i, 1)]
                              print x
                              }
                      }'`
		rcc=`grep $revrc $file | wc -l`

		echo $file","$j,"fwd="$fc"+"$frcc",rev="$rc"+"$rcc >> $outfile
		echo $file","$j,"fwd="$fc"+"$frcc",rev="$rc"+"$rcc
		done
	done
