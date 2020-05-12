#!/bin/sh

# Script that will check the fastq files for PCR primers 
# to see if the respective genes in the cojoined gene pair is expressed
# NOTE: CG link is not assumed
# This is a QC step to shortlist primers rather than an CG analysis step

# reads lib location from
library="/home/dyap/Projects/Takeda_T3/primerQC/lib"
# format SA ID<TAB>Seq Library

# This is the whole batch of primers designed
primerpath="/home/dyap/Projects/PrimerDesign/Splice/primer3"
primerfile=$primerpath"/hct116_htert_primer_order.txt"

#Formatting the path on beast
for id in `cat $library | awk -F"\t" '{print $1}'`
	do
	lib=`grep $id $library | awk -F"\t" '{print $2}'`	
        path="/share/lustre/archive/"$id"/illumina_wtss/"$lib"/sequence"
	outfile="/home/dyap/Projects/Takeda_T3/primerQC/batch/"$id

# If run the second time, uncomment this
#	rm -f $outfile	

	# Formatting the path on beast and the processed fastq file
	cd $path
	for file in `ls *2.fastq` 
		do
		echo "File query= "$file

	# Getting the left and right primers (AND reverse complement)
	# from the primer design file (which has a header)
        for j in `cat $primerfile | awk -F"," '{print $1}' | tail -n +2`
                do
                echo $j
                fwd=`grep $j $primerfile | awk -F"," '{print $6}'`
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

		# count the occurance of the forward (& rc) primers in fastq
                fc=`grep -c $fwd $file`
                frcc=`grep -c $fwdrc $file`


                rev=`grep $j $primerfile | awk -F"," '{print $8}'`
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

		# count the occurance of the reverse (& rc) primers in fastq
                rcc=`grep -c $revrc $file`
                rc=`grep -c $rev $file`

                echo $j,"fwd="$fc"+"$frcc",rev="$rc"+"$rcc >> $outfile
                echo $j,"fwd="$fc"+"$frcc",rev="$rc"+"$rcc

                done

		done
        done
