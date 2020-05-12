#!/bin/bash

# Script to summarise FACS output imported for FlowJo

# Damian Yap 
# 18 Aug 2014

# Example of format
# Aug 18 2014_57.fcs              4176
#   FSC-A, SSC-A subset  96.8    4042
#      <V450-40-A> subset        83      3355
#         <mCherry-A> subset     0       0
#            % RED EdU +ve               0
# Aug 18 2014_58.fcs              4188
#   FSC-A, SSC-A subset  96.4    4038
#      <V450-40-A> subset        96      3878
#        <mCherry-A> subset     0       0
#            % RED EdU +ve               0

wd="/home/dyap/Projects/Mammary"
fname=$wd"/FACS_180814"
outfile=$wd"/Edu+ve_180814.csv"
keyword="Aug"

rm -f $outfile
#unique=`grep $keyword $fname | awk -F. '{print $1}' | awk -F_ '{print $2}' | sort -n | uniq | tail -n +2`
unique=`grep $keyword $fname | awk -F" " '{print $3}'`

echo $unique

echo "+"

for i in $unique
	do
	echo $i
	name=`grep "$i" $fname | awk -F. '{print $1}' | awk -F_ '{print $2}' `
	edu=`grep -A4 "$i" $fname | grep EdU | awk -F"\t" '{print $2}'`

	echo $name","$edu >> $outfile

	done

exit

