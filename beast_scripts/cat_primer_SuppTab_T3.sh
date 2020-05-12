#!/bin/sh

# Script to summarize the Supplmental table for primers

wd="/home/dyap/Projects/PrimerDesign/Splice/primer3"
cd $wd

pat1="Splice-RT_SupplFig.csv"
pat2="Splice-RT2_SupplFig.csv"
ordered="ordered_RT.txt"
output="Suppl_Tab_CG_Expr.tmp"
finaloutput="Suppl_Tab_CG_Expr.csv"

rm -f $output

cat $pat1 | tail -n +2 > all_primers.tmp
cat $pat2 | tail -n +2 >> all_primers.tmp

for i in `cat $ordered`
	do
	transcript=`echo $i | awk -F"_" '{print $1}'`
	transcript2=`echo $i | awk -F"_" '{print $2}'`
	transcript3=`echo $i | awk -F"_" '{print $3}'`
	echo $transcript,$transcript2,$transcript3

	grep $transcript all_primers.tmp | grep $transcript2 | grep $transcript3 >> $output
	done

for j in `cat $output | awk -F"_" '{print $1}'`
	do
	gene=`grep -m1 $j $output | awk -F"_" '{print $1}'`
	transcript=`grep -m1 $j $output | awk -F, '{print $2}'`
	ampstart=`grep -m1 $j $output | awk -F, '{print $3}'`
	ampend=`grep -m1 $j $output | awk -F, '{print $4}'`
	leftpri=`grep -m1 $j $output | awk -F, '{print $5}'`
	rightpri=`grep -m1 $j $output | awk -F, '{print $6}'`
	amplen=`grep -m1 $j $output | awk -F, '{print $7}'`

	echo $gene","$leftpri","$rightpri","$transcript","$ampstart","$ampend","$amplen >> $output.tmp

	done

cat $output.tmp | sort -u  > $finaloutput

rm -f *.tmp


