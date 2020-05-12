#!/bin/sh

# This scripts takes the indel and SNV orders and combines them for combined order
# adds adaptors and labels them using informative notations from primerlist.txt

# Forward adaptor for Fluidigm
fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
ra="TACGGTAGCAGAGACTTGGTCT"

# Working directory & files
dir="/home/dyap/Projects/TNBC/primer3"

source1="TNBC-SNV_primer_order.txt"
source2="TNBC-indel_primer_order.txt"

tempdir="/home/dyap/dyap_temp"
tmp=$tempdir"/TNBC-SNV.tmp"
tmp3=$tempdir"/TNBC-SNV-id.tmp"
tmp2=$tempdir"/TNBC-indel.tmp"
tmp4=$tempdir"/TNBC-indel-id.tmp"

query1="TNBC-SNV_primerlist.txt"
query2="TNBC-indel_primerlist.txt"


dest="TNBC-combined_order.txt"
final="TNBC-combined_order.csv"

cd $dir

echo "Sample,PosID,left_primer,right_primer,Amplen:ProjID" > $dest

##########################
echo "Processings SNVs..."

cat $source1 | tr -d " " | awk -F"," '{print $1","$6","$8","$5}' > $tmp
cat $source1 | tr -d " " | awk -F"," '{print $1}' > $tmp3

	for i in `cat $tmp3`
		do
		{
		myid=`grep $i $query1 | awk -F"," '{print $2}' | sort -u`
		snvid=`grep $i $tmp | awk -F"," '{print $1}' | sort -u`
		left=`grep $i $tmp | awk -F"," '{print $2}'`
		right=`grep $i $tmp | awk -F"," '{print $3}'`
		leng=`grep $i  $tmp | awk -F"," '{print $4}'`
		lpri=$fa$left
		rpri=$ra$right
		sample=`echo $snvid | awk -F"-" '{print $1}'`
		
		echo $sample","$myid","$lpri","$rpri","$leng":"$snvid >> $dest
		}
	done

############################
echo "Processings Indels..."

cat $source2 | tr -d " " | awk -F"," '{print $10","$6","$8","$5}' > $tmp2
cat $source2 | tr -d " " | awk -F"," '{print $10}' > $tmp4

	for j in `cat $tmp4`
		do
		{
		myid=`grep $j $query2 | awk -F"," '{print $2}'| sort -u`
 			if [[ ${myid} -eq "" || ${myid} -eq "PosID"  ]];
	 			then
					continue 
 			fi
		indelid=`grep $j $tmp2 | awk -F"," '{print $1}'`
		left=`grep $j $tmp2 | awk -F"," '{print $2}'`
		right=`grep $j $tmp2 | awk -F"," '{print $3}'`
		leng=`grep $j $tmp2 | awk -F"," '{print $4}'`
		lpri=$fa$left
		rpri=$right$ra
		sample=`echo $indelid | awk -F"-" '{print $1}'`
		
		echo $sample","$myid","$lpri","$rpri","$leng":"$indelid >> $dest
		}
	done

echo "Done."

sams=`grep "^SA" $dest | awk -F"," '{print $1}' | sort -u -t1`
count=`grep "^SA" $dest | awk -F"," '{print $1}' | sort -u | wc -l`
echo "=============================="
echo "Summary for "$count" samples"
echo "=============================="

	for k in `echo $sams`
		do
		{
		scnt=`grep -c $k $dest`
		echo $k = $scnt
		echo Plates:
		echo "$scnt / 96" | bc -l
		echo "=============================="
		}
		done

cat $dest | awk -F"," '{print $2","$3","$4","$5}' > $final

# Count number of bases in the primers with adaptors
lc=`cat $final | awk -F, '{print $2}' | wc -c`
rc=`cat $final | awk -F, '{print $3}' | wc -c`
pp=`cat $final | awk -F, '{print $1}' | wc -l`
bases=`echo "$lc + $rc" | bc`
costperbase=0.08

echo "-----------------------------"

echo "Number of primer pairs = "$pp
echo "No of bases (Primers+adaptors) = "$bases
echo "Cost per base / 500pM scale synthesis = $"$costperbase 
cost=`echo "$bases * $costperbase" | bc`
echo "Cost of ordering primers = $ "$cost

echo "-----------------------------"
   
exit;

 
