#!/bin/sh

# This scripts takes the indel and SNV orders and combines them for combined order
# adds adaptors and labels them using informative notations from primerlist.txt

# Forward adaptor for Fluidigm
fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
ra="TACGGTAGCAGAGACTTGGTCT"

# Working directory & files
dir="/home/dyap/Projects/"$Project"/primer3"

#source1="TNBC-SNV_primer_order.txt"
source1=$name"_primer_order.txt"

tempdir="/home/dyap/dyap_temp/ITH"
tmp=$tempdir"/"$Project"-SNV.tmp"
tmp3=$tempdir"/"$Project"-SNV-id.tmp"

query1=$name"_primerlist.txt"


dest=$name"-final_order.txt"
final=$name"-final_order.csv"

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

echo "Done."

sams=`grep "^PAT" $dest | awk -F"," '{print $1}' | sort -u -t1`
count=`grep "^PAT" $dest | awk -F"," '{print $1}' | sort -u | wc -l`
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

 
