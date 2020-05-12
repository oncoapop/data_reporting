#!/bin/sh

# This script gets the mutations from mutation seq run in single sample mode
# This allows for germ line as well as somatic mutations to be listed

tumpath="/share/lustre/archive/BOB/ctDNA/OUTPUT_tumor/RUN/"
nompath="/share/lustre/archive/BOB/ctDNA/OUTPUT_normal/RUN/"
outpath="/home/dyap/Projects/ctDNA/AT094/validation/mutations"

# This first gets all the mutations into a directory
# In the first instance all are unfiltered

cd $tumpath
rm -f $outpath/*.txt

for i in `ls | grep F`
	do
	sam=`echo $i | awk -F"_" '{print $1"_"$2}'`

	cd $tumpath$i"/outputs/results"
	file=`ls | grep TASK_20`

	if [ -s "$file" ]
		then

		outfile=$outpath"/"$sam".txt"
		echo $outfile

		cat $file | awk -F"\t" '{print $2}' | tail -n +2 | sort >> $outfile
		
	fi


	done 


cd $nompath
echo $nompath

for i in `ls | grep G`
	do
	sam=`echo $i | awk -F"_" '{print $1"_"$2}'`

	cd $nompath$i"/outputs/results"
	file=`ls | grep TASK_20`

	if [ -s "$file" ]
		then

		outfile=$outpath"/"$sam".txt"
		echo $outfile

		cat $file | awk -F"\t" '{print $2}' | tail -n +2 | sort >> $outfile
		
	fi


	done 

cd $outpath

# now we have to make sure mutations are germline - ie in all samples
for j in `ls | grep .txt`
	do
	sam=`echo $j | awk -F"_" '{print $1}'`
	
	echo $sam >> sample_names.txt
	
	done

cd $outpath

for k in `cat sample_names.txt | sort -u`
	do
	echo $k
	match=`ls | grep $k`

	for l in `echo $match`
		do
		for m in `cat $l`
			do
			echo $m >> $k"_mut.txt"

			done
		
		done		


	done


for m in `ls *mut.txt`
	do
	for n in `cat $m`
		do

		pos=`grep $n $m`
		count=`echo $pos | wc -w`
		sample=`echo $m | awk -F"_" '{print $1}'`

			if [ $count -gt 1 ]
				then
				
				echo $sample" "$pos >> summary.txt

			fi 

		done

	done


# 
