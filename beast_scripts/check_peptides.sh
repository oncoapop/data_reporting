#!/bin/sh

cd /home/dyap/Projects/Takeda_T3/TiO2

for i in `ls Predicted_CLK*`
	do
	echo $i
	for j in `cat $i | awk -F"\t" '{print $1}' | sort -u`
		do
		if [[ $j =~ "Source=" ]]
			then continue
		fi
		count=`grep -w "$j" $i | wc -l`
		echo $count
		for k in  $(eval echo "{1..$count}")
			do
				gene=`grep -m$k -w "$j" $i | tail -n1 | awk -F"\t" '{print $1}'`
				seq=`grep -m$k -w "$j" $i | tail -n1 | awk -F"\t" '{print $2}'`
				grep $seq "PeptideSeq_TiO2Expt.csv"  | awk -F, '{print $2}' | awk -F" " '{print $3}'
			done
		done

	done
