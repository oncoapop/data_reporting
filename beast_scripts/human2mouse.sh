#!/bin/sh

# Script to get mouse homologue of the human genes

wd="/home/dyap/Projects/eIF4A3_NMD/comparisons"
cd $wd

match="Mao_2016_Eif4a3_human.csv"
# Convert using tr '\r' '\n' < macfile.txt > unixfile.txt
mu="Mao_2016_Eif4a3.txt"

# Get header and add extra column at the end for human gene
head -n1 $mu  | awk -F"," '{print $0",HumanGene"}' > $match

for i in `cat $mu | tail -n +2 |  awk -F, '{print $2}' | sort -u `
	do
	mou=`echo $i`
#	echo $mou
	hum=`grep -w -m1 -A1 "$i" HOM_MouseHumanSequence.rpt | tail -n +2 | awk -F"\t" '{print $4}'`

	if [[ "$hum" != "" && "$mou" != "" ]]
		then
	grep -w -m1 "$i" $mu |  awk -v var=$hum 'BEGIN{FS=" "}{print $0 "," var}' >> $match
#		echo $hum
		else
	grep -w -m1 "$i" $mu |  awk -v var="N/A" 'BEGIN{FS=" "}{print $0 "," var}' >> $match
#		echo "NOOOooooooooooooooooooo!!!!"
	fi

	done

exit;

	echo $hum","$mou >> $match

	done

   
