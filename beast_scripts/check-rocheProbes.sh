#!/bin/sh

# This Script was written by Damian Yap (Sep 2015)

dir="/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15/Roche_150807_custom/Iteration_1/starting_fastas_and_selected_probes"
out=$dir"/matches"
notarget=$dir"/no-match-probe"
noprobe=$dir"/no-match-target"

rm -f $out
rm -f $notarget
rm -f $noprobe

# Roadmap
# 1. get the target sequences from $targets
target="/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15/Roche_150807_custom/targets/selected_regions.fa"

# 2. From each of the probes from this file $name, map that unto the $targets
name=$dir"/regions_info_selected_probes.txt"

for i in `cat $name  | awk -F"\t" '{print $1}'`
	do
	echo "======================================================"
	echo $i
	echo "------------------------------------------------------"
	grep --color=always -B1 "$i" $target
	echo "======================================================"

	{
	echo "======================================================"
	echo $i
	echo "------------------------------------------------------"
	grep --color=always -B1 "$i" $target
	echo "======================================================"
	} >> $out

	test=`grep -B1 "$i" $target`
	if [[ -z "$test" ]]
		then
		{
		echo $i
		grep $i $name | awk -F"\t" '{print $2}'		
		echo $test
		echo "======================================================"
		} >> $notarget
	fi

	done

# 3. outputting the stats
# How many regions do not have any probes?
# a. get regions from $target
# b. grep that from $name
# c. output only if there are NO matches


for j in `grep ">" $target | sort -u`
	do
	echo $j
	match=`grep $j $out`
	echo $match
	if [[ -z "$match" ]]
		then
		{
		echo $j
		echo $match
		echo "======================================================"
		} >> $noprobe
	fi
	
	done
exit;


