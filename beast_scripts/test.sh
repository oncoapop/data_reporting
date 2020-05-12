#!/bin/sh

dir="/home/dyap/Projects/FFPE"
file=$dir"/gene.txt"
out=$dir"/output.txt"

rm -f $out

for i in `cat $file | awk -F"\t" '{print $1}' | sort -u | sed 's/-//g'`
	do
	echo $i
	cat $file | awk -F"\t" '{print $1}' | grep -o "$i" | wc -l

	echo $i >> $out 
#	cat $file | awk -F"\t" '{print $1}' | grep "$i" >> $out
	cat "$file" |  awk -F"\t" '{ if ($1 == "$i" )  print $o }' >> $out

	echo "=================="
	echo "==================" >> $out
	done

exit;



