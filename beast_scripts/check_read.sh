#/bin/sh

pwd="/share/scratch/amazloomian_temp/MiSeq/allAlignments"

cd $pwd

for i in `ls *.counts`
	do
	first=`cat $i | awk -F" " '{print $4}' |  awk '{ SUM += $1} END { print SUM }'`
	second=`cat $i | awk -F" " '{print $5}' |  awk '{ SUM += $1} END { print SUM }'`
	third=`cat $i | awk -F" " '{print $6}' |  awk '{ SUM += $1} END { print SUM }'`
	
	sum=`echo "$first+$second+$third" | bc`
	echo -e $i "\t=\t"$first"\t" $second"\t"$third"\t"$sum


	done

exit;

