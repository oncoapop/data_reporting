#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis/"$Project

for line in `cat table`

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo $sample "Converting bdg->bedgraph." >> $dir/macsjobs.log

	echo $line

	expt=`echo $line | awk -F, '{print $1}'`
	treat=`echo $line | awk -F, '{print $2}'`
	cont=`echo $line | awk -F, '{print $3}'`
	name=`echo $line | awk -F, '{print $4"_"$1}'`


	path=$dir"/EXPT-"$expt
	cd $path

        out=`ls | grep $name | grep "_treat_pileup.bdg" | sed 's/bdg/bedgraph/g'`
	echo $name
	header="track type=bedGraph name="$name" description="$name
	echo "$header"  > /home/dyap/dyap_temp/$out.tmp
	cat "$f" >> /home/dyap/dyap_temp/$out.tmp
    	mv /home/dyap/dyap_temp/$out.tmp "$out"

	rm -f /home/dyap/dyap_temp/$out.tmp

	echo "==========="

echo "done"	
done

exit

##############
