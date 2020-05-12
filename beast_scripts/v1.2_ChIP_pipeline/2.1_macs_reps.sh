#!/bin/sh

# This script does NOT appear to work as you cannot specify command delimited replicates on the commandline

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis3/"$Project
for line in `cat /home/dyap/Scripts/v1.2_ChIP_pipeline/rep_table`

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo $sample "started processing." >> $dir/macsjobs.log

	echo "Line = " $line

	expt=`echo $line | awk -F";" '{print $1}'`
	treat=`echo $line | awk -F";" '{print $2}'| sed 's/\,/\\\|/g'`
	cont=`echo $line | awk -F";" '{print $3}' | sed 's/\,/\\\|/g'`
	name=`echo $line | awk -F";" '{print $4"_"$1}'`


	path=$dir"/EXPT-"$expt
	cd $path

	# get only the sorted.bam files NOT the bam or sorted.bam.bai files

	treatname=`ls | grep $treat | grep -v bai | grep sorted | xargs | sed -e 's/ /,/g'`
	contname=`ls | grep $cont | grep -v bai | grep sorted | xargs | sed -e 's/ /,/g'`
	echo "Treatname = "$treatname
	echo "Control   = "$contname
	echo "Name      = "$name
	echo "==================="

# macs2 peak calling command for ChIP-seq

# Sam's modelling
#I recommend using the --SPMR flag in MACS2 which normalises the coverage plot values by millions of reads.  This allows you to compare visually between tracks if the samples are very different in terms of total reads.

#/home/dyap/bin/macs2 callpeak --broad --nomodel -t $treatname -c $contname -f BAM -n $name -B --SPMR -g hs --extsize=250 
echo "macs2 callpeak --broad --nomodel -t $treatname -c $contname -f BAM -n $name -B --SPMR -g hs --extsize=250" > $dir/$name.log 


#############
echo "Renaming bedgraphs..."

# rename .bdg to .bedgraph and add respective track names
for f in `ls *.bdg`
	do 
        out=`echo $f | sed 's/bdg/bedgraph/g'`
	name=`echo $out | sed 's/_treat_pileup.bedgraph//g'`
	echo $name
	header="track type=bedGraph name="$name" description="$name
	echo "$header"  > /home/dyap/dyap_temp/$f.tmp
	cat "$f" >> /home/dyap/dyap_temp/$f.tmp
    	mv /home/dyap/dyap_temp/$f.tmp "$out"
	done

	echo "==========="

echo "done"	
done

exit

##############
