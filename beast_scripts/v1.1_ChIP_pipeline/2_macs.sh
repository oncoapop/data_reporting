#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis2/"$Project

for line in `cat /home/dyap/Scripts/v1.1_ChIP_pipeline/table`

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo $sample "started processing." >> $dir/macsjobs.log

	echo $line

	expt=`echo $line | awk -F, '{print $1}'`
	treat=`echo $line | awk -F, '{print $2}'`
	cont=`echo $line | awk -F, '{print $3}'`
	name=`echo $line | awk -F, '{print $4"_"$1}'`


	path=$dir"/EXPT-"$expt
	cd $path
	treatname=`ls | grep $treat".sorted.bam\$"`
	contname=`ls | grep $cont".sorted.bam\$"`
	echo $treatname
	echo $contname
	echo $name


# macs2 peak calling command for ChIP-seq

# Damian's command

#/home/dyap/bin/macs2 callpeak --nomodel --extsize 250 -t $path/$treatname -c $path/$contname --broad -g hs -n $name -B 
#echo "macs2 callpeak --keep-dup 1 --shift 73 --nomodel --extsize 250 -t $treatname -c $contname --broad -g hs -n $name -B" >> $dir/macsjobs.log

#-nomodel -shift 37 - extsize 73

# Sam's modelling
#I recommend using the --SPMR flag in MACS2 which normalises the coverage plot values by millions of reads.  This allows you to compare visually between tracks if the samples are very different in terms of total reads.

/home/dyap/bin/macs2 callpeak --broad --nomodel -t $path/$treatname -c $path/$contname -f BAM -n $name -B --SPMR -g hs --extsize=290 
echo "macs2 callpeak --broad --nomodel -t $treatname -c $contname -f BAM -n $name -B --SPMR -g hs --extsize=290" > $dir/$name.log 

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
