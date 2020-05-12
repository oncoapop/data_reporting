#!/bin/sh

# This script does NOT appear to work as you cannot specify command delimited replicates on the commandline

# Name of Project
Project="CX5461"
# Project Directory

clear

dir="/home/dyap/dyap_temp/ChIPseqAnalysis3/"$Project
for line in `cat /home/dyap/Scripts/v1.3_ChIP_pipeline/rep_table`

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo $sample "started processing." >> $dir/macsjobs.log

	echo "Line = " $line

	unset treatname
	unset contname
	unset expt
	unset treat

	expt=`echo $line | awk -F";" '{print $1}'`
	treat=`echo $line | awk -F";" '{print $2}'| sed 's/\,/\\\|/g'`
	cont=`echo $line | awk -F";" '{print $3}' | sed 's/\,/\\\|/g'`
	name=`echo $line | awk -F";" '{print $4"_"$1}'`

	echo $cont

	path=$dir"/EXPT-"$expt
	cd $path

	# get only the sorted.bam files NOT the bam or sorted.bam.bai files or S10 when I want only S1

	treatname=`ls | grep $treat | grep -v bai | grep -v "S1[0-9]" | grep sorted`
	contname=`ls | grep $cont | grep -v bai | grep -v "S1[0-9]" | grep sorted`
	
	t1=`echo $treatname | awk -F" " '{print $1}'`
	t2=`echo $treatname | awk -F" " '{print $2}'`
	t3=`echo $treatname | awk -F" " '{print $3}'`

	c1=`echo $contname |  awk -F" " '{print $1}'`
	c2=`echo $contname |  awk -F" " '{print $2}'`
	c3=`echo $contname |  awk -F" " '{print $3}'`

	echo "Name      = "$name
	echo "==================="

# macs2 diffbdg calling command for ChIP-seq difference of bedgraphs

macs2 bdgdiff --t1 $t1 --c1 $c1 \
              --t2 drug_treat_pileup.bdg \
              --c2 drug_control_lambda.bdg \
              --d1 (depth of seq in cond1) \
              --d2 (depth of seq in cond2) \
              --max-gap 1000 (default:100, larger for broad peaks) \
              --OUTDIR \
              --o-prefix PREFIX (3 files output:PREFIX_cond1.bed, PREFIX_cond2.bed, PREFIX_common.bed) \

############

done
exit

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
