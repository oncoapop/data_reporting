#!/bin/sh

# script to count the number of reads per sample
# of fastq.gz files in this directory

path="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/140718_M00897_0139_000000000-A92U5/Data/Intensities/BaseCalls"

###################################

outpath="/home/dyap/Projects/Takeda_T3/CG"
runid=`echo $path | awk -F- '{print $2}' | awk -F"/" '{print $1}'`
outfile=$outpath"/"$runid"-counts"

rm -f $outfile

cd $path

#for i in `ls -I "Undetermined*" | grep ".fastq.gz"` 
for i in `ls *.fastq.gz` 
	do
	echo $i
#	total=`zgrep -c '^@' $i`
	zcat $i | echo $((`wc -l` /4)) >> $outfile
	done

cat $outfile | paste -sd+ - | bc >> $outfile
