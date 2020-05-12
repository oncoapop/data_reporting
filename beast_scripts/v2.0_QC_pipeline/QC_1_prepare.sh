#!/bin/sh

# Script to copy all the vcfs into the QC directory

######################
# ENTER RUN IDs here #
miseq="AMU63"
#####################

QCdir="/home/dyap/Projects/HCT116"

outdir=$QCdir"/"$miseq

if [ ! -d "$outdir" ]; then
	mkdir $outdir
fi

#if [ ! -d $outdir"/miseq" ]; then
#	mkdir $outdir"/miseq"
#fi


path1="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files"
path2="/Data/Intensities/BaseCalls/Alignment"

cd $path1
msname=`ls | grep $miseq`
echo $msname

msc=`ls $path1"/"$msname$path2 | grep ".vcf"$ | grep "_S"`
echo $msc

for i in $msc
	do
	
	cd $outdir"/"
	cp $path1"/"$msname$path2"/"$i $outdir"/"
	echo here
	done

cp $path1"/"$msname$path2"/SampleSheetUsed.csv" $outdir"/"
cp $path1"/"$msname"/"*".AmpliconManifest" $outdir"/"
cp $path1"/"$msname"/SampleSheet.csv" $outdir"/"

exit;

