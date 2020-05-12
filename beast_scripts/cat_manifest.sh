#!/bin/sh

# Script to cat two manifest files (less headers)
# into one new manifest file
# Generate the AmpliconManifest for MiSeq run

temp="/home/dyap/dyap_temp"

Project="manual"
Sample="SA501"
set="combined"

echo $Project,$Sample,$set

basedir="/home/dyap/Projects/PrimerDesign/"$Project
file1=$basedir"/SA501-set1.AmpliconManifest"
file2=$basedir"/SA501-set2.AmpliconManifest"

manifestfile=$basedir"/"$Sample"-"$set".AmpliconManifest.tmp"
final=$basedir"/"$Sample"-"$set".AmpliconManifest"

manfile=$temp"/"$Sample"-manifest.tmp"

rm -fr $manfile

# Amplicon Manifest Header
 		echo "[Header]" > $manifestfile
                echo $Sample"-"$set"  Manifest Version,1" >> $manifestfile
                echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manifestfile
                echo "  " >> $manifestfile
                echo "[Regions]" >> $manifestfile
                echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manifestfile

# Get positions from selected file and generate manifest and orderfile

inputfile=$file1
#cat $inputfile

	for i in `cat $inputfile | awk -F"\t" '{print $1}'`

		do
		id=`grep -m1 $i $inputfile | awk -F"\t" '{print "SA501_"$1}' | tr "C" "c"`
		echo $id
		chr=`grep -m1 $i $inputfile | awk -F"\t" '{print $2}'`
		sta=`grep -m1 $i $inputfile | awk -F"\t" '{print $3}'`
		end=`grep -m1 $i $inputfile | awk -F"\t" '{print $4}'`
		llen=`grep -m1 $i $inputfile | awk -F"\t" '{print $5}'`
		rlen=`grep -m1 $i $inputfile | awk -F"\t" '{print $6}'`
		comm=`grep -m1 $i $inputfile | awk -F"\t" '{print $7}'`

		 # Generate the Amplicon Manifest file 
		 echo $id","$chr","$sta","$end","$llen","$rlen","$comm >> $manfile

		done

echo "+++++++++++"

inputfile=$file2
	for i in `cat $inputfile | awk -F"\t" '{print $1}'`

		do
		id=`grep -m1 $i $inputfile | awk -F"\t" '{print $1}'`
		echo $id
		chr=`grep -m1 $i $inputfile | awk -F"\t" '{print $2}'`
		sta=`grep -m1 $i $inputfile | awk -F"\t" '{print $3}'`
		end=`grep -m1 $i $inputfile | awk -F"\t" '{print $4}'`
		llen=`grep -m1 $i $inputfile | awk -F"\t" '{print $5}'`
		rlen=`grep -m1 $i $inputfile | awk -F"\t" '{print $6}'`
		comm=`grep -m1 $i $inputfile | awk -F"\t" '{print $7}'`

		 # Generate the Amplicon Manifest file 
		 echo $id","$chr","$sta","$end","$llen","$rlen","$comm >> $manfile

		done



cat $manifestfile | tr "," "\t" > $final
cat $manfile | tr "," "\t" | sort -k 1 >> $final


exit;

