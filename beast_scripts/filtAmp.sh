#!/bin/sh

# Script to filter the ordered primers 
# or list of positions from the AmpliconManifests

# Important for filtering from primer order file

# Forward adaptor for Fluidigm
fa="ACACTGACGACATGGTTCTACA"
# Reverse adaptor for Fluidigm (5'->3')
ra="TACGGTAGCAGAGACTTGGTCT"

# Illlumina Adaptors (5'->3')
#fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
#ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

temp="/home/dyap/dyap_temp"

Project="Tumour_Xenograft_Rev-re"
Sample="SA530"

#basedir="/home/dyap/Projects/"$Project"/"$Sample
basedir="/home/dyap/Projects/PrimerDesign/"$Project"/primer3"

inputfile=$basedir"/"$Sample"_designed_positions.txt"

cat $basedir"/"$Sample"_pos.txt" | awk -F"\t" '{print "chr"$1"_"$2}' > $inputfile

#manifestfile=$basedir"/"$Sample"-"$set".AmpliconManifest"
#orderfile=$basedir"/"$Sample"-"$set"-primer-order.csv"

oldmanifestfile=$basedir"/"$Sample".AmpliconManifest"
#orderfile=$basedir"/"$Sample"-"$set"-primer-order.csv"

manifestfile=$basedir"/"$Sample"-filtered.AmpliconManifest"
#orderfile=$basedir"/"$Sample"-"$set"-primer-order.csv"

manfile=$temp"/"$Sample"-manifest.tmp"
#ordfile=$temp"/"$Sample"-order.tmp"

# Amplicon Manifest Header
 		echo "[Header]" > $manfile
                echo $Sample"  Manifest Version,2" >> $manfile
                echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
                echo "  " >> $manfile
                echo "[Regions]" >> $manfile
                echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manfile

# Primer Order file
#		echo "WellPosition,Name,Sequence,Notes" > $ordfile
#		count=0

# Get positions from selected file and generate manifest and orderfile


for i in `cat $inputfile | awk -F, '{print $1}'`

	do

	id=`grep $i $inputfile | awk -F, '{print $1}'`
	echo $id

	 # Generate the Amplicon Manifest file 
	 grep $i $oldmanifestfile >> $manfile

	done

cat $manfile | tr "," "\t" > $manifestfile


exit;

