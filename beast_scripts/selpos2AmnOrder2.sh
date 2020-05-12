#!/bin/sh

# Script to take the order form which contains selected positions
# from isPCR output
# Generate the AmpliconManifest for MiSeq run
temp="/home/dyap/dyap_temp"

Project="manual"
Sample="DAH55_56_Shared"

basedir="/home/dyap/Projects/PrimerDesign/"$Project
inputfile=$basedir"/"$Sample"_primers.txt"
isPCR=$basedir"/"$Sample"_isPCR-output.txt"

manifestfile=$basedir"/"$Sample"-2.AmpliconManifest"

manfile=$temp"/"$Sample"-manifest.tmp"

# Amplicon Manifest Header
 		echo "[Header]" > $manfile
                echo $Sample"  Manifest Version,2" >> $manfile
#                echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
                echo "ReferenceGenome,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
                echo "  " >> $manfile
                echo "[Regions]" >> $manfile
                echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manfile


# Get positions from selected file and generate manifest and orderfile

for i in `cat $inputfile | awk -F" " '{print $1}' | sort -u`

	do
		echo $i
		id=`grep -m1 $i $isPCR | awk -F" " '{print $2}'`
		echo $id

		snv=`echo $id | awk -F"_" '{print $NF}'`
		chr=`echo $id | awk -F"_" '{print $NF-1}'`

                ampchr=`grep -m1 $i $isPCR | grep ">chr" | awk -F":" '{print $1}' | sed 's/>//g'`

                plussta=`grep -m1 $i $isPCR | grep ">chr" | awk -F":" '{print $2}'| awk -F" " '{print $1}' | awk -F"+" '{print $1}'`
                plusend=`grep -m1 $i $isPCR | grep ">chr" | awk -F":" '{print $2}'| awk -F" " '{print $1}' | awk -F"+" '{print $2}'`

                minussta=`grep -m1 $i $isPCR | grep ">chr" | awk -F":" '{print $2}'| awk -F" " '{print $1}' | awk -F"-" '{print $1}'`
                minusend=`grep -m1 $i $isPCR | grep ">chr" | awk -F":" '{print $2}'| awk -F" " '{print $1}' | awk -F"-" '{print $2}'`

                if [[ $plussta > 0 ]] && [[ $plusend > 0 ]];
                        then
                        ampsta=$plussta
                        ampend=$plusend
                fi

                if [[ $minussta > 0 ]] && [[ $minusend > 0 ]];
                        then
                        ampsta=$minussta
                        ampend=$minusend
                fi

                fwdpri=`grep -m1 $i $isPCR  | grep ">chr" | awk -F" " '{print $4}'`
                revpri=`grep -m1 $i $isPCR  | grep ">chr" | awk -F" " '{print $5}'`


                llen=`echo $fwdpri | wc -c`
                leftlen=`echo "$llen - 1" | bc`
                rlen=`echo $revpri | wc -c`
                rightlen=`echo "$rlen - 1" | bc`

	 # Generate the Amplicon Manifest file 
	 echo $id","$ampchr","$ampsta","$ampend","$leftlen","$rightlen >> $manfile

	done

cat $manfile | tr "," "\t" > $manifestfile


exit;

