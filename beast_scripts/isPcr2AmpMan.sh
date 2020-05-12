#!/bin/bash
# Script to get the top hits of isPCR
# Generated from the commandline isPCR on beast
# Generate an order file for those primers - done
# Generate Manifest - done
# Generate HTML view file - done
# scp that to godel.cluster.bccrc.ca - done

# IF run as pipeline (comment this out)
Project="MitoVar"
# Required
name="mito"
# basedir
#dir="/home/dyap/Projects/PrimerDesign/"$Project
#dir="/home/dyap/Projects/"$Project"/"$name
dir="/home/dyap/Projects/Single_Cell/"$Project

reads=200

#==================================================
# Project Directory
# dir is imported from pipeline script in the environment
htmloutpath="/home/dyap/public_html"

# positions
#posdir=$dir"/positions"
posdir=$dir

# primer3 output
#p3dir=$dir"/primer3"
p3dir=$dir

# Annotation output
#annodir=$dir"/Annotate"
annodir=$dir

# Tmp file
tmp="/home/dyap/dyap_temp"

echo "Running Primer_Summary version v4.01..."

# from isPCR command line output
# fasta formatted
# infilesuffix="_isPCR-output.fa"

# from isPCR automated web iterface
# already parsed
# preinfilesuffix="_isPCR_chk.csv"
infilesuffix="_isPCR-output.txt"

outfilesuffix="_primer_order.txt"

# Name of the input file
inputfile=$p3dir"/"$name$infilesuffix

# Name of the output file
orderfile=$p3dir"/"$name$outfilesuffix
manifestfile=$p3dir"/"$name".AmpliconManifest"
qcfile=$tmp"/"$name"_QC.txt"
manfile=$tmp"/"$name"_manifest"
ordfile=$tmp"/"$name"_order"
suppfile=$tmp"/"$name"_supp"
manver=2

# Checking to see if required files exist in correct paths

if [ -f $manifestfile ]
	then
	echo "Overwriting "$manifestfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
fi
		echo "[Header]" > $manfile
		echo $name"  Manifest Version,"$manver >> $manfile
		echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
#		echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\rCRS-fasta" >> $manfile
#		echo "ReferenceGenome,Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
#		echo "ReferenceGenome,Custom\gp140" >> $manfile
		echo "	" >> $manfile
		echo "[Regions]" >> $manfile
		echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manfile

if [ -f $qcfile ]
	then
	echo "Overwriting "$qcfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
	rm -f $qcfile
fi
		echo "Amplicon ID,Chromosome,Amplicon Start,Amplicon End,Left Primer,Right Primer,Amplicon Length" > $suppfile


for i in `grep ">" $inputfile | awk -F" " '{print $2}' | tr -d ">" | tr -d '"' | tail -n +2`
	do
		# QC to test if grep is successful
		test=`grep $i $inputfile | wc -l`

		if [ $test -eq 0 ]
			then
			echo "Pattern not found"
			failed=`echo "$failed + 1" | bc`
			continue

		fi
		if [ $test -gt 5 ]
			then
			echo "Possible non-specific amplicon"
			grep $i $inputfile
			note="5 primers pairs - "$test" amplicons"
			else
			note=" "
		fi

	# This module is specific for command-line isPcr output

	ampchr=`grep -m1 $i $inputfile | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`

        plusstrand=`grep -m1 $i $inputfile | grep [0-9]+[0-9] | wc -c`
        minusstrand=`grep -m1 $i $inputfile | grep [0-9]-[0-9] | wc -c`


     		if [[ $plusstrand == 0 ]] && [[ $minusstrand == 0 ]];
                		then
				break
       		fi

        if [[ $plusstrand > 0 ]] && [[ $minusstrand == 0 ]];
                then
                ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
                ampend=`grep -m1 $i $inputfile |  awk -F" " '{print $1}' | awk -F":" '{print $2}' | awk -F"+" '{print $2}'`
        fi

        if [[ $minusstrand > 0 ]]  && [[ $plusstrand == 0 ]];
                then
                ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
                ampend=`grep -m1 $i $inputfile |  awk -F" " '{print $1}' | awk -F":" '{print $2}' | awk -F"-" '{print $2}'`
        fi

        myid=`grep -m1 $i $inputfile | awk -F" " '{print $2}'`
        amplen=`grep -m1 $i $inputfile | awk -F" " '{print $3}' | sed 's/bp//'`
        leftpri=`grep -m1 $i $inputfile | awk -F" " '{print $4}'`
        rightpri=`grep -m1 $i $inputfile | awk -F" " '{print $5}'`

echo $myid
echo "==================================="

# QC for MiSeq run (specific to input)
# Exported from the $Miseq
	reads=150
# assume width=20


	llim=`echo "$ampstart + $reads - 10 " | bc`
	rlim=`echo "$ampend - $reads + 10 " | bc`

        if [[ $llim -gt $rlim ]];
		then
		overlap="overlap"
		else
		overlap="NO OVERLAP"
	fi
        
	echo -e $myid"\t"$ampchr"\t"$ampstart"\t"$ampend"\t"$amplen"\t QC: "$reads"bp: "$Miseq >> $qcfile

		# Generate manifest file of all that pass

		llen=`echo $leftpri | wc -c`
		leftlen=`echo "$llen - 1" | bc`
		rlen=`echo $rightpri | wc -c`
		rightlen=`echo "$rlen - 1" | bc`
		echo $myid,$ampchr,$ampstart,$ampend,$leftlen,$rightlen,$amplen" bp","leftrange="$ampstart"-"$llim,"rightrange="$rlim"-"$ampend","$overlap >> $manfile
#		echo $myid,$ampchr,$ampstart,$ampend,$leftlen,$rightlen,$amplen" bp" >> $manfile

		# Supplemental Figure information
		echo $myid","$ampchr","$ampstart","$ampend","$leftpri","$rightpri","$amplen"bp" >> $suppfile

	done
echo "DONE."
echo "#####"


cat $manfile | tr "," "\t" > $manifestfile"_QC"
echo "Manifest File generated: " $manifestfile"_QC"

rm -f $p3dir"/*.tmp"

exit;


