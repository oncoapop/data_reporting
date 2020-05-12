#!/bin/sh
# Script to get the top hits of isPCR
# Generate an order file for those primers - done
# Generate Manifest - done
# Generate HTML view file - not done
# scp that to godel.cluster.bccrc.ca - not done

# IF run as pipeline (comment this out)
# Project="TOV"
# type="indel"
# type="SNV2"
# Required
# name=$Project"-"$type

# Project Directory
dir="/home/dyap/Projects/PrimerDesign/"$Project

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

# Tmp file
tmp="/home/dyap/dyap_temp"

# QC for 300bp  (2x150bp PE MiSeq run)
# reads=150

infilesuffix="_isPCR-output.txt"
outfilesuffix="_primer_order.txt"

# Name of the input file
inputfile=$p3dir"/"$name$infilesuffix
# raw=$posdir"/primerIn-TNBC-"$type"-fix-List-altRefSeq-AccountGermMut.txt"
# raw=$posdir"/primerIn-TNBC-"$type"-fix-List-AccountGermMut.txt"
raw=$posdir"/"$name"_positions.txt"

# Name of the output file
orderfile=$p3dir"/"$name$outfilesuffix
manifestfile=$p3dir"/"$name".AmpliconManifest"
qcfile=$tmp"/"$name"_QC.txt"
manfile=$tmp"/"$name"_manifest"
ordfile=$tmp"/"$name"_order"

# Checking to see if required files exist in correct paths

if [ -f $inputfile ]
	then
	echo "Reading from file: "$inputfile
	else
	echo $inputfile": File not found. Exiting...."
	exit 1;
fi

if [ -f $raw ]
	then
	echo "Reading from file: "$raw
	else
	echo $raw": File not found. Exiting...."
	exit 1;
fi

if [ -f $ordfile ]
	then
	echo "Overwriting "$ordfile". Press ENTER to continue or Ctrl-C to EXIT."
	read ans
		echo "ID, Chr, Start, End, Amplicon Length, Left Primer, Left length, Right Primer, Right Length, UID" > $ordfile
fi

if [ -f $manifestfile ]
	then
	echo "Overwriting "$manifestfile". Press ENTER to continue or Ctrl-C to EXIT."
	read ans
		echo "[Header]" > $manfile					
		echo $name"  Manifest Version,1" >> $manfile				
		echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
		echo "	" >> $manfile
		echo "[Regions]" >> $manfile					
		echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manfile
fi

if [ -f $qcfile ]
	then
	echo "Overwriting "$qcfile". Press ENTER to continue or Ctrl-C to EXIT."
	read ans
	rm -f $qcfile
fi

# Read one position at a time from the position file
# Match it with the first instance from the isPCR output file
# If conditions are good, accept it or else take the next best position

# This is the first primer pair
# IF fails, n=2 for the next primer pair
failed=0

for i in `cat $raw  | awk -F, '{print $2}' | tr -d '"' | sort -u`
	do
	echo $i
		# QC to test if grep is successful
		test=`grep $i $inputfile | wc -l`
		echo $test
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
		fi

	ampchr=`grep -m1 $i $inputfile | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`


	plusstrand=`grep -m1 $i $inputfile | grep [0-9]+[0-9] | wc -c`
	minusstrand=`grep -m1 $i $inputfile | grep [0-9]-[0-9] | wc -c`

	if [[ $plusstrand > 0 ]] && [[ $minusstrand == 0 ]];
		then
		ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
		ampend=`grep -m1 $i $inputfile |  awk -F"+" '{print $2}' | awk -F" " '{print $1}'`
	fi

	if [[ $plusstrand == 0 ]] && [[ $minusstrand > 0 ]];
		then
		ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
		ampend=`grep -m1 $i $inputfile |  awk -F"-" '{print $2}' | awk -F" " '{print $1}'`
	fi

	myid=`grep -m1 $i $inputfile | awk -F" " '{print $2}'`
	amplen=`grep -m1 $i $inputfile | awk -F" " '{print $3}' | sed 's/bp//'`
	leftpri=`grep -m1 $i $inputfile | awk -F" " '{print $4}'`
	rightpri=`grep -m1 $i $inputfile | awk -F" " '{print $5}'`

	lines=`echo "($amplen/51)+1" | bc`
	ampseq=`grep -m1 -A$lines $i $inputfile  | sed 's/^>.*$//' | tr -d "\n"`
	
	chr=`echo $myid | awk -F"_" '{print $1}' | sed 's/chr//'`
	pos=`echo $myid | awk -F"_" '{print $2}'`

	llim=`echo "$ampstart + $reads - 10 " | bc`
	rlim=`echo "$ampend - $reads + 10 " | bc`

	if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
		then
		Miseq="PASS"
		else
			for n in 2 3 4 5 6 7 8 9 10
			do
			ampchr=`grep -m$n $i $inputfile | tail -n1 | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`

			myid=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $2}'`
			amplen=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $3}' | sed 's/bp//'`
			leftpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $4}'`
			rightpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $5}'`

			lines=`echo "($amplen/51)+1" | bc`
			ampseq=`grep -m$n -A$lines $i $inputfile  | tail -n$lines | sed 's/^>.*$//' | tr -d "\n"`
	
			plusstrand=`grep -m$n $i $inputfile | tail -n1 | grep [0-9]+[0-9] | wc -c`
			minusstrand=`grep -m$n $i $inputfile | tail -n1 | grep [0-9]-[0-9] | wc -c`

	if [[ $plusstrand > 0 ]] && [[ $minusstrand == 0 ]];
		then
		ampstart=`grep -m$n $i $inputfile | awk -F":" '{print $2}' | tail -n1 | awk -F"+" '{print $1}'`
		ampend=`grep -m$n $i $inputfile |  awk -F"+" '{print $2}' | tail -n1 | awk -F" " '{print $1}'`
	fi

	if [[ $plusstrand == 0 ]] && [[ $minusstrand > 0 ]];
		then
		ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | tail -n1 | awk -F"-" '{print $1}'`
		ampend=`grep -m1 $i $inputfile |  awk -F"-" '{print $2}' | tail -n1 | awk -F" " '{print $1}'`
	fi

	chr=`echo $myid | awk -F"_" '{print $1}' | sed 's/chr//'`
	pos=`echo $myid | awk -F"_" '{print $2}'`

	# QC for 300bp  (2x150bp PE MiSeq run)
	reads=150
	llim=`echo "$ampstart + $reads - 10 " | bc`
	rlim=`echo "$ampend - $reads + 10 " | bc`
 
			if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
				then
				Miseq="PASS-"$n
				break
				else
				Miseq="FAIL-"$n
				failed=`echo "$failed + 1" | bc`
			fi
			done
				
	fi
 
	echo -e $myid"\t"$ampchr"\t"$ampstart"\t"$ampend"\t"$amplen"\t QC: "$reads"bp: "$Miseq >> $qcfile

	if [[ $Miseq =~ "PASS" ]]
		then
		# Generate manifest file of all that pass

		llen=`echo $leftpri | wc -c`
		leftlen=`echo "$llen - 1" | bc` 
		rlen=`echo $rightpri | wc -c` 
		rightlen=`echo "$rlen - 1" | bc` 
		echo $myid,"chr"$ampchr,$ampstart,$ampend,$leftlen,$rightlen" : isPCR= "$amplen"bp" >> $manfile

		# Generate order file of all that pass
		echo $myid, "chr"$ampchr, $ampstart, $ampend, $amplen" bp",$leftpri, $leftlen, $rightpri, $rightlen >> $ordfile

	fi

	done

cp $ordfile $orderfile
cat $manfile | tr "," "\t" > $manifestfile
	

echo "Number that failed and iterated:"
echo $failed

exit;


