#!/bin/sh
# Script to get the top hits of isPCR
# Generate an order file for those primers - done
# Generate Manifest - done
# Generate HTML view file - not done
# scp that to godel.cluster.bccrc.ca - not done

# IF run as pipeline (comment this out)
Project="TNBC"
type="indel"
# type="SNV"

# Required
name=$Project"-"$type

# Project Directory
dir="/home/dyap/Projects/"$Project
tmp="/home/dyap/dyap_temp/TNBC"

# positions
posdir=$dir"/positions"

# IF run as pipeline (comment this out)
raw=$posdir"/primerIn-TNBC-"$type"-fix-List-altRefSeq-AccountGermMut.txt"
#raw=$posdir"/primerIn-TNBC-"$type"-fix-List-AccountGermMut.txt"

# primer3 output
p3dir=$dir"/primer3"

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

if [ -f $orderfile ]
	then
	echo "Overwriting "$orderfile". Press ENTER to continue or Ctrl-C to EXIT."
	read ans
		echo "ID, Chr, Start, End, Amplicon Length, Left Primer, Left length, Right Primer, Right Length, UID" > $ordfile
	rm -f $outputfile
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
	rm -f $outputfile
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

for i in `cat $raw | awk -F"\t" '{print $1"_"$2"_"$3}'`
	do
	echo $i
		# QC to test if grep is successful
		test=`grep -m1 $i $inputfile | wc -l`
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
			grep -m1 $i $inputfile
		fi

	ampchr=`grep -m1 $i $inputfile | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`
	ampstart=`grep -m1 $i $inputfile | awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
	ampend=`grep -m1 $i $inputfile |  awk -F"+" '{print $2}' | awk -F" " '{print $1}'`

	myid=`grep -m1 $i $inputfile | awk -F" " '{print $2}'`
	amplen=`grep -m1 $i $inputfile | awk -F" " '{print $3}' | sed 's/bp//'`
	leftpri=`grep -m1 $i $inputfile | awk -F" " '{print $4}'`
	rightpri=`grep -m1 $i $inputfile | awk -F" " '{print $5}'`
	ampseq=`grep -m1 -A12 $i $inputfile  | sed 's/^>.*$//' | tr -d "\n"`
	
	sam=`echo $myid | awk -F"_" '{print $1}'`
	chr=`echo $myid | awk -F"_" '{print $2}'`
	pos=`echo $myid | awk -F"_" '{print $3}' | awk -F"-" '{print $1}'`

        uid=`grep -m1 $pos $raw | awk -F"\t" '{print $5}'` 

	# QC for 500bp  (2x250bp PE MiSeq run) OR
	# QC for 300bp  (2x150bp PE MiSeq run)
	reads=220
	llim=`echo "$ampstart + $reads - 12 " | bc`
	rlim=`echo "$ampend - $reads + 12 " | bc`
 
	if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
		then
		Miseq="PASS"
		else
			for n in 2 3 4 5 6 7 8 
			do
			ampchr=`grep -m$n $i $inputfile | tail -n1 | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`
			ampstart=`grep -m$n $i $inputfile | tail -n1 | awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
			ampend=`grep -m$n $i $inputfile |  tail -n1 | awk -F"+" '{print $2}' | awk -F" " '{print $1}'`

			myid=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $2}'`
			amplen=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $3}' | sed 's/bp//'`
			leftpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $4}'`
			rightpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $5}'`
			ampseq=`grep -m$n -A12 $i $inputfile  | tail -n1 | sed 's/^>.*$//' | tr -d "\n"`
	
			# QC  MiSeq run
			llim=`echo "$ampstart + $reads - 12 " | bc`
			rlim=`echo "$ampend - $reads + 12 " | bc`
 
			if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]] && [[ $amplen > $reads ]];
				then
					Miseq="PASS-"$n
					break
			fi

			if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
				then
					Miseq="Condi-PASS-"$n
						
				else
					Miseq="FAIL-"$n
					failed=`echo "$failed + 1" | bc`
				

			fi
			done
				
	fi
 
	echo -e $uid"\t"$myid"\t"$ampchr"\t"$ampstart"\t"$ampend"\t"$amplen"\t QC: "$reads"bp: "$Miseq >> $qcfile

	if [[ $Miseq =~ "PASS" ]] && [[ $UID != "ID" ]];
		then
		# Generate manifest file of all that pass

		llen=`echo $leftpri | wc -c`
		leftlen=`echo "$llen - 1" | bc` 
		rlen=`echo $rightpri | wc -c` 
		rightlen=`echo "$rlen - 1" | bc` 
		echo $myid,"chr"$ampchr,$ampstart,$ampend,$leftlen,$rightlen,$uid" : isPCR= "$amplen"bp" >> $manfile

		# Generate order file of all that pass
		echo $myid, "chr"$ampchr, $ampstart, $ampend, $amplen" bp",$leftpri, $leftlen, $rightpri, $rightlen, $uid >> $ordfile

	fi

	done

cp $ordfile $orderfile
cat $manfile | tr "," "\t" > $manifestfile
	

echo "Number that failed and iterated:"
echo $failed

fail=`cat $qcfile | grep FAIL | awk -F"\t" '{print $2}'`

for i in $fail
	do
	echo $i
	grep $i $inputfile
	done
exit;


