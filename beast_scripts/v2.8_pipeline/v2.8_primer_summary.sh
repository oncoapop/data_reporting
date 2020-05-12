#!/bin/bash
# Script to get the top hits of isPCR
# Generated from the commandline isPCR on Pleiades
# Generate an order file for those primers - done
# Generate Manifest - done
# Generate HTML view file - done
# scp that to godel.cluster.bccrc.ca - not done

# IF run as pipeline (comment this out)
# Project="Tumour_Xenograft"
# Required
# name="SA999"
# basedir
# dir="/home/dyap/Projects/PrimerDesign/"$Project
# reads=150

#==================================================
# Project Directory
# dir is imported from pipeline script in the environment
htmloutpath="/home/dyap/public_html"

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

# Annotation output
annodir=$dir"/Annotate"

# Tmp file
tmp="/home/dyap/temp"

echo "Running Primer_Summary version 3.1..."

# from isPCR command line output
# fasta formatted
infilesuffix="_isPCR-output.fa"

# from isPCR automated web iterface
# already parsed
# preinfilesuffix="_isPCR_chk.csv"
# infilesuffix="_isPCR-output.txt"

outfilesuffix="_primer_order.txt"

# Name of the input file
inputfile=$p3dir"/"$name$infilesuffix
raw=$posdir"/"$name"_positions.txt"
designfile=$p3dir"/"$name"_p3_design.txt"

# Name of the output file
failfile=$p3dir"/"$name"_fail_QC"
suppfile=$p3dir"/"$name"_SupplFig.csv"
orderfile=$p3dir"/"$name$outfilesuffix
manifestfile=$p3dir"/"$name".AmpliconManifest"
qcfile=$tmp"/"$name"_QC.txt"
manfile=$tmp"/"$name"_manifest"
ordfile=$tmp"/"$name"_order"
annofile=$annodir"/"$name"_Annotate.csv"
annofile2=$annodir"/"$name"_anno.txt"
annofile3=$annodir"/"$name".variant_function"
htmfile=$tmp"/"$name"_summary"
htmlfile=$annodir"/"$name"_summary.html"

# Checking to see if required files exist in correct paths

rm -f $failfile"*"

if [ -f $inputfile ]
	then
	echo "Overwriting "$inputfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
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
	# read ans
fi
	echo "ID, Chr, Start, End, Amplicon Length, Left Primer, Left length, Right Primer, Right Length, Comments" > $ordfile


if [ -f $manifestfile ]
	then
	echo "Overwriting "$manifestfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
fi
		echo "[Header]" > $manfile
		echo $name"  Manifest Version,1" >> $manfile
		echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $manfile
		echo "	" >> $manfile
		echo "[Regions]" >> $manfile
		echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $manfile

if [ -f $htmlfile ]
	then
	echo "Overwriting "$htmlfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
fi
############ PASS FILE
	echo "Amplicons for "$Project" Project" > $htmfile"1.tmp"
	echo "Sequences for Project: "$name >> $htmfile"1.tmp"
	date >> $htmfile"1.tmp"
	echo " " >> $htmfile"1.tmp"
	echo "The First row is the USCS Amplicon with the primers in Green and the SNV context in blue." >> $htmfile"1.tmp"
	echo "The Second Row contains SNP/SNV/indel masking (red Ns). To help, the SNV context is red if there are Ns present and blue if there are no Ns present" >> $htmfile"1.tmp"
	echo " " >> $htmfile"1.tmp"

	echo "Amplicons generated by Primer3, checked by in-silico PCR, QC'ed for Miseq run of "$reads >> $htmfile"1.tmp"
	echo "ID ---------------------------------------------------------------  ANNOVAR annotation ----------------------------------------------------   UCSC coordinates (SNV)" >> $htmfile"1.tmp"
	echo "==========  LEFT PRIMER ====================================== context of SNV (red or blue) ================================================= RIGHT PRIMER =====" | GREP_COLOR="1;32" grep "PRIMER"  >> $htmfile"1.tmp"

	echo " " >> $htmfile"1.tmp"
	echo " " >> $htmfile"1.tmp"

############ FAIL FILE
	echo "############### FAILED ############## Amplicons for "$Project" Project" > $failfile"1.tmp"
	echo "Sequences for Project: "$name >> $failfile"1.tmp"
	date >> $failfile"1.tmp"
	echo " " >> $failfile"1.tmp"

	echo "#### FAILED ##### Amplicons failed QC'ed for Miseq run of "$reads >> $failfile"1.tmp"
	echo "ID ---------------------------------------------------------------  ANNOVAR annotation ----------------------------------------------------   UCSC coordinates (SNV)" >> $failfile"1.tmp"
	echo "================================================================== context of SNV (red) =========================================================================="  >> $failfile"1.tmp"

	echo " " >> $failfile"1.tmp"
	echo " " >> $failfile"1.tmp"



if [ -f $qcfile ]
	then
	echo "Overwriting "$qcfile". Press ENTER to continue or Ctrl-C to EXIT."
	# read ans
	rm -f $qcfile
fi
		echo "Amplicon ID,Chromosome,Amplicon Start,Amplicon End,Left Primer,Right Primer,Amplicon Length" > $suppfile


# Read one position at a time from the position file
# Match it with the first instance from the isPCR output file
# If conditions are good, accept it or else take the next best position

# This is the first primer pair
# IF fails, n=2 for the next primer pair
failed=0
echo "Head of primer file:"
head $raw

for i in `cat $raw  | awk -F"," '{print $4}' | tr -d '"' | sort -u`
	do
	echo "i="$i
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
			note="5 primers pairs - "$test" amplicons"
			else
			note=" "
		fi

	# This module is specific for command-line isPcr output

	ampchr=`grep -m1 $i $inputfile | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`


        plusstrand=`grep -m1 $i $inputfile | grep [0-9]+[0-9] | wc -c`
        minusstrand=`grep -m1 $i $inputfile | grep [0-9]-[0-9] | wc -c`

     		if [[ $plusstrand == 1 ]] && [[ $minusstrand == 1 ]];
                		then
				break
       		fi

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
        rightprir=`grep -m1 $i $inputfile | awk -F" " '{print $5}' | awk  'BEGIN {
                          j = n = split("A C G T", t)
                                for (i = 0; ++i <= n;)
                                map[t[i]] = t[j--]
                                        }
                        {
                                if (/LEFT/) print
                        else {
                                for (i = length; i; i--)
                                printf "%s", map[substr($0, i, 1)]
                                print x
                                }
                        }' `



        chr=`echo $myid | awk -F"_" '{print $(NF-1)}' | sed 's/chr//'`
        pos=`echo $myid | awk -F"_" '{print $NF}'`


	ann=`grep -m1 $pos $annofile3 | awk -F"\t" '{print $1}'`
	gene=`grep -m1 $pos $annofile3 | awk -F"\t" '{print $2}'`
	snv=`grep -m1 $pos $annofile3 | awk -F" " '{print $6}'`

	snvchk=`grep -m1 $i $designfile | awk -F"," '{print $6}' | tr -d '"'`
	idchk=`grep -m1 $i $designfile | awk -F"," '{print $2}' | tr -d '"'`
	cxt=`grep -m1 $i $designfile | awk -F"," '{print $7}' | tr -d '"'`

	lines=`echo "($amplen/51)+1" | bc`
        ampseq=`grep -m1 -A$lines $i $inputfile  | sed 's/^>.*$//' | tr -d "\n" | tr "a-z" "A-Z" | GREP_COLOR="1;34" grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`
	maskseq=`grep -m1 $i $designfile | awk -F, '{print $8}' | grep --color=always $cxt | grep --color=always N |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`

	amplen=`echo $ampseq | wc -c`
	masklen=`echo $maskseq | wc -c`

        		if [[ $amplen == 1 ]];
                		then
				moreline=`echo "$lines+1" | bc`
			        ampseq=`grep -m1 -A$moreline $i $inputfile  | sed 's/^>.*$//' | tr -d "\n" | tr "a-z" "A-Z" | GREP_COLOR="1;34" grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`
			fi

        		if [[ $masklen == 1 ]];
                		then
				maskseq=`grep -m1 $i $designfile | awk -F, '{print $8}' |  GREP_COLOR="1;34" grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`
				
        		fi

echo "==================================="
echo $ampseq
echo $maskseq

# QC for MiSeq run (specific to input)
# Exported from the $Miseq
#	reads=150
# assume width=20

echo "+++++++++++++++++++++++++++++++"

	llim=`echo "$ampstart + $reads - 10 " | bc`
	rlim=`echo "$ampend - $reads + 10 " | bc`

 
	if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
		then
		echo "PASS!"
		Miseq="PASS"
		else
			for n in 2 3 4 5
			do

			ampchr=`grep -m$n $i $inputfile | tail -n1 | awk -F">" '{print $2}' | awk -F":" '{print $1}' | sed 's/chr//'`

        		plusstrand=`grep -m$n $i $inputfile | tail -n1 | grep [0-9]+[0-9] | wc -c`
        		minusstrand=`grep -m$n $i $inputfile | tail -n1 | grep [0-9]-[0-9] | wc -c`

        		if [[ $plusstrand == 1 ]] && [[ $minusstrand == 1 ]];
                		then
				break
        		fi

         		if [[ $plusstrand > 0 ]] && [[ $minusstrand == 0 ]];
                		then
                		ampstart=`grep -m$n $i $inputfile | tail -n1 | awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
                		ampend=`grep -m$n $i $inputfile |  tail -n1 | awk -F"+" '{print $2}' | awk -F" " '{print $1}'`
        		fi

        		if [[ $plusstrand == 0 ]] && [[ $minusstrand > 0 ]];
                		then
                		ampstart=`grep -m$n $i $inputfile | tail -n1 | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
                		ampend=`grep -m$n $i $inputfile |  tail -n1 | awk -F"-" '{print $2}' | awk -F" " '{print $1}'`
        		fi


                        myid=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $2}'`
                        amplen=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $3}' | sed 's/bp//'`
                        leftpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $4}'`
                        rightpri=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $5}'`
        		rightprir=`grep -m$n $i $inputfile | tail -n1 | awk -F" " '{print $5}' | awk  'BEGIN {
                          j = n = split("A C G T", t)
                                for (i = 0; ++i <= n;)
                                map[t[i]] = t[j--]
                                        }
                        {
                                if (/LEFT/) print
                        else {
                                for (i = length; i; i--)
                                printf "%s", map[substr($0, i, 1)]
                                print x
                                }
                        }' `


			lines=`echo "($amplen/51)+1" | bc`
        		ampseq=`grep -m$n -A$lines $i $inputfile  | tail -n$lines | sed 's/^>.*$//' | tr -d "\n" | tr "a-z" "A-Z" | GREP_COLOR="1;34" grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`
			maskseq=`grep -m$n $i $designfile | tail -n1 | awk -F, '{print $8}' | grep --color=always $cxt | grep --color=always N |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`

			amplen=`echo $ampseq | wc -c`
			masklen=`echo $maskseq | wc -c`

        		if [[ $amplen == 1 ]];
                		then
				moreline=`echo "$lines+1" | bc`
			        ampseq=`grep -m1 -A$moreline $i $inputfile  | tail -n$moreline | sed 's/^>.*$//' | tr -d "\n" | tr "a-z" "A-Z" | grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`
			fi

        		if [[ $masklen == 1 ]];
                		then
				maskseq=`grep -m$n $i $designfile | tail -n1 | awk -F, '{print $8}' |  GREP_COLOR="1;40" grep --color=always $cxt |  GREP_COLOR="1;32" grep --color=always $leftpri | GREP_COLOR="1;32" grep --color=always $rightprir`

				
        		fi

			if [[ $chr == $ampchr ]] && [[ $llim > $pos ]] && [[ $rlim < $pos ]];
				then
				Miseq="PASS-"$n
				break
				else
				Miseq="FAIL-"$n

		echo "#######  " $Miseq "  ########  " $myid" ---------------------------------------   "$ann"-"$gene"   -----------------------------------------------    "$ampchr":"$ampstart"-"$ampend" ("$snv")" >> $failfile"1.tmp"
		test=`grep $i $designfile | awk -F, '{print $8}' | grep --color=always $cxt | grep --color=always N  | wc -c`
        		if [[ $test == 1 ]];
				then
				grep $i $designfile | awk -F, '{print $8}' | GREP_COLOR="1;32" grep --color=always $cxt  >> $failfile"1.tmp"
				else
				grep $i $designfile | awk -F, '{print $8}' | grep --color=always $cxt | grep --color=always N  >> $failfile"1.tmp"
			fi
		
		echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $failfile"1.tmp"
		echo  "######  Reason for QC FAIL #####   "$note"  ########### "  >> $failfile"1.tmp"	
		echo  "  " >> $failfile"1.tmp"
		echo  "  " >> $failfile"1.tmp"

		llen=`echo $leftpri | wc -c`
		leftlen=`echo "$llen - 1" | bc`
		rlen=`echo $rightpri | wc -c`
		rightlen=`echo "$rlen - 1" | bc`
		echo $myid,"chr"$ampchr,$ampstart,$ampend,$leftlen,$rightlen,$ann"-"$gene"-"$amplen"bp-"$idchk >> $failfile"_man.tmp"

		echo $myid,"chr"$ampchr,$ampstart,$ampend,$amplen" bp",$leftpri,$leftlen,$rightpri,$rightlen,$ann,$gene,$idchk >> $failfile"_ord.tmp"

		# Supplemental Figure information
		echo $myid",chr"$ampchr","$ampstart","$ampend","$leftpri","$rightpri","$amplen"bp" >> $failfile"_supfig.tmp"

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
		echo $myid,"chr"$ampchr,$ampstart,$ampend,$leftlen,$rightlen,$ann"-"$gene"-"$amplen"bp-"$idchk >> $manfile

		# Generate order file of all that pass
		echo $myid,"chr"$ampchr,$ampstart,$ampend,$amplen" bp",$leftpri,$leftlen,$rightpri,$rightlen,$ann,$gene,$idchk >> $ordfile

		# Supplemental Figure information
		echo $myid",chr"$ampchr","$ampstart","$ampend","$leftpri","$rightpri","$amplen"bp" >> $suppfile

		# Generate the HTML File for view of all that pass
		echo $myid" ----------------------------------------------------   "$ann"-"$gene"   -----------------------------------------------    "$ampchr":"$ampstart"-"$ampend" ("$snv") - "$note >> $htmfile"1.tmp"
		echo  $ampseq  >> $htmfile"1.tmp"
		echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  Amplicon only ------------------------------------------------------------------------------------SNP/SNV/masked sequence vvvvvvvvvvvvvvvvvvvvvvvvvvvv _____________________________________________________" >> $htmfile"1.tmp"
		echo  $maskseq  >> $htmfile"1.tmp"

		echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $htmfile"1.tmp"
		echo  "  " >> $htmfile"1.tmp"
		echo  "  " >> $htmfile"1.tmp"

	fi

	done
echo "DONE."
echo "#####"

cp $ordfile $orderfile
echo "Order file generated: " $orderfile

cat $manfile | tr "," "\t" > $manifestfile
echo "Manifest File generated: " $manifestfile

cat $htmfile"1.tmp" | /home/dyap/INSTALL/aha-master/aha  > $htmlfile
sed "s/stdin/$sample/g" $htmlfile
cat $failfile"1.tmp" | /home/dyap/INSTALL/aha-master/aha  > $failfile
sed "s/stdin/$sample/g" $failfile

echo "Amplicon Check File generated: " $htmlfile
mkdir $htmloutpath"/"$Project

sample=$name

mkdir $htmloutpath"/"$Project"/"$sample
cp $htmlfile $htmloutpath"/"$Project"/"$sample"/"$name"_summary.html"
cp $orderfile $htmloutpath"/"$Project"/"$sample"/"$name"_primer_order.txt"
cp $manifestfile $htmloutpath"/"$Project"/"$sample"/"$name".AmpliconManifest.txt"
cp $qcfile $htmloutpath"/"$Project"/"$sample"/"$name"_QC.txt"
cp $suppfile $htmloutpath"/"$Project"/"$sample"/"$name"_SuppleFigFile.csv"

mkdir $htmloutpath"/"$Project"/"$sample"/failed"
cp $failfile $htmloutpath"/"$Project"/"$sample"/failed/"$name"_FAILED_positions.html"
cp $failfile"_man.tmp" $htmloutpath"/"$Project"/"$sample"/failed/"$name"_FAILED_manifest.txt"
cp $failfile"_ord.tmp" $htmloutpath"/"$Project"/"$sample"/failed/"$name"_FAILED_order.txt"
cp $failfile"_supfig.tmp" $htmloutpath"/"$Project"/"$sample"/failed/"$name"_FAILED_SuppFig.txt"

cat $htmloutpath"/"$Project"/"$sample"/failed/"$name"_FAILED_order.txt" | awk -F, '{print $1}' | awk -v sample="$sample" -F_ '{print sample","$(NF-1)":"$NF"-"$NF}' | 
sort -u > $htmloutpath"/"$Project"/"$sample"/failed/"$name"_pos-to-redesign-primers.txt"

echo "Number that failed and iterated:"
echo $failed

echo "Files can be found here:"
echo $htmloutpath"/"$Project"/"$sample
echo " or online at http://godel.cluster.bccrc.ca/workflow/primer3check/v3_pipeline/"$Project"/"$sample


rsync -vr --progress $htmloutpath/$Project/$sample dyap@godel.cluster.bccrc.ca:/var/www/html/workflow/primer3check/v3_pipeline/$Project/

rm -f $p3dir"/*.tmp"

exit;


