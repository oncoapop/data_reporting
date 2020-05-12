#!/bin/sh

# Script to demultiplex nuclear barcodes from single cell experiments
# These barcodes are at the beginning of the reads

# Directory
	dir="/home/dyap/Projects/Single_Cell/Nuc_Barcode"

# Run
	run="A8YFD"

# Working Directories
	wd=$dir"/"$run
	mkdir $wd"/output"
	mkdir $wd"/output/demux2"
	out="/home/dyap/Projects/Single_Cell/Nuc_Barcode/"$run"/output"
	demux="/home/dyap/Projects/Single_Cell/Nuc_Barcode/"$run"/output/demux2"


# Making a list of sample names
	cd $wd
	list=`ls *.fastq | awk -F_ '{print $1"_"$2}' | uniq`

# Command to combine pair-end reads into a single fastq
	cmd="/home/dyap/INSTALL/FLASH-1.2.10/flash"

	for i in $list
		do
		cd $wd
		names=`ls $i* | grep fastq`
		echo $cmd $names
		$cmd -t 8 --output-directory=$out --output-prefix=$i --max-overlap=200 $names
		done

# These are locus specific sequences
# Chr11_104972190	ATCACG
# Chr12_62104219	CGATGT
# Chr13_110844721	TTAGGC
# Chr17_26695832	TGACCA
# Chr19_19446936	GGCTAC
# Chr20_41100774	GCCAAT
# Chr22_18655860	CTTGTA
# Chr6_32188383		ACTTGA

LOC="ATCACG CGATGT TTAGGC TGACCA GGCTAC GCCAAT CTTGTA ACTTGA"

for Lo in $LOC
	do

	# These are called "nucodes" - which specify nuclei
	C1=$Lo"AAGCTA"
	C2=$Lo"GTATAG"
	C3=$Lo"GGAACT"
	C4=$Lo"CCTTGC"
	C5=$Lo"AGCATC"
	C6=$Lo"TCTGAG"
	C7=$Lo"CGGCCT"
	C8=$Lo"GAATGA"

	# Illlumina Adaptors (5'->3')
	fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
	ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
	rar="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

	# short adaptors in new strategy
	# fas="GTCAGATGTGTATAAGAGACAG"


	# use these commands to put files in wd
	# cp *.fastq.gz /home/dyap/Projects/Single_Cell/Nuc_Barcode/$run
	# gunzip *.fastq.gz


	# Summary file
	sumfile=$out"/summary_"$Lo".txt"

	echo "Locus = (Nucodes in order)"$C1","$C2","$C3","$C4","$C5","$C6","$C7","$C8 > $sumfile

	for i in $list

		do
		cd $out
		B1=`grep ^$C1 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C1 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B1.fastq"
		B2=`grep ^$C2 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C2 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B2.fastq"
		B3=`grep ^$C3 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C3 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B3.fastq"
		B4=`grep ^$C4 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C4 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B4.fastq"
		B5=`grep ^$C5 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C5 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B5.fastq"
		B6=`grep ^$C6 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C6 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B6.fastq"
		B7=`grep ^$C7 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C7 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B7.fastq"
		B8=`grep ^$C8 $i.extendedFrags.fastq | wc -l`
		grep -B1 -A2 ^$C8 $i.extendedFrags.fastq | sed '/--/d' > $demux"/"$i"_B8.fastq"

		echo $i" = "$B1","$B2","$B3","$B4","$B5","$B6","$B7","$B8 >> $sumfile

		done

	echo $Lo " complete."
	done
