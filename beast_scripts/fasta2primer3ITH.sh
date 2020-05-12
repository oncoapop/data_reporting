#!/bin/sh

# This Script was written by Dr Damian Yap (Apr 2013)
# updated Jul 2013 for ITH pipeline V2.0

# Primer3 is installead on beast
# and needs a specific input
# This script takes the output of GetRanges.R
# and parses it into primer3 input format

echo "This is the exported name" $name

# Working Directory
dir="/home/dyap/Projects/ITH/primer3"

# Source and Output directories where Barcoded files stored
sourcedir=$dir
outdir="/"
outfile="primer3_input"
file=$dir$outdir$outfile

echo "Output file = "$file".txt"

# get the name of source file
source=$name"_design.txt"


cd $sourcedir
rm -f $dir$outdir$outfile
rm -f $dir$outdir$outfile.txt
rm -f $dir$outdir*.tmp

echo "Preparing Prime3 input file..."
# Get the Sample_ID from the file

	for i in `cat $source | awk -F"," '{print $2}' | tr -d '"' ` 
	do
echo i=$i	
	# This section extracts the required information from the file
	# Sequence_ID must contain the Sample ID_Chr_start-end

	chr=`grep $i $source | awk -F, '{print $3}' | tr -d '"' `
	sta=`grep $i $source | awk -F, '{print $4}' | tr -d '"' `
	end=`grep $i $source | awk -F, '{print $5}' | tr -d '"' `
 	sam=`echo $i | awk -F"\t" '{print $1}'`
	range=`echo "$end - $sta" | bc`
echo range=$range
	if [[ $sta =~ $end ]];
		then	ID=$sam"_"$chr"_"$sta;
		else	ID=$sam"_"$chr"_"$sta"+"$range"bp";
	fi
	# Primer3 does not read any DNA sequence other than ATCG - clean up SNP masks with N
	seq=`grep $i $source | awk -F, '{print $8}' | tr -d '"' |  tr -c 'ATCGatcg' 'N' `

	seqlen=`echo $seq | wc -c`

echo seq=$seq,seqlen=$seqlen

		# The SNV or Indel start is always the middle position of seq
		# <start>,<length>
		width=20  
		# 10 bp on either side of SNV
		rstart=`echo "($seqlen/2)-($width/2)" | bc`
echo $rstart,$width

	        if [[ $range < "1" ]];
                then    target=$rstart","$width
                else    length=`echo "$width + $range" | bc ` 
			target=$rstart","$length
        	fi	

	if [[ $seq =~ "Design" ]];
		then	continue;
	fi


	{ 
	echo "PRIMER_SEQUENCE_ID="$ID;
	echo "SEQUENCE_TEMPLATE="$seq;
        echo "SEQUENCE_TARGET="$target;
	echo "P3_COMMENT="$sam"_"$chr"-"$sta;
	echo "PRIMER_OPT_SIZE=22";
	echo "PRIMER_MIN_SIZE=18";
	echo "PRIMER_MAX_SIZE=25";
	echo "PRIMER_NUM_NS_ACCEPTED=0";
	echo "PRIMER_PRODUCT_SIZE_RANGE=260-280 250-300";
	echo "PRIMER_FILE_FLAG=0";
	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
	echo "PRIMER_MIN_TM=58.0";
	echo "PRIMER_OPT_TM=60.0";
	echo "PRIMER_MAX_TM=62.0";	
	echo "PRIMER_GC_CLAMP=2";
	echo "PRIMER_PAIR_NUM_RETURNED=1";
	echo "PRIMER_EXPLAIN_FLAG=1";
	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
	echo "PRIMER_MISPRIMING_LIBRARY=/home/dyap/Projects/Genomes/RepBase18.06.fasta/primer3.ref.fa";
	echo "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0";
	echo "=";
	} >> $dir$outdir$outfile
	
echo "################################"
	
	done

# Clean up comment for testing
rm -f $dir$outdir/*.tmp
cp $file $file".txt"

echo "Number of inputs to Primer3:"
grep ^= $file".txt" | wc -l

exit;
