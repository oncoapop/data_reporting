#!/bin/sh

# This Script was written by Dr Damian Yap (Apr 2013)
# updated Jul 2013 for TNBC pipeline V2.0

# Primer3 is installed on beast
# and needs a specific input
# This script takes the output of GetRanges.R
# and parses it into primer3 input format

echo "This is the exported name" $name

# Working Directory
dir="/home/dyap/Projects/TNBC"

# Source and Output directories where Barcoded files stored
sourcedir=$dir"/positions"
outdir=$dir"/primer3/"
outfile="primer3_input"
file=$outdir$outfile

echo "Output file = "$file".txt"

# get the name of source file
	if [[ $name == "indel" ]];
		then 	source=$sourcedir"/primerIn-TNBC-indel-fix-List-altRefSeq-AccountGermMut.txt"
		else    source=$sourcedir"/primerIn-TNBC-SNV-fix-List-AccountGermMut.txt"
	fi


echo $source
read ans

rm -f $outdir$outfile
rm -f $outdir$outfile.txt
rm -f $outdir*.tmp

# Get the Sample_ID from the file

	for i in `cat $source | awk -F"\t" '{print $5}' | tr -d '"' ` 

	do
	echo $i
	# This section extracts the required information from the file
	# Sequence_ID must contain the Sample ID_Chr_start-end

	chr=`grep $i $source | awk -F"\t" '{print $2}' | tr -d '"' `
	sta=`grep $i $source | awk -F"\t" '{print $3}' | tr -d '"' `
	end=`grep $i $source | awk -F"\t" '{print $4}' | tr -d '"' `
	ref=`grep $i $source | awk -F"\t" '{print $6}' | tr -d '"' | wc -c`
	alt=`grep $i $source | awk -F"\t" '{print $7}' | tr -d '"' | wc -c`
 	sam=`echo $i | awk -F"-" '{print $1}'`
	range=`echo "$end - $sta" | bc`
	# echo $range
	if [[ $sta =~ $end ]];
		then	ID=$sam"_"$chr"_"$sta;
		else	
			if [[ $ref < $alt ]];
			# insertion
			then  
				ID=$sam"_"$chr"_"$sta"_insert_"$alt;
				range=$alt
				offset=0
			fi

			if [[ $ref > $alt ]];
			# deletion
			then  
				ID=$sam"_"$chr"_"$sta"_deletion_"$alt;
				range=0
				offset=$alt
			fi
 
	fi
	# Primer3 does not read any DNA sequence other than ATCG - clean up SNP masks with N
	if [[ $name == "indel" ]];
		then 	seq=`grep $i $source | awk -F"\t" '{print $10}' | tr -d '"' |  tr -c 'ATCGatcg' 'N' `;
		else 	seq=`grep $i $source | awk -F"\t" '{print $8}' | tr -d '"' |  tr -c 'ATCGatcg' 'N' `;
	fi
	
		# The SNV or Indel is startP, endP
		# For TNBC project, startP = mutationStartPos - 301
		# For TNBC project, endP =  mutationEndPos + 351
		# <start>,<width, where width=20, 10bp on either side>
		nchar=`echo $seq | wc -c`
		start=`echo "($nchar-351-10-$offset+$range)" | bc` 

	        if [[ $range < 1 ]];
                then    target=$start",20";
                else    length=`echo "$start + $range + 20" | bc ` ;
			target=$start","$length;
        	fi	

	if [[ $seq =~ "seqRegion" ]];
		then	continue;
	fi


	{ 
	echo "PRIMER_SEQUENCE_ID="$ID;
	echo "SEQUENCE_TEMPLATE="$seq;
        echo "SEQUENCE_TARGET="$target;
	echo "P3_COMMENT="$ID
	echo "PRIMER_OPT_SIZE=22";
	echo "PRIMER_MIN_SIZE=18";
	echo "PRIMER_MAX_SIZE=34";
	echo "PRIMER_NUM_NS_ACCEPTED=0";
# For 2 x 150 bp PE MiSeq run
#	echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190";
# For 2 x 250 bp PE MiSeq run
	echo "PRIMER_PRODUCT_SIZE_RANGE=240-260 230-300 220-420";
	echo "PRIMER_FILE_FLAG=0";
	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
	echo "PRIMER_MIN_TM=58.0";
	echo "PRIMER_OPT_TM=60.0";
	echo "PRIMER_MAX_TM=62.0";	
	echo "PRIMER_GC_CLAMP=0";
	echo "PRIMER_NUM_RETURN=1";
	echo "PRIMER_MAX_POLY_X=5";
	echo "PRIMER_EXPLAIN_FLAG=1";
	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
	echo "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0";
	echo "=";
	} >> $outdir$outfile
	
	
	done

# Clean up comment for testing
rm -f $outdir/*.tmp
cp $file $file".txt"

echo "Number of inputs to Primer3:"
grep ^= $file".txt" | wc -l

exit;
