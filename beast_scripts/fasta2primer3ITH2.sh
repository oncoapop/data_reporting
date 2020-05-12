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
outfile="primer3_input"
file=$dir"/"$outfile

echo "Output file = "$file".txt"

# get the name of source file
source=$name"_design.txt"

cd $sourcedir
rm -f $dir"/"$outfile
rm -f $dir"/"$outfile.txt
rm -f $dir"/*.tmp"

echo "Preparing Prime3 input file..."
# Get the Sample_ID from the file

	for i in `cat $source | awk -F"\t" '{print $2}'` 
	do
	echo i=$i	
	# This section extracts the required information from the file
	# Sequence_ID must contain the Sample ID_Chr_start-end

	chr=`grep -m1 $i $source | awk -F"\t" '{print $2}' | awk -F":" '{print $1}' `
	sta=`grep -m1 $i $source | awk -F"\t" '{print $2}' | awk -F":" '{print $2}' | awk -F"-" '{print $1}' `
	end=`grep -m1 $i $source | awk -F"\t" '{print $2}' | awk -F":" '{print $2}' | awk -F"-" '{print $2}' `
	sam=`grep -m1 $i $source | awk -F"\t" '{print $1}' `
	range=`echo "$end - $sta" | bc`
echo end=$end
echo sta=$sta
echo range=$range
	if [[ $sta =~ $end ]];
		then	ID=$sam"_"$chr"_"$sta;
		else	ID=$sam"_"$chr"_"$sta"+"$range"bp";
	fi
	# Primer3 does not read any DNA sequence other than ATCG - clean up SNP masks with N
	seq=`grep $i $source | awk -F"\t" '{print $3}' |  tr -c 'ATCGatcg' 'N' `

                nchar=`echo $seq | wc -c`
                seqlen=`echo "($nchar-1)" | bc`

echo seq=$seq
echo seqlen=$seqlen

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


# IF you change the conditions here, please also change the reiteration conditions
# in the primer3ITH script!
 
echo "======================================="
	{ 
	echo "PRIMER_SEQUENCE_ID="$ID;
	echo "SEQUENCE_TEMPLATE="$seq;
        echo "SEQUENCE_TARGET="$target;
	echo "P3_COMMENT="$sam"_"$chr"-"$sta;
	echo "PRIMER_OPT_SIZE=22";
	echo "PRIMER_MIN_SIZE=18";
	echo "PRIMER_MAX_SIZE=25";
	echo "PRIMER_NUM_NS_ACCEPTED=0";
	echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190 130-220";
	echo "PRIMER_FILE_FLAG=0";
	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
	echo "PRIMER_MIN_TM=58.0";
	echo "PRIMER_OPT_TM=60.0";
	echo "PRIMER_MAX_TM=63.0";	
	echo "PRIMER_GC_CLAMP=2";
	echo "PRIMER_PAIR_NUM_RETURNED=3";
	echo "PRIMER_EXPLAIN_FLAG=1";
	echo "PRIMER_PRODUCT_OPT_SIZE=150";
	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
	echo "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0";
	echo "=";
	} >> $dir"/"$outfile
	
echo "################################"
	
	done

# Clean up comment for testing
cp $file $file".txt"

echo "Number of inputs to Primer3:"
grep ^= $file".txt" | wc -l

exit;
