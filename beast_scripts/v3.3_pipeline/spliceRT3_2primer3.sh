#!/bin/sh

# This is the script to read the files from the splicing project
# http://wiki.mo.bccrc.ca/x/QaBC

temp="/home/dyap/dyap_temp"
indir="/home/dyap/Projects/PrimerDesign/manual"
outdir="/home/dyap/Projects/PrimerDesign/manual"
tempfile="/home/dyap/dyap_temp/newfactor4CG_input.temp"

################
name="Splice-RT3"
###############

outfile=$outdir"/"$name"_primer3_output"
infile=$indir"/"$name"_primer3_input.txt"

# This module selects only a subset of the cluster ID to work with.
# If you want to use the all the results comment out this section

	tfas=$temp"/"$name"_primer3_input.fa"

rm -f $tfas
rm -f $tempfile
rm -f $infile
rm -f $outfile

#<<<-------SELECTED SVs--------->>>
# This is the file has we have used to short-list the ones that we want to work on
dir="/home/dyap/Projects/PrimerDesign/manual"
sel=$dir"/new-factors_seq.txt"
sltd=$dir"/new-factors_seq_formatted.txt"


# This is in MAC format so we have to translate that into unix format
# awk '{ gsub("\r", "\n"); print $0;}' $sel > $sltd
# echo "converted"

#format: of $sel
# fasta

# To extract only the primer_id (since sample primer set for all samples used)
# cat $sel | awk '{ gsub("\r", "\n"); print $0;}' | awk -F"\t" '{print ","$1"_"$3"_"$4","$4}' | tr -d " " | sort -u > $sltd

# This concaternates the seq on the second line to the rest of the information on the first line
awk 'NR%2{printf $0" ";next;}1' $sel | awk -F"_" '{print ","$1","$2","$3}' | tr -d ">"  > $sltd


for i in `cat $sltd |  awk -F, '{print $2}' | sort -u ` 
	do

	if [[ $i == "___" ]];
		then 	{
		echo "Skipping header..."
		continue
			 }
	fi

        last=`grep $i $sltd | wc -l`
        first="1"

        for j in  $(eval echo "{$first..$last}")
                do

		# String (position) Field
		# $1	Sample
		# $2	gene_transcriptID_proteinID (non unique)
		# $3	sequence

		#	sa=`grep $i $sltd | awk -F, '{print $1}'` 
		# all primers to be run on all samples 

		id=`grep -m$j $i $sltd | tail -n1 | awk -F, '{print $2}'` 
		seq=`grep -m$j $i $sltd | tail -n1 | awk -F" " '{print $2}'` 
		seqlen=`echo $seq | wc -c` 

		if [[ $seq =~ "NA" ]];
			then 	{
				echo "Skipping failed position..."
				continue
			 	}
		fi

# Version 1.1 (with manual reiterations)
# The domain enriched in CG formation that is coded from is enclosed by the "[ ]"
# This counts the first instance for any other character other than A-Z which happens to be `[`
# Hence this leftdist represents the target start
# right is where the target sequence ends
		win=10
		leftdist=`echo $seq | grep -P -o '(?<='$match')[a-zA-Z]*(?='$match')' | wc -c`
		right=`echo "$leftdist + $win" | bc`

		cxt=`echo $seq | cut -c "$leftdist"-"$right"`


# New case: IF $len > $cutoff then I would want to use ONLY The region within the brackets for my design space!
	cutoff=1000
# The min is the min length of the seq in brackets above which we can only ensure that the 5' junction is targeted for primer design
# if it above min, amplicons will all fail because they are too long.
	min=100

# The window length does not give enough design space then only use first junction.
#	if	[[ $win -gt $min ]]
#		then
#		win=20
#	fi
	
		echo $i,$leftdist,$right,$len

 		echo $i","$seqlen",======,"$leftdist","$right","$len >> $tempfile


	if [[ -z "$cxt" ]];
		then 	{
		echo "Skipping... No context"
		continue
			 }
	fi

#	if [[ $len -gt $cutoff ]]
#		then
		{
		echo "PRIMER_SEQUENCE_ID="$id; 
		echo "SEQUENCE_TEMPLATE="$seq ;
		echo "P3_COMMENT="$cxt;
		echo "PRIMER_OPT_SIZE=21";
        	echo "PRIMER_MIN_SIZE=18";
        	echo "PRIMER_MAX_SIZE=26";
        	echo "PRIMER_NUM_NS_ACCEPTED=0";
        	echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-180"
        	echo "PRIMER_FILE_FLAG=0";
        	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
        	echo "PRIMER_MIN_TM=58.0";
        	echo "PRIMER_OPT_TM=60.0";
        	echo "PRIMER_MAX_TM=62.0";
        	echo "PRIMER_GC_CLAMP=2";
        	echo "PRIMER_NUM_RETURN=1";
		echo "PRIMER_MAX_POLY_X=4";
        	echo "PRIMER_EXPLAIN_FLAG=1";
        	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
        	echo "=";
        	} >> $tfas
#		else
#		{
#		echo "PRIMER_SEQUENCE_ID="$id; 
#		echo "SEQUENCE_TEMPLATE="$seq ;
#		echo "COMMENT="$leftdist","$right","$seqlen;
#		echo "P3_COMMENT="$cxt;
#		echo "SEQUENCE_TARGET="$leftdist","$win;
#		echo "PRIMER_OPT_SIZE=21";
 #       	echo "PRIMER_MIN_SIZE=18";
  #      	echo "PRIMER_MAX_SIZE=26";
   #     	echo "PRIMER_NUM_NS_ACCEPTED=0";
#        	echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-180"
 #       	echo "PRIMER_FILE_FLAG=0";
  #      	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
   #     	echo "PRIMER_MIN_TM=58.0";
    #    	echo "PRIMER_OPT_TM=60.0";
     #   	echo "PRIMER_MAX_TM=62.0";
#        	echo "PRIMER_GC_CLAMP=2";
 #       	echo "PRIMER_NUM_RETURN=1";
#		echo "PRIMER_MAX_POLY_X=4";
#        	echo "PRIMER_EXPLAIN_FLAG=1";
#        	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
#        	echo "=";
#        	} >> $tfas

#	fi

	done

	done
echo "Completed" 

#<<<-------SELECTED SVs--------->>>

# clean up to make unrecognised characters and spaces

cat $tfas  | tr -d "[]" | sed 's/=\ /=/' > $infile

echo "Running Primer3..."

cd $outdir

#primer3_core -format_output -output=primer3_view.txt < $infile
primer3_core -output=$outfile < $infile

echo "done."

echo "Press enter to continue..."
read ans

clear

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $infile

echo "Number of outputs in " $outfile
grep "^=" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile


#echo "Press Enter to show..."

#read ans
# grep -B26 "PAIR_NUM_RETURNED=0" $outfile | more
grep -B26 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

echo "Next iteration..."

###################################
# reiteration module goes in here #
###################################

# Part of the redesign pipeline
rename=$temp"/"$name"_p3_redesign.txt"

cat $outfile".failed" | sed 's/PRIMER_GC_CLAMP=2/PRIMER_GC_CLAMP=0/' | sed 's/140-180/140-180 130-190 120-200/'  | sed 's/SEQUENCE_TARGET=*.*$//' | sed '/^$/d' | grep -A18 "PRIMER_SEQUENCE_ID=" | sed 's/--/=/' > $rename

echo "=" >> $rename

reoutfile=$temp"/"$name"_p3_output2"
primer3_core -output=$reoutfile  < $rename

echo "Number of sequences in input file"
grep "SEQUENCE_TEMPLATE=" -c $rename

echo "Number of outputs in " $reoutfile
grep "^=" -c $reoutfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $reoutfile

grep -B25 "PAIR_NUM_RETURNED=0" $reoutfile > $reoutfile.failed

######################################

cat $outfile > $outfile".txt"
cat $reoutfile >> $outfile".txt"

echo "Total failed"
grep -c "PAIR_NUM_RETURNED=0" $reoutfile.failed
total=`grep -c "PAIR_NUM_RETURNED=0" $outfile`

exit;
	
