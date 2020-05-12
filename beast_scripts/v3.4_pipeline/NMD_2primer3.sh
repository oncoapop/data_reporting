#!/bin/sh

# This is the script to read the files from the eIF4A3 NMD project
#      Input is from /home/dyap/Projects/eIF4A3/design.csv
#        Generated manually using the R-script, Plot_NMD-compare.R
#        Format of csv: (no header)
#POLR1A_ENST00000409024,ATTCCTCCCGAATTTCAGAGGCAGAGGGATCG...........ATGTTCTTGGAGATCAACAT,127

temp="/home/dyap/dyap_temp"
indir="/home/dyap/Projects/PrimerDesign/manual"
outdir="/home/dyap/Projects/PrimerDesign/manual"


################
name="eIF4A3_NMD"
# no of primers to design per target
nrprimers="10"
###############

tempfile="/home/dyap/dyap_temp/"$name"_input.temp"

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
# Or can be used for formatting of the files
dir="/home/dyap/Projects/PrimerDesign/manual"
sel=$dir"/"$name"_seq.txt"
cp /home/dyap/Projects/eIF4A3/design.csv $sel
sltd=$dir"/"$name"_seq_formatted.txt"


# This is in MAC format so we have to translate that into unix format
# awk '{ gsub("\r", "\n"); print $0;}' $sel > $sltd
# echo "converted"

#format: of $sel (no headers)
#POLR1A_ENST00000409024,ATTCCTCCCGAATTTCAGAGGCAGAGGGATCG...........ATGTTCTTGGAGATCAACAT,127

# To extract only the primer_id (since sample primer set for all samples used)
# cat $sel | awk '{ gsub("\r", "\n"); print $0;}' | awk -F"\t" '{print ","$1"_"$3"_"$4","$4}' | tr -d " " | sort -u > $sltd
# This concaternates the seq on the second line to the rest of the information on the first line
# awk 'NR%2{printf $0" ";next;}1' $sel | awk -F"\t" '{print ","$1"_"$4"_"$5"_"$6","$13}' | tr -d " " > $sltd
# This just outputs the correct columns from an already correctly formatted input file (minus the header information)
cat $sel |  awk -F"," '{print ","$1","$2","$3}' | tr -d " " > $sltd

for i in `cat $sltd |  awk -F, '{print $2}' | sort -u ` 
	do

	if [[ $i == "_" ]];
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
		# $2	gene_transcriptID (unique)
		# $3	sequence

		#	sa=`grep $i $sltd | awk -F, '{print $1}'` 
		# all primers to be run on all samples 

		id=`grep -m$j $i $sltd | tail -n1 | awk -F, '{print $2}'` 
		seq=`grep -m$j $i $sltd | tail -n1 | awk -F, '{print $3}'` 
		seqlen=`grep -m$j $i $sltd | tail -n1 | awk -F, '{print $4}'` 

		if [[ $seq =~ "NA" ]];
			then 	{
				echo "Skipping failed position..."
				continue
			 	}
		fi

# Version 1.1 (with manual reiterations)

		{
		echo "PRIMER_SEQUENCE_ID="$id; 
		echo "SEQUENCE_TEMPLATE="$seq ;
		echo "P3_COMMENT="$seqlen;
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
        	echo "PRIMER_NUM_RETURN="$nrprimers;
		echo "PRIMER_MAX_POLY_X=4";
		echo "PRIMER_MIN_THREE_PRIME_DISTANCE=10";
        	echo "PRIMER_EXPLAIN_FLAG=1";
        	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
        	echo "=";
        	} >> $tfas

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
grep "^=$" -c $outfile

echo ================================

echo "Number of failed sequences where there are no primers"
grep -c "PAIR_NUM_RETURNED=0" $outfile
nrfail=`grep -c "PAIR_NUM_RETURNED=0" $outfile`

rm -f $outfile.failed

echo "Next iteration..."

###################################
# reiteration module goes in here #
###################################

first=1
for nr in $(eval echo "{$first..$nrfail}")
	do
        # This gets each record between the first ^=$ and the next ^=$
	awk -v n=$nr '/^=$/{l++} l>n{exit} l==n' $outfile > $outfile.failed.tmp

        failid=`grep "SEQUENCE_ID=" $outfile.failed.tmp`
        failseq=`grep "SEQUENCE_TEMPLATE=" $outfile.failed.tmp`

		{
		echo $failid;
		echo $failseq;
        	echo "PRIMER_MIN_SIZE=17";
        	echo "PRIMER_MAX_SIZE=26";
        	echo "PRIMER_NUM_NS_ACCEPTED=0";
        	echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-180 120-200 100-220"
        	echo "PRIMER_FILE_FLAG=0";
        	echo "PRIMER_PICK_INTERNAL_OLIGO=0";
        	echo "PRIMER_MIN_TM=55.0";
        	echo "PRIMER_OPT_TM=60.0";
        	echo "PRIMER_MAX_TM=65.0";
        	echo "PRIMER_GC_CLAMP=0";
        	echo "PRIMER_NUM_RETURN="$nrprimers;
		echo "PRIMER_MIN_THREE_PRIME_DISTANCE=10";
        	echo "PRIMER_EXPLAIN_FLAG=1";
        	echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/apps/primer3/primer3-2.3.5/src/primer3_config/";
        	echo "=";
        	} >> $outfile.failed

	done


# Part of the redesign pipeline
rename=$temp"/"$name"_p3_redesign.txt"

cat $outfile".failed"  > $rename

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
	
