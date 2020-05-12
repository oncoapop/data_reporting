#!/bin/sh

# This is the script to read the files from the destruct direction and then generate primer3 design 
# space from them using the python script that is written by Andrew McPherson

indir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/destruct_results/"
outdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/destruct_results/"

python="/share/apps/install/bin/python"
script="/home/dyap/Scripts/primer_sequences.py"

genome="/share/lustre/backup/dyap/Projects/Genomes/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa"
offset="200"

temp="/home/dyap/dyap_temp"
cd $indir
rm -f nohup.out

for i in `ls SA* | sed 's/_.*$//'` 

	do

	echo $i
	file=$indir$i"_filtered_breakpoints.tsv"

	# call the python script for each file name

	cd $temp

	echo "Processing "$i

	nohup $python $script $genome $offset < $file

	cp nohup.out $outdir$i"_designspace.txt"

	echo "done."

	rm nohup.out

	done

	
# This module selects only a subset of the cluster ID to work with.
# If you want to use the all the results comment out this section

	tfas=$temp/SV_primer3_input.fa

rm -f $tfas

#<<<-------SELECTED SVs--------->>>
# This is the file has we have used to short-list the ones that we want to work on
sel="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/SVs_selected.csv"
sltd="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/SVs_selected.txt"


# This is in MAC format so we have to translate that into unix format
awk '{ gsub("\r", "\n"); print $0;}' $sel > $sltd

echo "converted"

for i in `cat $sltd |  awk -F, '{print $2}'` 
	do

	if [[ $i =~ "cluster_id" ]];
		then 	{
		echo "Skipping header..."
		continue
			 }
	fi

	echo $i


# String (position) Field
# $1	Sample
# $2	cluster_id (unique)
# $3	chromosome_1
# $5	start_1
# $6	end_1
# $8 	gene_name1
# $9	gene_location1
# $10	break_1
# $11	chromosome_2
# $13	start_2
# $14	end_2
# $16	gene_name2
# $17	gene_location2
# $18	break_2
# $24	type (of SV)
# $32	sequence

# 
	sa=`grep $i $sltd | awk -F, '{print $1}'` 
	id=`grep $i $sltd | awk -F, '{print $2}'` 

	cr1=`grep $i $sltd | awk -F, '{print $3}'` 
	s1=`grep $i $sltd | awk -F, '{print $5}'` 
	e1=`grep $i $sltd | awk -F, '{print $6}'` 
	gene1=`grep $i $sltd | awk -F, '{print $8}'` 
	loc1=`grep $i $sltd | awk -F, '{print $9}'` 
	bk1=`grep $i $sltd | awk -F, '{print $10}'` 

	cr2=`grep $i $sltd | awk -F, '{print $11}'` 
	s2=`grep $i $sltd | awk -F, '{print $13}'` 
	e2=`grep $i $sltd | awk -F, '{print $14}'` 
	gene2=`grep $i $sltd | awk -F, '{print $16}'` 
	loc2=`grep $i $sltd | awk -F, '{print $17}'` 

	type=`grep $i $sltd | awk -F, '{print $24}'` 
	seq=`grep $i $sltd | awk -F, '{print $32}'` 

	cd $outdir
	# The first seq is the chimeric seq that contains the breakpoint
	seq1=`grep $i *_designspace.txt | awk -F"\t" '{print $2}'`
	seq2=`grep $i *_designspace.txt | awk -F"\t" '{print $3}'`
	seq3=`grep $i *_designspace.txt | awk -F"\t" '{print $4}'`


# Input and output file for Primer3 
outfile="primer3_output.txt"
infile="primer3_input.txt"

# Version 1.1 (with manual reiterations)
# The breakpoint is denoted by the "[]" encoded by the python script (only first seq is chimeric)
# The breakpoint is in the middle, if offset=200 then it is position 201
# <start>,<length>

	{
	echo "PRIMER_SEQUENCE_ID="$id"_"$sa"_chimeric_"$type; 
	echo "SEQUENCE_TEMPLATE="$seq1 ;
	echo "P3_COMMENT="$id"_"$cr1":"$s1"-"$e1"@"$type$"@"$cr2":"$s2"="$e2 ;
        echo "SEQUENCE_TARGET=175,50";
	echo "PRIMER_OPT_SIZE=21";
        echo "PRIMER_MIN_SIZE=17";
        echo "PRIMER_MAX_SIZE=32";
        echo "PRIMER_NUM_NS_ACCEPTED=0";
        echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190"
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

	{
	echo "PRIMER_SEQUENCE_ID="$id"_"$sa"_WT_1";
	echo "SEQUENCE_TEMPLATE="$seq2;
	echo "P3_COMMENT="$id"_"$cr1":"$s1"-"$e1;
        echo "SEQUENCE_TARGET=175,50";
	echo "PRIMER_OPT_SIZE=21";
        echo "PRIMER_MIN_SIZE=17";
        echo "PRIMER_MAX_SIZE=32";
        echo "PRIMER_NUM_NS_ACCEPTED=0";
        echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190"
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

	{
	echo "PRIMER_SEQUENCE_ID="$id"_"$sa"_WT_2";
	echo "SEQUENCE_TEMPLATE="$seq3 ;
	echo "P3_COMMENT="$id"_"$cr2":"$s2"-"$e2; 
        echo "SEQUENCE_TARGET=175,50";
	echo "PRIMER_OPT_SIZE=21";
        echo "PRIMER_MIN_SIZE=17";
        echo "PRIMER_MAX_SIZE=32";
        echo "PRIMER_NUM_NS_ACCEPTED=0";
        echo "PRIMER_PRODUCT_SIZE_RANGE=150-170 140-190"
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

	
	done

echo "Completed" 

#<<<-------SELECTED SVs--------->>>

# clean up to make unrecognised characters and spaces

cat $tfas | tr -d "[]" | sed 's/=\ /=/' > $outdir$infile

echo "Running Primer3..."

cd $outdir

primer3_core -format_output -output=primer3_view.txt < $infile
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
#grep -B25 "PAIR_NUM_RETURNED=0" $outfile | more
grep -B25 "PAIR_NUM_RETURNED=0" $outfile > $outfile.failed

echo Press Return to continue...
read ans

echo "Total failed"
grep -c "PAIR_NUM_RETURNED=0" $outfile
total=`grep -c "PAIR_NUM_RETURNED=0" $outfile`

echo "WT-1"
grep "_1" $outfile.failed | grep ID | wc -l
wt1=`grep "_1" $outfile.failed | grep ID | wc -l`

echo "WT-2"
grep "_2" $outfile.failed | grep ID | wc -l
wt2=`grep "_2" $outfile.failed | grep ID | wc -l`

echo "Chimeric"
echo "( $total - ($wt1 + $wt2))" | bc


exit;
	
