#!/bin/sh

# This is the script to read the files from the destruct direction and then generate primer3 design 
# space from them using the python script that is written by Andrew McPherson

indir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/destruct_results/"
outdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/destruct_results/"

temp="/home/dyap/dyap_temp"
cd $indir
rm -f nohup.out

cd $outdir
outfile="primer3_output.txt"
infile="primer3_input.txt"

clear
echo "Reinterating conditions and resubmitted to Primer3..."

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
	
