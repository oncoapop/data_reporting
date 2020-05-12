#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# This file takes the primers in fasta format and greps them from the amplicon file
# to generate the amplicon for PCR

# Working Directory
dir="/home/dyap/Projects/Tumour_Evol/Primer3_outputs/"

clear 
cd $dir

name="PE"

fname="primer3_"$name"_primer3_view.txt"

# Name processing - Do not change

# Source and Output directories where working files are stored
sourcedir=$dir
outdir=$dir

tmpfile=`echo $fname | awk -F_ '{print $2}'`
# contains the list of positions
infile=$sourcedir$tmpfile"_pos.txt"
# contains the output of this script
outfile=$outdir$tmpfile"_amplicons"

# For debugging
echo tmp=$tmpfile
echo in=$infile
echo out=$outfile

rm -f $outdir*.tmp
echo 

echo "Sequences for this tumour sample: " $tmpfile > $outfile"2.tmp"
echo "Starting..."
		for i in `grep "chr" $infile | awk -F- '{print $1}' | tr "c" "C" | tr ":" "_"`
		do
		echo $i

# Only reverse perimers are reversed
				grep -A1 $i $sourcedir$fname >> $outfile"2.tmp"


		grep -A28 $i $sourcedir$fname | grep "PRIMER_LEFT_0_SEQUENCE="  >> $outfile"2.tmp"
		
		grep -A28 $i $sourcedir$fname | grep "PRIMER_RIGHT_0_SEQUENCE=" | sed 's/^.*[=]//' | awk  'BEGIN {
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
  			}' >> $outfile"2.tmp"


		echo "@" >> $outfile"2.tmp"

		done


	echo --------------------------------


echo Completed.

echo "Number of records in summary file:"
grep "^@" -c $outfile"2.tmp"


echo Counterchecking with templates...
# 2.tmp contains the reverse complement of the primers (ie for correct alignment to template)

echo "Sequences for this tumour sample: " $tmpfile > $outfile"3.tmp"

		for k in `cat $sourcedir$fname | grep "Chr" | sed 's/^.*[=]//'`

               do
		
			echo $k

			record=`grep -A4 $k $outfile"2.tmp"`
			echo $record
			echo $k >> $outfile"3.tmp"
			left=`grep -A4 $k $outfile"2.tmp" | grep "PRIMER_LEFT_0_SEQUENCE=" | sed 's/^.*[=]//'`
			right=`grep -A4 $k $outfile"2.tmp" | grep -A1 "PRIMER_LEFT_0_SEQUENCE=" | sed 's/^PRIMER.*$//' | tr -d "\n"`

# This special case gets the intervening sequence between the primers left and right primers but not the primers themselves
					invseq=`grep -A1 $k $outfile"2.tmp" | grep "SEQUENCE_TEMPLATE=" | sed 's/SEQUENCE_TEMPLATE=//' | grep -P -o '(?<='$left')[A-Z]*(?='$right')'`
			echo "Left="$left"        Right="$right
			echo "Invseq="$invseq
			 				        echo "AMPLICON:     "$left$invseq$right
			 				        echo $left$invseq$right >> $outfile"3.tmp"
#			 				        echo $left | wc -c >> $outfile"3.tmp"
#			 				        echo $invseq | wc -c >> $outfile"3.tmp"
#			 				        echo $right | wc -c >> $outfile"3.tmp"
			 				        amplicon=`echo $left$invseq$right | wc -c`

								A=`echo $left$invseq$right | grep "A" -o | wc -l` 
								T=`echo $left$invseq$right | grep "T" -o | wc -l` 
								C=`echo $left$invseq$right | grep "C" -o | wc -l` 
								G=`echo $left$invseq$right | grep "G" -o | wc -l` 
							tm=`echo "(( 64.9 + (41 * ($C+$G)-(16.4))/($A+$T+$C+$G) ))" | bc`
						echo "Amplicon len= "$amplicon"         Tm ˚C= "$tm >> $outfile"3.tmp"
				echo $left","$tm >> $outfile"4.tmp"	 

		echo "#################" >> $outfile"3.tmp"
		echo "#################" 
	
	
			done;

cp $outfile"3.tmp" $outfile"_temp.txt"
cp $outfile"4.tmp" $outfile"_left_amplicontemp.txt"
 
exit;

