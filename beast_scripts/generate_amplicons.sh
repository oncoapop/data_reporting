#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# This script takes the output of primer3 (commandline)
# appended in SAnnn_p3_output.txt (from Mplex_primer3.sh)
# summaries all the top hits PRIMER_0 (zero) and PRIMER_1 (top two hits for example)
# into the file SAnnn_p3_summary

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/"

clear 
cd $dir

# Input name of INPUT FILE(s)
# $name is exported from fasta2primer3.sh
# if this script is run independently,
# uncomment this and input name here
# name="SA029"

fname=$name"_p3_output.txt"

# Name processing - Do not change

# Source and Output directories where working files are stored
sourcedir=$dir"SNV/"
outdir=$dir"SNV/"

tmpfile=`echo $fname | awk -F_ '{print $1}'`
# contains the list of positions
infile=$sourcedir$tmpfile"_pos.txt"
# contains the output of this script
outfile=$outdir$tmpfile"_p3"

# For debugging
echo tmp=$tmpfile
echo in=$infile
echo out=$outfile
#exit;

rm -f $outdir*.tmp

echo "Sequences for this tumour sample: " $tmpfile > $outfile"2.tmp"
echo "Check primers before ordering - use this file for insilico-PCR" >> $outfile"2.tmp"

echo "Sequences for this tumour sample: " $tmpfile > $outfile"3.tmp"
echo "DO NOT USE THIS FILE FOR ORDERING PRIMERS - for alignment use only" >> $outfile"3.tmp"

		for j in `cat $infile  | awk -F- '{print $1}' | tr ":" "_" | tr "c" "C"`

		do
		echo .

# This tests to see if there are any primers returned, if not, the record is skipped
		test=`grep -A19 $j $sourcedir$fname | grep "PRIMER_PAIR_EXPLAIN="` 
		
			if [[ $test =~ "ok 0" ]];
			
				then
					echo $j "------------- SKIPPED.";
			
				else 
					echo $j >> $outfile"2.tmp";
					echo $j >> $outfile"3.tmp";
					echo $j "included.";
			fi
		
		grep -A2 $j $sourcedir$fname | grep "SEQUENCE_TEMPLATE=" >> $outfile"3.tmp"

# gets the reverse complement of RIGHT primer but NOT LEFT primer -> 3.tmp
		grep -A30 $j $sourcedir$fname | grep "_0_SEQUENCE=" | sed s'/PRIMER_RIGHT_0_SEQUENCE=//' | awk  'BEGIN {
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
  			}' >> $outfile"3.tmp"
		grep -A30 $j $sourcedir$fname | grep "_0_SEQUENCE=" >> $outfile"2.tmp"
		grep -A55 $j $sourcedir$fname | grep "_1_SEQUENCE=" >> $outfile"2.tmp"

# gets the reverse complement of RIGHT primer but NOT LEFT primer -> 3.tmp
                grep -A55 $j $sourcedir$fname | grep "_1_SEQUENCE=" | sed s'/PRIMER_RIGHT_1_SEQUENCE=//' | awk  'BEGIN {
                          j = n = split("A C G T", t)
                                for (i = 0; ++i <= n;)
                                map[t[i]] = t[j--]
                                        }
                        {
# leaves the left primer alone (just prints it as it is)
                                if (/LEFT/) print
                        else {
                                for (i = length; i; i--)
                                printf "%s", map[substr($0, i, 1)]
                                print x
                                }
                        }' >> $outfile"3.tmp"		
		echo "@" >> $outfile"2.tmp"
		echo "@" >> $outfile"3.tmp"
		done


	echo --------------------------------


echo Completed.

echo "Number of records in summary file:"
grep "^@" -c $outfile"2.tmp"
grep "^@" -c $outfile"3.tmp"


echo Counterchecking with templates...
# 3.tmp contains the reverse complement of the primers (ie for correct alignment to template)
# 4.tmp contains the output - fwd and rev primers in color! (uses html format to view)
# 22.tmp contains the exact amplicon given by the primers

echo "<HTML>" > $outfile"4.tmp"
echo "<H1>Sequences for this tumour sample: "$tmpfile" </H1>" >> $outfile"4.tmp"
echo "<H2><i>File checking only - Right primers are in REV COMPL for alignment not for PCR!</i></H2>" >> $outfile"4.tmp"
echo "<H4><i>There is a bug the <b>second time</b> the right primer appears in the sequence, it can be ignored!!</i></H2>" >> $outfile"4.tmp"

echo ";File contains amplicons" > $outfile"22.tmp"

                for k in `cat $infile  | awk -F- '{print $1}' | tr ":" "_" | tr "c" "C"`

                do
		
	                echo "<H2>"$k"</H2><p>" >> $outfile"4.tmp"
			echo $k

	               test=`grep -A19 $k $sourcedir$fname | grep "PRIMER_PAIR_EXPLAIN="`

			# IF the output does not have any primers, skip it
                        if [[ $test =~ "ok 0" ]];

                                then
                                        echo "-------------------------------- SKIPPED.";

                                else
				{
					right=`grep -A3 $k $outfile"3.tmp" | grep -A1 "_0_" | sed 's/^P.*$//' | tr -d "\n"`
					left=`grep -A2 $k $outfile"3.tmp" | grep -A1 "_0_" | sed 's/^P.*=//'`

# This special case gets the intervening sequence between the primers left and right primers but not the priemrs themselves!
					invseq=`grep -A2 $k $outfile"3.tmp" | grep -P -o '(?<='$left')[A-Z]*(?='$right')'`

					echo "		Left(5'->3')	= " $left
					echo "		Right(3'->5')	= " $right
					echo "<H3>Left(5'->3')	= " $left"</H3>" >> $outfile"4.tmp"
					echo "<H3>Right(3'->5')	= " $right"</H3><p>" >> $outfile"4.tmp"

 		 	               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=always $right >> $outfile"4.tmp"
 			               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=auto $right
				       echo ">"$k"_0_" >> $outfile"22.tmp"
				       echo "AMPLICON:  "$left$invseq$right
				       echo $left$invseq$right >> $outfile"22.tmp" 
	
					echo "</p>" >> $outfile"4.tmp"

					right1=`grep -A5 $k $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*$//' | tr -d "\n"`
					left1=`grep -A4 $k $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*=//'`
					invseq1=`grep -A2 $k $outfile"3.tmp" | grep -P -o '(?<='$left1')[A-Z]*(?='$right1')'`

		                       	
					 	if [[ $left1 != "" ]];
                		                	then {
                                        			                                   			
								echo "          	Left1(5'->3')   = " $left1
                                        			echo "<H3>Left1(5'->3')         = " $left1" </H3>" >> $outfile"4.tmp"
								lflag=1
								}
					 	fi
                        		
					
						if [[ $right1 != "" ]];
							then { 
								
								echo "		Right1(3'->5')	= " $right1
								echo "<H3>Right1(3'->5')	= " $right1" </H3><p>" >> $outfile"4.tmp"
								rflag=1
								}
						fi

						
						if [[ $left1 != "" && $right != "" ]];
							then    {
						                grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=always $right1 >> $outfile"4.tmp"
                						grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=auto $right1
								echo ">"$k"_1_" >> $outfile"22.tmp"
			 				        echo "AMPLICON:     "$left1$invseq1$right1
			 				        echo $left1$invseq1$right1 >> $outfile"22.tmp" 
								}
						fi

					}
			fi
			

	

		echo "</p><HR>" >> $outfile"4.tmp"
		echo "#################" 
	
		echo " " 
		echo "<HR>" >> $outfile"4.tmp"
		
	
			done;

cp $outfile"2.tmp" $outfile"_order.txt"

cp $outfile"22.tmp" $outfile"_amplicons.fa"
cat $outfile"4.tmp" | tr -d "\^\[" > $outfile"5.tmp"
cat $outfile"5.tmp" | sed 's/00m/\<\/B\>/g' > $outfile"6.tmp"
cat $outfile"6.tmp" | sed 's/01\;31m/\<B\>/g' > $outfile"7.tmp"
cat $outfile"7.tmp" | tr -d "\^\[K" | tr -d "^[" > $outfile"_check.html"


rm -f $outdir*.tmp

exit;

