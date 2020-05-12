#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

# This script takes the output of primer3 (commandline)
# appended in SAnnn_p3_output
# summaries all the top hits PRIMER_0 (zero) and PRIMER_1 (top two hits for example)
# into the file SAnnn_p3_summary

# Working Directory
# exported from CGI script
# export dir=$dir
# export tmp=$tmp
# export Project=$Project
# export name=$name
# export basedir=$basedir
# export posdir=$posdir
# export p3dir=$p3dir
# export anndir=$anndir

clear 
cd $dir

# Input name of INPUT FILE(s)
# $name is exported from v0_primer3.sh

fname=$tmp"/"$name"_p3_output"

# Name processing - Do not change

# Source and Output directories where working files are stored
sourcedir=$tmp
outdir=$p3dir

# contains the list of positions
infile=$posdir"/"$name
# contains the output of this script
outfile=$outdir"/"$name"_p3"
wtpos=$posdir"/"$name"_WT_positions.txt"

# For debugging
echo tmp=$tmpfile
echo in=$infile
echo out=$outfile
#exit;

rm -f $outdir*.tmp

echo "Sequences for this tumour sample: " $sample > $outfile"2.tmp"
echo "Check primers before ordering - use this file for insilico-PCR" >> $outfile"2.tmp"

echo "Sequences for this tumour sample: " $tmpfile > $outfile"3.tmp"
echo "DO NOT USE THIS FILE FOR ORDERING PRIMERS - for alignment use only" >> $outfile"3.tmp"

		for j in `cat $infile  | awk -F- '{print $1}' | tr ":" "_"`

		do
		echo .

# This tests to see if there are any primers returned, if not, the record is skipped
		test=`grep -A19 $j $fname | grep "PRIMER_PAIR_EXPLAIN="` 
		
			if [[ $test =~ "ok 0" ]];
			
				then
					echo $j "------------- SKIPPED.";
			
				else 
					echo $j >> $outfile"2.tmp";
					echo $j >> $outfile"3.tmp";
					echo $j "included.";
			fi
		
		grep -A2 $j $fname | grep "SEQUENCE_TEMPLATE=" >> $outfile"3.tmp"

# gets the reverse complement of RIGHT primer but NOT LEFT primer -> 3.tmp
		grep -A30 $j $fname | grep "_0_SEQUENCE=" | sed s'/PRIMER_RIGHT_0_SEQUENCE=//' | awk  'BEGIN {
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
		grep -A30 $j $fname | grep "_0_SEQUENCE=" >> $outfile"2.tmp"
		grep -A55 $j $fname | grep "_1_SEQUENCE=" >> $outfile"2.tmp"

# gets the reverse complement of RIGHT primer but NOT LEFT primer -> 3.tmp
                grep -A55 $j $fname | grep "_1_SEQUENCE=" | sed s'/PRIMER_RIGHT_1_SEQUENCE=//' | awk  'BEGIN {
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

echo Processing files
cat $wtpos | awk -F, '{print $2"_"$3,$4}' | tr -d '"' > $outfile"66.tmp"

echo Counterchecking with templates...
# 3.tmp contains the reverse complement of the primers (ie for correct alignment to template)
# 4.tmp contains the output - fwd and rev primers in color! (uses html format to view)
# 22.tmp contains the exact amplicon given by the primers

echo "Sequences for this tumour sample: "$sample > $outfile"4.tmp"

echo ";File contains amplicons" > $outfile"22.tmp"

readfile=$anndir"/"$name"_anno.txt"

                for k in `cat $infile  | awk -F- '{print $1}' | tr ":" "_"`

                do
		
			grep $k $readfile >> $outfile"4.tmp"

			grep $k $readfile

		        match=`grep $k $outfile"66.tmp" | awk -F" " '{print $2}'`

	               test=`grep -A19 $k $fname | grep "PRIMER_PAIR_EXPLAIN="`

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

					echo "Amp " $left$invseq$right > $outtemp

					echo "		Left(5'->3')	= " $left
					echo "		Right(3'->5')	= " $right
					echo "Left(5'->3')	= " $left >> $outfile"4.tmp"
					echo "Right(3'->5')	= " $right >> $outfile"4.tmp"

					
						if [ ! $left = " " ] && [ ! $right = " " ];
							then    {

 		 	               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=always $right | grep --color=always $match >> $outfile"4.tmp"
 			               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=always $right | grep --color=always $match
				       echo ">"$k"_0" >> $outfile"22.tmp"
									}

						fi

				       echo "AMPLICON:  "$left$invseq$right
				       echo $left$invseq$right >> $outfile"22.tmp" 
	

					right1=`grep -A5 $k $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*$//' | tr -d "\n"`
					left1=`grep -A4 $k $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*=//'`
					invseq1=`grep -A2 $k $outfile"3.tmp" | grep -P -o '(?<='$left1')[A-Z]*(?='$right1')'`

		                       	
					echo "Amp " $left1$invseq1$right1 > $outtemp



					 	if [[ $left1 != "" ]];
                		                	then {
                                        			                                   			
								echo "          	Left1(5'->3')   = " $left1
                                        			echo "Left1(5'->3')         = " $left1>> $outfile"4.tmp"
								lflag=1
								}
					 	fi
                        		
					
						if [[ $right1 != "" ]];
							then { 
								
								echo "		Right1(3'->5')	= " $right1
								echo "Right1(3'->5')	= " $right1 >> $outfile"4.tmp"
								rflag=1
								}
						fi

						
						if [[ $left1 = "" ]] && [[ $right1 = "" ]];
							then    {
								echo "here!"
								continue
								} else {

						                grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=always $right1 | grep --color=always $match >> $outfile"4.tmp"
                						grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=always $right1 | grep --color=always $match
			   				        echo ">"$k"_1" >> $outfile"22.tmp"


			 				        echo "AMPLICON:     "$left1$invseq1$right1
			 				        echo $left1$invseq1$right1 >> $outfile"22.tmp" 



								}
						fi

					}
			fi
			

	

		echo "----" >> $outfile"4.tmp"
		echo "#################" 
	
		echo " " 
		
	
			done;

cp $outfile"2.tmp" $outfile"_order.txt"

cp $outfile"22.tmp" $outfile"_amplicons.fa"

cat $outfile"4.tmp" | /home/dyap/INSTALL/aha-master/aha > $outfile"_primer_check.html"

rm -f $outdir"/"*.tmp

# Copy to web directory for viewing
html="/home/dyap/public_html"
cd $html
mkdir $Project
cd $Project
mkdir $sample
htmlout=$html"/"$Project"/"$sample

cp $outfile"_primer_check.html"  $htmlout"/"
cp $outfile"_amplicons.fa"  $htmlout"/"
cp $outfile"_order.txt"  $htmlout"/"
cp $outfile"_failed.txt" $htmlout"/"

exit;

