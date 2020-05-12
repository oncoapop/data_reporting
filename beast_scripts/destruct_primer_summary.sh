#!/bin/sh

# This Script was written by Damian Yap (Jun 2013)

# This script takes the output of primer3 (commandline)
# appended in primer3_output.txt (from destruct2primer3.sh)
# summaries all the top hits PRIMER_0 (zero) 
# into the file destruct_p3_summary

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements"

clear 
cd $dir

fname="primer3_output.txt"

# Source and Output directories where working files are stored
sourcedir=$dir"/destruct_results"
temp="/home/dyap/dyap_temp"
outdir=$temp

tmpfile="destruct_p3"

# Getting the cluster_id of the selected regions
# also get the breakpoint context for mapping

# contains the list of cluster_id (which are unique per SV)
infile=$temp"/SVs_selected.txt"

# contains the output of this script
outfile=$outdir"/"$tmpfile
wtpos=$temp"/selected_wt_pos.csv"
bkpos=$temp"/selected_bk_pos.csv"

rm -f $wtpos
rm -f $bkpos

# temp files for calling Amplicon calculation script
outtemp="/home/dyap/dyap_temp/Temp_cal-input"
intemp="/home/dyap/dyap_temp/Temp_cal-output"

# temp files for generating context and matches
fasta=$temp"/SV_primer3_input.fa"

countid=1
countseq=1

# Test to check for unique cluster_ids in $infile
countall=`cat $infile | awk -F, '{print $2}' | wc -l`
countuniq=`cat $infile | awk -F, '{print $2}' | sort -u | wc -l`

 if [[ $countall =~ $countuniq ]];

                                then
                                        echo $j "All Cluster IDs are unique";

                                else
                                        echo "Cluster IDs are not unique";
                                        echo "Exiting....";
                                        exit;
 fi


# makes the file with amplicon name and breakpoint context
for i in `cat $infile | awk -F, '{print $2}'`
	do
	# This gets the unique cluster_id (which matches 2 WT and one chimeric sequence)
	ids=`grep $i $fasta | grep ID= | sed 's/PRIMER_SEQUENCE_ID=//'`
	# This gets the 6bp surrounding the breakpoint marked by "[]" by python script
	seqs=`grep -A1 $i $fasta | grep -o -P '.{0,6}\[.{0,7}'` 

		array=($ids)
		ID1=`printf "%s\n" "${array[0]}"`
		ID2=`printf "%s\n" "${array[1]}"`
		ID3=`printf "%s\n" "${array[2]}"`
# echo These are the individual IDS $ID1, this is one $ID2, this is the final one $ID3

		array=($seqs)
		SEQ1=`printf "%s\n" "${array[0]}"`
		SEQ2=`printf "%s\n" "${array[1]}"`
		SEQ3=`printf "%s\n" "${array[2]}"`

# echo THese are the indiv seqs $SEQ1, another one $SEQ2, final one $SEQ3

		echo $ID1","$SEQ1 >> $bkpos
		echo $ID2","$SEQ2 >> $wtpos
		echo $ID3","$SEQ3 >> $wtpos

	done		

rm -f $outdir*.tmp

echo "Sequences for this tumour sample: " $tmpfile > $outfile"2.tmp"
echo "Check primers before ordering - use this file for insilico-PCR" >> $outfile"2.tmp"

echo "Sequences for this tumour sample: " $tmpfile > $outfile"3.tmp"
echo "DO NOT USE THIS FILE FOR ORDERING PRIMERS - for alignment use only" >> $outfile"3.tmp"

echo "----------NEXT-----------"


		for l in `grep "PRIMER_SEQUENCE_ID=" $sourcedir/$fname | sed 's/PRIMER_SEQUENCE_ID=//'`

		do
		echo .

	# This tests to see if there are any primers returned, if not, the record is skipped
		test=`grep -A19 $l $sourcedir/$fname | grep "PRIMER_PAIR_EXPLAIN="` 
		
			if [[ $test =~ "ok 0" ]];
			
				then
					echo $l "------------- SKIPPED.";
			
				else 
					echo $l >> $outfile"2.tmp";
					echo $l >> $outfile"3.tmp";
					echo $l "included.";
			fi
		

		grep -A2 $l $sourcedir/$fname | grep "SEQUENCE_TEMPLATE=" >> $outfile"3.tmp"


# gets the reverse complement of RIGHT primer but NOT LEFT primer -> 3.tmp
		grep -A30 $l $sourcedir/$fname | grep "_0_SEQUENCE=" | sed s'/PRIMER_RIGHT_0_SEQUENCE=//' | awk  'BEGIN {
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

# This is for ordering so no need to reverse complement the primers
		grep -A30 $l $sourcedir/$fname | grep "_0_SEQUENCE=" >> $outfile"2.tmp"

		echo "@" >> $outfile"2.tmp"
		echo "@" >> $outfile"3.tmp"
		done


	echo --------------------------------


echo Completed.

echo "Number of records in summary file:"
grep "^@" -c $outfile"2.tmp"
grep "^@" -c $outfile"3.tmp"

echo Processing files
cat $bkpos | sed 's/^,*$//' | sed '/^ *$/d' | tr -d "[" | tr -d "]" > $outfile"66.tmp"
cat $wtpos | sed 's/^,*$//' | sed '/^ *$/d' | tr -d "[" | tr -d "]" >> $outfile"66.tmp"

echo Counterchecking with templates...
# 3.tmp contains the reverse complement of the primers (ie for correct alignment to template)
# 4.tmp contains the output - fwd and rev primers in color! (uses html format to view)
# 22.tmp contains the exact amplicon given by the primers
# 5.tmp contains the makings of the manifest file
# 6.tmp is the custom genome

# Headers for manifest and html files

echo "[Header]" > $outfile"5.tmp"					
echo "Destruct Manifest Version	1" >> $outfile"5.tmp" 		
echo "ReferenceGenome	C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $outfile"5.tmp"
echo "[Regions]" >> $outfile"5.tmp"					
echo "Name	Chromosome	Amplicon Start	Amplicon End	Upstream Probe Length	Downstream Probe Length" >> $outfile"5.tmp"

echo "Custom genome for Destruct alignment"  > $outfile"6.tmp"

echo "<HTML>" > $outfile"4.tmp"
echo "<H1>Sequences for "$tmpfile" </H1>" >> $outfile"4.tmp"
echo "<H2><i>File checking only - Right primers are in REV COMPL for alignment not for PCR!</i></H2>" >> $outfile"4.tmp"

echo ";File contains amplicons" > $outfile"22.tmp"


		for m in `grep "PRIMER_SEQUENCE_ID=" $sourcedir/$fname | sed 's/PRIMER_SEQUENCE_ID=//'`

                do
	                echo "<HR><p>"$count"<p><H2>"$m"<H2><p>" >> $outfile"4.tmp"

			echo $m

		        match=`grep $m $outfile"66.tmp" | awk -F"," '{print $2}'`

	               test=`grep -A25 $m $sourcedir/$fname | grep "PRIMER_PAIR_EXPLAIN="`

			# IF the output does not have any primers or no SA number they are useless, skip them!
                        if [[ $test =~ "ok 0" && $m =~ "_SA" ]];

                                then
                                        echo "-------------------------------- SKIPPED.";

                                else
				{
					right=`grep -A3 $m $outfile"3.tmp" | grep -A1 "_0_" | sed 's/^P.*$//' | tr -d "\n"`
					left=`grep -A2 $m $outfile"3.tmp" | grep -A1 "_0_" | sed 's/^P.*=//'`

# This special case gets the intervening sequence between the primers left and right primers but not the priemrs themselves!
					invseq=`grep -A2 $m $outfile"3.tmp" | grep -P -o '(?<='$left')[A-Z]*(?='$right')'`

					echo "Amp " $left$invseq$right > $outtemp
					/home/dyap/Scripts/Temp_cal.sh
					atemp=`grep "Amp" $intemp | awk -F"_" '{print $2}'`
					aleng=`grep "Amp" $intemp | awk -F"_" '{print $3}'`

					echo "		Left(5'->3')	= " $left
					echo "		Right(3'->5')	= " $right
					echo "<H5>Left(5'->3')	= " $left"</H5>" >> $outfile"4.tmp"
					echo "<H5>Right(3'->5')	= " $right"</H5><p>" >> $outfile"4.tmp"

					echo "<PRE><p>" >> $outfile"4.tmp"

					 if [[ $left != "" && $right != "" && $match != "" ]];
                                                        then    {

 		 	               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=always $right | grep --color=always $match >> $outfile"4.tmp"
 			               grep -B1 --color=always $left $outfile"3.tmp" | grep --color=always $right | grep --color=always $match
				       echo ">"$k"_0_"$atemp"_"$aleng >> $outfile"22.tmp"

				       echo "AMPLICON:  "$left$invseq$right
				       echo $left$invseq$right >> $outfile"22.tmp" 
	
					echo "</PRE></p>" >> $outfile"4.tmp"
					echo "<p>Amplicon " $atemp" C   "$aleng" bp </p><hr>" >> $outfile"4.tmp"
					echo "Amplicon " $atemp" C   "$aleng" bp"
					len=`echo $aleng | sed s'/Length=//'`

					right1=`grep -A5 $m $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*$//' | tr -d "\n"`
					left1=`grep -A4 $m $outfile"3.tmp" | grep -A1 "_1_" | sed 's/^P.*=//'`
					invseq1=`grep -A2 $m $outfile"3.tmp" | grep -P -o '(?<='$left1')[A-Z]*(?='$right1')'`

		                       	
					echo "Amp " $left1$invseq1$right1 > $outtemp
					/home/dyap/Scripts/Temp_cal.sh
					atemp=`grep "Amp" $intemp | awk -F"_" '{print $2}'`
					aleng=`grep "Amp" $intemp | awk -F"_" '{print $3}'`

					# This generates the manifest file & "genome" for breakpoint only
 					if [[ $m =~ "chimeric" ]];
                                                       then {
								clussam=`echo $m | awk -F"-" '{print $1}'`
								chr=$clussam
								uppos=1
								downpos=$len
								upprim=`echo $left | wc -c`
								downprim=`echo $right | wc -c`
								echo $clussam"	"$chr"	"$uppos"	"$downpos"	"$upprim"	"$downprim >> $outfile"54.tmp"
								
								echo ">"$clussam":1-"$len >> $outfile"6.tmp"
								echo $left$invseq$right >> $outfile"6.tmp"
                                                                }
                                                fi


					
								}
					fi
					 	if [[ $left1 != "" ]];
                		                	then {
                                        			                                   			
								echo "          	Left1(5'->3')   = " $left1
                                        			echo "<H5>Left1(5'->3')         = " $left1" </H5>" >> $outfile"4.tmp"
								lflag=1
								}
					 	fi
                        		
					
						if [[ $right1 != "" ]];
							then { 
								
								echo "		Right1(3'->5')	= " $right1
								echo "<H5>Right1(3'->5')	= " $right1" </H5><p>" >> $outfile"4.tmp"
								rflag=1
								}
						fi

						
						if [[ $left1 != "" && $right1 != "" && $match != "" ]];
							then    {
			
								echo "<PRE><p>" >> $outfile"4.tmp"

						                grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=always $right1 | grep --color=always $match >> $outfile"4.tmp"
                						grep -B1 --color=always $left1 $outfile"3.tmp" | grep --color=always $right1 | grep --color=always $match
			   				        echo ">"$k"_1_"$atemp"_"$aleng >> $outfile"22.tmp"


			 				        echo "AMPLICON:     "$left1$invseq1$right1
			 				        echo $left1$invseq1$right1 >> $outfile"22.tmp" 

								echo "</PRE></p>" >> $outfile"4.tmp"

								echo "<p>Amplicon " $atemp" C   "$aleng" bp </p><hr>" >> $outfile"4.tmp"
								echo "Amplicon " $atemp" C   "$aleng" bp"

								}
						fi

					 count=`echo $count+1 | bc`
					}
			fi
			

	

		echo "</p><HR>" >> $outfile"4.tmp"
		echo "#################" 
	
		echo " " 
		echo "<HR>" >> $outfile"4.tmp"
		
	
			done;

cp $outfile"2.tmp" $outfile"_order.txt"
cp $outfile"_order.txt" $dir"/destruct_primer_order.txt"

# Generate the header and manifest and genome files for breakpoint data only
rm -f $dir"/Destruct.AmpliconManifest.txt"

cat $outfile"5.tmp" > $dir"/Destruct.AmpliconManifest.txt"
cat $outfile"54.tmp" >> $dir"/Destruct.AmpliconManifest.txt"
cat $outfile"6.tmp" > $dir"/Destruct_genome.fa"

cp $outfile"22.tmp" $outfile"_amplicons.fa"
cat $outfile"4.tmp" | tr -d "\^\[" > $outfile"5.tmp"
cat $outfile"5.tmp" | sed 's/00m/\<\/FONT\>/g' > $outfile"6.tmp"
cat $outfile"6.tmp" | sed 's/01\;31m/\<FONT\ COLOR=#FF0000\>/g' > $outfile"7.tmp"
cat $outfile"7.tmp" | tr -d "\^\[K" | tr -d "^[" > $outfile"8.tmp"

# This nifty command removes duplicate lines
awk '!x[$0]++' $outfile"8.tmp" > $outfile"_check.html"

rm -f $outdir*.tmp

godelpath="/var/www/html/workflow/primer3check/"

expname=$outfile"_check.html"

fname=`echo $expname| sed 's/^.*\///'`

scp $expname dyap@godel.cluster.bccrc.ca:$godelpath$fname

exit

