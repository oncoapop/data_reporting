#!/bin/bash
# (CONTINUE FROM RTN - part II for debugging purposes)

# Script to get the top hits of isPCR
# Generated from the commandline isPCR on beast
# Generate an order file for those primers - done
# Generate Manifest - done
# Generate HTML view file - done
# scp that to godel.cluster.bccrc.ca - done

# IF run as pipeline (comment this out)
# Project="eIF4A3"
# Required
# sample="eIF4A3_NMD"
# name=$sample
# basedir
# dir="/home/dyap/Projects/PrimerDesign/"$Project
# reads=150

#==================================================
# Project Directory
# dir is imported from pipeline script in the environment
htmloutpath="/home/dyap/public_html"

# positions
posdir=$dir"/positions"

# primer3 output
p3dir=$dir"/primer3"

# Annotation output
annodir=$dir"/Annotate"

# Tmp file
tmp="/home/dyap/temp"

echo "Running Primer_Summary version v3.4..."

# Formatting html
bold=$(tput bold)
normal=$(tput sgr0)

# from isPCR command line output
# fasta formatted
infilesuffix="_isPCR-output.fa"

outfilesuffix="_primer_order.txt"

# Name of the input file
inputfile=$p3dir"/"$name$infilesuffix
raw=$posdir"/"$name"_positions.txt"
designfile=$p3dir"/"$name"_p3_design.txt"
chkfile=$p3dir"/"$name"_isPCR-input"
genomefile=$p3dir"/"$name"_isPCR-genomic.fa"

# Name of the output file
failfile=$p3dir"/"$name"_fail_QC"
suppfile=$p3dir"/"$name"_SupplFig.csv"
orderfile=$p3dir"/"$name$outfilesuffix
manifestfile=$p3dir"/"$name".AmpliconManifest"
qcfile=$tmp"/"$name"_QC.txt"
manfile=$tmp"/"$name"_manifest"
ordfile=$tmp"/"$name"_order"
annofile=$annodir"/"$name"_Annotate.csv"
annofile2=$annodir"/"$name"_anno.txt"
annofile3=$annodir"/"$name".variant_function"
htmfile=$tmp"/"$name"_summary"
htmlfile=$annodir"/"$name"_summary.html"

# This file contains all selected positions in the html file (for selection later)
allfile=$p3dir"/short-listed_positions.txt"

rm -f $failfile"*"
rm $p3dir/set2.csv
rm -f $allfile

# 1. Primers pairs which amplify unique transcripts
# 2. Primer pairs which amplify transcripts which can be distinguished either by length
# 3. Primer pairs which amplify transcripts which can be distinguished either by sequence (TODO)

# Check this file
# This file contains all the correctly (identified by length) amplified transcripts by the respective primers
selfile=$p3dir"/selected.tmp"

first=1

for m in `cat $selfile | awk -F"," '{print $6}' | sort -u `
        do
		echo $bold$m$normal
                # Split each line into respective components
                gene=`grep "$m" $selfile | awk -F"," '{print $2}'`
                enst=`grep "$m" $selfile | awk -F"," '{print $1}'`
                left=`grep "$m" $selfile | awk -F"," '{print $3}'`
                right=`grep "$m" $selfile | awk -F"," '{print $4}'`
                amplen=`grep "$m" $selfile | awk -F"," '{print $5}'`
                coords=`grep "$m" $selfile | awk -F"," '{print $6}' | awk -F: '{print $2}' | tr "+" "," | tr "-" ","`
                pos=`grep "$m" $selfile | awk -F"," '{print $6}' | awk -F: '{print $2}'`

		amptag=$enst":"$pos
		gentag=$gene"_"$enst
		rightrev=`echo $right | awk  'BEGIN {
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
                        }' `		


		llen=`echo $left | wc -c`
                leftlen=`echo "$llen - 1" | bc`
                rlen=`echo $right | wc -c`
                rightlen=`echo "$rlen - 1" | bc`

		# Count how many transcripts are amplifed by the primer pair combination
                count=`grep "$enst" $inputfile | grep "$left" | grep "$right" | wc -l`
		echo $count
		lim=`echo "$count-1" | bc`
		lines=`echo "($amplen/51)+1" | bc`
		morelines=`echo "$lines+1" | bc`

			echo "================="

			# use of dynamic variable assignment to keep track 
                        for l in  $(eval echo "{$first..$count}")
                                do
                                match=`grep -m$l "$right" $inputfile | grep "$left" | tail -n1`

                                trans=`echo $match | awk -F" " '{print $1}' | tr -d ">"`
                                ENST=`echo $match | awk -F" " '{print $1}' | awk -F":" '{print $1}' | tr -d ">"`
                                len=`echo $match | awk -F" " '{print $3}' | tr -d "bp"`

				# declaration of dynamic variables
			 	declare "match$l=$match"
			 	declare "trans$l=$trans"
			 	declare "ENST$l=$ENST"
			 	declare "len$l=$len"

				done

			not=0
			
			#  In this section, we check each of the assigned lengths
                        for l in  $(eval echo "{$first..$count}")
                                do
				
			        var="len$l"

				# Condition where length matches correct amplicon length and it is the only one
				# These are the unique trancripts for the primer set
	                        if [[ "${!var}" == "$amplen" ]] && [[ "$count" == "1" ]]
        	                        then
						echo $gene,$enst,$left,$right,$coords >> $p3dir"/set2.csv"
						echo "Unique: "$gene,$enst,$left,$right,$coords
						note="Unique"

						# Order File Generation
						echo $gene"_"$enst","$left","$right >> $ordfile

						# Amplicon Manifest
						echo $gene"_"$enst,$enst,$coords,$leftlen,$rightlen,"unique" >> $manfile 

						# Supplemental Figure information
                				echo $gene"_"$enst","$enst","$coords","$left","$right","$amplen"bp" >> $suppfile

        					ampseq=`grep -A$lines -m1 "$amptag" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |   GREP_COLOR="1;31" grep --color=always $left | GREP_COLOR="1;31" grep --color=always $rightrev`
						amplength=`echo $ampseq | wc -c`


						if [[ $amplength == 1 ]]
                                			then
        					ampseq=`grep -A$morelines -m1 "$amptag" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |  GREP_COLOR="1;31" grep --color=always $left | GREP_COLOR="1;31" grep --color=always $rightrev`
						fi

		# Genomic testing
        					genome=`grep "$enst" $genomefile | grep "$left" | grep "$right"`  
						genlen=`echo $genome | wc -c`
						echo $enst,$left,$right,"("$rightrev")"

						if [[ $genlen == 1 ]]
                                			then
        					genomeinfo=`echo $left"  - no match on genomic - " $rightrev  |  GREP_COLOR="1;34" grep --color=always "$left" | GREP_COLOR="1;34" grep --color=always "$rightrev"`
							else
						genomeinfo=`echo $genome | GREP_COLOR="1;34" grep --color=always "$left" | GREP_COLOR="1;34" grep --color=always "$right"`
						fi

	                # Generate the HTML File for view of all that pass
                		echo $bold$amptag" ----------------------------------------------------   "$gene" : "$enst" -----------------------------------------------  "$pos " ("$amplen") - "$note$normal >> $htmfile"1.tmp"
                		echo  $ampseq  >> $htmfile"1.tmp"
                		echo  $genomeinfo  >> $htmfile"1.tmp"
				echo  $genomeinfo
#                echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  Amplicon only " >> $htmfile"1.tmp"
               echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $htmfile"1.tmp"
                echo  "  " >> $htmfile"1.tmp"
                echo  "  " >> $htmfile"1.tmp"
				
				echo $amptag >> $allfile

				fi

				# Condition where lengths does NOT match correct amplicon length
				# Counts the mismatches
	                        if [[ "${!var}" != "$amplen" ]]
					then
						not=`echo "$not+1" | bc`
						echo $not
				fi

				# By elimination, if all of the non matches are different amplicon lengths
				# And there is a correct assigment of amplicon match
				# This means that we are able to use this primer set to distinguish correct transcript.
        	                if  [[ "${!var}" != "$amplen" ]] && [[ "$not" == "$lim"  ]]
   		      	                then
						echo $gene,$enst,$left,$right,$coords >> $p3dir"/set2.csv"
						echo $gene,$enst,$left,$right,$coords
						note="Discriminate_by_length"
 
						# Order File Generation
						echo $gene"_"$enst","$left","$right >> $ordfile
				
						# Amplicon Manifest
						echo $gene"_"$enst,$enst,$coords,$leftlen,$rightlen,"discriminate_by_length" >> $manfile 

						# Supplemental Figure information
                				echo $gene"_"$enst","$enst","$coords","$left","$right","$amplen"bp" >> $suppfile

			ampseq=`grep -A$lines -m1 "$amptag" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |  GREP_COLOR="1;31" grep --color=always $left | GREP_COLOR="1;31" grep --color=always $rightrev`
						amplength=`echo $ampseq | wc -c`

						if [[ $amplength == 1 ]]
                                			then
        		ampseq=`grep -A$morelines -m1 "$amptag" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |  GREP_COLOR="1;31" grep --color=always $left | GREP_COLOR="1;31" grep --color=always $rightrev`
						fi

		# Genomic testing
        					genome=`grep "$enst" $genomefile | grep "$left" | grep "$right"`  
						genlen=`echo $genome | wc -c`
						echo $enst,$left,$right,"("$rightrev")"

						if [[ $genlen == 1 ]]
                                			then
        					genomeinfo=`echo $left"  - no match on genomic - " $rightrev  |  GREP_COLOR="1;34" grep --color=always "$left" | GREP_COLOR="1;34" grep --color=always "$rightrev"`
							else
						genomeinfo=`echo $genome | GREP_COLOR="1;34" grep --color=always "$left" | GREP_COLOR="1;34" grep --color=always "$right"`
						fi

	                # Generate the HTML File for view of all that pass
                		echo $bold$amptag" ----------------------------------------------------   "$gene" : "$enst" -----------------------------------------------  "$pos " ("$amplen") - "$note$normal >> $htmfile"1.tmp"
                		echo  $ampseq  >> $htmfile"1.tmp"
                		echo  $genomeinfo  >> $htmfile"1.tmp"
				echo  $genomeinfo
#                echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  Amplicon only " >> $htmfile"1.tmp"
#                echo "                                                                                  Other Amplicons (that differ by length) vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" >> $htmfile"1.tmp"
				echo $amptag >> $allfile

                        for p in  $(eval echo "{$first..$count}")
                                do
				var="len$p"
				nonlines=`echo "("${!var}"/51)+1" | bc`
				nonmorelines=`echo "$nonlines+1" | bc`
			
				var2="match$p"
	nonseq=`grep -A$nonlines -m1 "${!var2}" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |  GREP_COLOR="1;32" grep --color=always $left | GREP_COLOR="1;32" grep --color=always $rightrev`
				nonlength=`echo $nonseq | wc -c`

					if [[ $nonlength == 1 ]]
                                			then
        nonseq=`grep -A$nonmorelines -m1 "${!var2}" $inputfile  | sed 's/^>.*$//' | tr -d "\n"  |  GREP_COLOR="1;32" grep --color=always $left | GREP_COLOR="1;32" grep --color=always $rightrev`
					fi
                	
			 	var3="ENST$p"
					if [[ "${!var3}" != "$enst" ]]
						then
						echo  $nonseq " ( "${!var3}" - "${!var}" )" >> $htmfile"1.tmp"
					fi
				done

                echo "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> $htmfile"1.tmp"
                echo  "  " >> $htmfile"1.tmp"
                echo  "  " >> $htmfile"1.tmp"
        					
                	        fi

				done
        done

echo "DONE."
echo "#####"

cp $ordfile $orderfile
echo "Order file generated: " $orderfile

cat $manfile | tr "," "\t" > $manifestfile
echo "Manifest File generated: " $manifestfile

cat $htmfile"1.tmp" | /home/dyap/INSTALL/aha-master/aha  > $htmlfile
sed "s/stdin/$sample/g" $htmlfile
cat $failfile"1.tmp" | /home/dyap/INSTALL/aha-master/aha  > $failfile
sed "s/stdin/$sample/g" $failfile

echo "Amplicon Check File generated: " $htmlfile
mkdir $htmloutpath"/"$Project

mkdir $htmloutpath"/"$Project"/"$sample
cp $htmlfile $htmloutpath"/"$Project"/"$sample"/"$name"_summary.html"
cp $orderfile $htmloutpath"/"$Project"/"$sample"/"$name"_primer_order.txt"
cp $manifestfile $htmloutpath"/"$Project"/"$sample"/"$name".AmpliconManifest.txt"
cp $qcfile $htmloutpath"/"$Project"/"$sample"/"$name"_QC.txt"
cp $suppfile $htmloutpath"/"$Project"/"$sample"/"$name"_SuppleFigFile.csv"

echo "Number that failed and iterated:"
echo $failed

echo "Files can be found here:"
echo $htmloutpath"/"$Project"/"$sample

rsync -vr --progress $htmloutpath/$Project/$sample dyap@pleione.myseqtools.com:/var/www/html/output/$Project/

rm -f $p3dir"/*.tmp"

exit;

