#!/bin/bash
# 1. Primers pairs which amplify unique transcripts
# 2. Primer pairs which amplify transcripts which can be distinguished either by length
# 3. Primer pairs which amplify transcripts which can be distinguished either by sequence (TODO)

# Check this file
infilesuffix="_isPCR-output.fa"
p3dir="/home/dyap/Projects/PrimerDesign/eIF4A3/primer3"
inputfile=$p3dir"/eIF4A3_NMD"$infilesuffix
selfile=$p3dir"/selected.tmp"
first=1
rm $p3dir/set.txt
rm $p3dir/set2.txt

for m in `cat $selfile | awk -F"," '{print $6}' | sort -u `
        do
		echo $m
                # create a unique pair pair name combination
                gene=`grep "$m" $selfile | awk -F"," '{print $2}'`
                enst=`grep "$m" $selfile | awk -F"," '{print $1}'`
                left=`grep "$m" $selfile | awk -F"," '{print $3}'`
                right=`grep "$m" $selfile | awk -F"," '{print $4}'`
                amplen=`grep "$m" $selfile | awk -F"," '{print $5}'`
                coords=`grep "$m" $selfile | awk -F"," '{print $6}' | awk -F: '{print $2}' | tr "+" "," | tr "-" ","`

                count=`grep "$enst" $inputfile | grep "$left" | grep "$right" | wc -l`
		echo $count
		lim=`echo "$count-1" | bc`

			echo "================="

                        for l in  $(eval echo "{$first..$count}")
                                do
                                match=`grep -m$l "$right" $inputfile | grep "$left" | tail -n1`

                                trans=`echo $match | awk -F" " '{print $1}' | tr -d ">"`
                                ENST=`echo $match | awk -F" " '{print $1}' | awk -F":" '{print $1}' | tr -d ">"`
                                len=`echo $match | awk -F" " '{print $3}' | tr -d "bp"`

			 	declare "match$l=$match"
			 	declare "trans$l=$trans"
			 	declare "ENST$l=$ENST"
			 	declare "len$l=$len"

				done

			not=0

                        for l in  $(eval echo "{$first..$count}")
                                do
				
			        var="len$l"

	                        if [[ "${!var}" == "$amplen" ]] && [[ "$count" == "1" ]]
        	                        then
						echo $gene,$enst,$left,$right,$coords >> $p3dir"/set2.txt"
						echo "Unique: "$gene,$enst,$left,$right,$coords 
				fi

	                        if [[ "${!var}" != "$amplen" ]]
					then
						not=`echo "$not+1" | bc`
						echo $not
				fi

        	                if  [[ "${!var}" != "$amplen" ]] && [[ "$not" == "$lim"  ]]
   		      	                then
						echo $gene,$enst,$left,$right,$coords >> $p3dir"/set2.txt"
						echo $gene,$enst,$left,$right,$coords 
				fi
                	        

				done
        done

