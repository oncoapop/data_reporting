#!/bin/sh

# Script to copy vcf into directory to be read
# input directory
#source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/"
source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/VariantCallingLogs"

# output directory
dest="/home/dyap/Projects/ctDNA/AT094/variants"

# Find somatic mutations
# vcf present on the FFPE but not in the Genome (G -WT, buffy coat)

cd $dest

rm germline*
rm somatic*

#######################
# Compare for somatic


# get positions from sample #24
# match with #1 -> germline
# no match in #1 -> somatic


for i in `seq 24 2 48`
	do
	
	for j in `grep "^chr" "FFPE"$i"_common" | awk -F"\t" '{print $2}'`

		do
		match=`grep -m1 $j "FFPE"$i"_common" | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'`
        	next=`echo "$i - 23" | bc`
		match2=`grep -m1 $j "Buffy"$next"_common" | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'`
		echo "next"$next
		echo "match2="$match2


                        if [[ $match == $match2 ]];

                               then
				echo "----------------------------------"
                                       echo $match2 " matches " $match;
					echo $match2 >> "germline_"$i"_"$next
                                else 
					echo "++++++++++++++++++++++++++++++++++"
                                        echo $match2 " does not match " $match;
					echo $match >> "somatic_"$i"_"$next
                        fi

		done


	done

for i in `ls somatic*`
	do
	cat $i | sort -u > $i.txt
	cat $i | sed 's/chr//' | awk -F"\t" 
/home/dyap/bin/ANNOVAR/annovar/annotate_variation.pl
	done
