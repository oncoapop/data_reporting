#!/bin/sh

# Scipt to map qPCR plate with primer plate and also with the amplicon name, length as well as cal tm

# Enter the name of sample to run 
name="SA500B"

platedir="/home/dyap/Projects/Tumour_Evol/"$name"/QC"

amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"B_p3_amplicons.fa"

check="/home/dyap/dyap_temp/Plate_mapping"

# For Plate 1
count=1

echo "" > $check

for i in A B C D E F G H 
        do
                for j in 01 03 05 07 09 11 13 15 17 19 21 23
                do
		
		pcrplate=$i$j

		p1=$count

		        if [[ $p1 -lt "10" ]];
               		 then   {
                        	plate1="0"$p1
                        	}

			else plate1=$p1

                	fi

                primerplate1=$i$plate1

		echo $pcrplate maps to $primerplate1 on plate 1 >> $check
		
		count=`echo $count+1 | bc`
		
		if [[ $count =~ "13" ]];
		then 	{
			count=1
			}
		fi
                done

        echo $i "series done."
        done

#For Plate 2
count=1

for i in A B C D E F G H
        do
                for j in 02 04 06 08 10 12 14 16 18 20 22 24
                do

                pcrplate=$i$j

                p2=$count

                        if [[ $p2 -lt "10" ]];
                         then   {
                                plate2="0"$p2
                                }

                        else plate2=$p2

                        fi

                primerplate2=$i$plate2

                echo $pcrplate maps to $primerplate2 on plate 2 >> $check        

                count=`echo $count+1 | bc`

                if [[ $count =~ "13" ]];
                then    {
                        count=1
                        }
                fi
                done

        echo $i "series done."
        done

echo "Check 384 well qPCR plate mapping to Primer plate 1 & 2" > $check.txt
cat $check | sort >> $check.txt

exit;

