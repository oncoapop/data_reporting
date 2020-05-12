#!/bin/sh

# Scipt to map qPCR plate with primer plate and also with the amplicon name, length as well as cal tm

#################################
# Enter the name of sample to run 
name="SA499"
################################

platedir="/home/dyap/Projects/Tumour_Evol/"$name"/QC"

amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_amplicons.fa"

out="/home/dyap/dyap_temp/Temp_"
chk="/home/dyap/dyap_temp/check"

# For Plate 1
count=1
##################################
# Enter position on 384 well plate
# Rows A-H of 384 well plate
top="Tumor X2 X4 NSG"

# Rows I-P of 384 well plate
bottom="Normal X3 X5"
#################################


# For Rows A-H (of 384 well plate)

for m in `echo $top`
do

pcrfile=$platedir"/"$name"-"$m
primerplate=$platedir"/"$name"_primer_plate1"

check=$chk"_"$name"-"$m".txt"
output=$out$name"-"$m".txt"

echo "PCR_well,Name_"$m",Actual_Temp,Cal_Temp,CT,Amp_Len" > $output

echo $pcrfile
echo $pcrfile > $check

echo $primerplate

for k in A B C D E F G H 
        do
                for l in 1 3 5 7 9 11 13 15 17 19 21 23
                do
		
		pcrplate=$k$l

		p1=$count

		        if [[ $p1 -lt "10" ]];
               		 then   {
                        	plate1="0"$p1
                        	}

			else plate1=$p1

                	fi

                primerplate1=$k$plate1

		echo $pcrplate maps to $primerplate1 on plate 1 >> $check

		actemp=`grep -w $pcrplate $pcrfile | awk -F" " '{print $10}'` 
		CT=`grep -w $pcrplate $pcrfile | awk -F" " '{print $8}'` 
		leftprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
		rightprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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
                        }'`

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`

#		echo $primerplate1,$primerplate,$leftprimer,$rightprimer,$ampname,$actemp,$caltemp,$ampleng,$CT

#		echo "PCR_well,Name,Actual_Temp,Cal_Temp"
#		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng
		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

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
primerplate=$platedir"/"$name"_primer_plate2"
echo $primerplate
count=1

for q in A B C D E F G H 
        do
                for r in 2 4 6 8 10 12 14 16 18 20 22 24
                do
		
		pcrplate=$q$r

		p2=$count

		        if [[ $p2 -lt "10" ]];
               		 then   {
                        	plate2="0"$p2
                        	}

			else plate2=$p2

                	fi

                primerplate2=$q$plate2

		echo $pcrplate maps to $primerplate2 on plate 2 >> $check

		actemp=`grep -w $pcrplate $pcrfile | awk -F" " '{print $10}'` 
		CT=`grep -w $pcrplate $pcrfile | awk -F" " '{print $8}'` 
		leftprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
		rightprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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
                        }'`

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`

#		echo $primerplate2,$primerplate,$leftprimer,$rightprimer,$ampname,$actemp,$caltemp,$ampleng,$CT

#		echo "PCR_well,Name,Actual_Temp,Cal_Temp"
#		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng
		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

		count=`echo $count+1 | bc`
		
		if [[ $count =~ "13" ]];
		then 	{
			count=1
			}
		fi

# For incomplete second plates
# uncomment the following conditions
#                if [[ $primerplate2 =~ "C09" ]];
#                then    {
#                        exit;
#                        }
#                fi
#
                
                done

        echo $i "series done."
        done
done
# exit;

# For the second half of the 384 well plate

# For Plate 1
count=1

 
# For rows (I-P) of 384 well plate

for m in `echo $bottom`
do

pcrfile=$platedir"/"$name"-"$m
primerplate=$platedir"/"$name"_primer_plate1"

check=$chk"_"$name"-"$m".txt"
output=$out$name"-"$m".txt"

echo "PCR_well,Name_"$m",Actual_Temp,Cal_Temp,CT,Amp_Len" > $output

#echo $pcrfile 
echo $pcrfile > $check

#echo $primerplate

for z in I J K L M N O P 
        do
                for y in 1 3 5 7 9 11 13 15 17 19 21 23
                do
		
		pcrplate=$z$y

		p1=$count

		        if [[ $p1 -lt "10" ]];
               		 then   {
                        	plate1="0"$p1
                        	}

			else plate1=$p1

                	fi

		if [[ $z =~ "I" ]];
		then 	{
			g="A"
			}
		fi
		if [[ $z =~ "J" ]];
		then 	{
			g="B"
			}
		fi
		if [[ $z =~ "K" ]];
		then 	{
			g="C"
			}
		fi
		if [[ $z =~ "L" ]];
		then 	{
			g="D"
			}
		fi
		if [[ $z =~ "M" ]];
		then 	{
			g="E"
			}
		fi
		if [[ $z =~ "N" ]];
		then 	{
			g="F"
			}
		fi
		if [[ $z =~ "O" ]];
		then 	{
			g="G"
			}
		fi
		if [[ $z =~ "P" ]];
		then 	{
			g="H"
			}
		fi

                primerplate1=$g$plate1

		echo $pcrplate maps to $primerplate1 on plate 1 >> $check

#		echo $pcrplate

		actemp=`grep -w $pcrplate $pcrfile | awk -F" " '{print $10}'` 
		CT=`grep -w $pcrplate $pcrfile | awk -F" " '{print $8}'` 
		leftprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
		rightprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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
                        }'`

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`

#		echo $primerplate1,$primerplate,$leftprimer,$rightprimer,$ampname,$actemp,$caltemp,$ampleng,$CT

#		echo "PCR_well,Name,Actual_Temp,Cal_Temp"
#		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng
		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

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
primerplate=$platedir"/"$name"_primer_plate2"
echo $primerplate
count=1

for z in I J K L M N O P
        do
                for y in 2 4 6 8 10 12 14 16 18 20 22 24
                do
		
		pcrplate=$z$y

		p2=$count

		        if [[ $p2 -lt "10" ]];
               		 then   {
                        	plate2="0"$p2
                        	}

			else plate2=$p2

                	fi

		if [[ $z =~ "I" ]];
		then 	{
			g="A"
			}
		fi
		if [[ $z =~ "J" ]];
		then 	{
			g="B"
			}
		fi
		if [[ $z =~ "K" ]];
		then 	{
			g="C"
			}
		fi
		if [[ $z =~ "L" ]];
		then 	{
			g="D"
			}
		fi
		if [[ $z =~ "M" ]];
		then 	{
			g="E"
			}
		fi
		if [[ $z =~ "N" ]];
		then 	{
			g="F"
			}
		fi
		if [[ $z =~ "O" ]];
		then 	{
			g="G"
			}
		fi
		if [[ $z =~ "P" ]];
		then 	{
			g="H"
			}
		fi

                primerplate2=$g$plate2

		echo $pcrplate maps to $primerplate2 on plate 2 >> $check

		actemp=`grep -w $pcrplate $pcrfile | awk -F" " '{print $10}'` 
		CT=`grep -w $pcrplate $pcrfile | awk -F" " '{print $8}'` 
		leftprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
		rightprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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
                        }'`

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`

#		echo $primerplate2,$primerplate,$leftprimer,$rightprimer,$ampname,$actemp,$caltemp,$ampleng,$CT

#		echo "PCR_well,Name,Actual_Temp,Cal_Temp"
#		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng
		echo $pcrplate,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

		count=`echo $count+1 | bc`
		
		if [[ $count =~ "13" ]];
		then 	{
			count=1
			}
		fi

# For incomplete second plates
# uncomment the following conditions
#                if [[ $primerplate2 =~ "C09" ]];
#                then    {
#                        exit;
#                        }
#                fi
#
                
                done

        echo $i "series done."
        done
done
exit;

