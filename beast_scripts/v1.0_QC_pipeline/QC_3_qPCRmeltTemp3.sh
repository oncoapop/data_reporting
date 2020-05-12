#!/bin/sh

# This is the one that was modified to take a generic plate qPCR
# Script to map qPCR plate with primer plate and also with the amplicon name, length as well as cal tm

#################################
# Enter the name of sample to run 
name="DAH"
################################

platedir="/home/dyap/Projects/CellMix2"

Calculated="/home/dyap/dyap_temp/Temp_cal-output"
out="/home/dyap/dyap_temp/Temp_"
chk="/home/dyap/dyap_temp/check"

# For Plate 1
count=1
##################################
# Enter position on 384 well plate
# Rows A-H of 384 well plate
top="A"

# Rows I-P of 384 well plate
bottom="B"
######################################################################
# This new pipeline now uses primer sequences sans adaptor sequences #
######################################################################

###########################################################
# For Primer Plate 1 on top (odd) of 384 well qPCR plate # 
###########################################################
# For Rows A-H (of 384 well plate)

for m in `echo $top`
do

	pcrfile=$platedir"/"$name"-"$m
	primerplate=$platedir"/"$name"_primer_plate1"
	pname=`echo $primerplate | awk -F"/" '{print $NF}'`
	
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

			# This is the well on the 384-well qPCR plate		
			well=$k$l

			p1=$count

		        if [[ $p1 -lt "10" ]];
               		 then   {
                        	plate1="0"$p1
                        	}

			else plate1=$p1

                	fi

			# This gives the mapped well on the 96-well primerplate
        	        primerplate1=$k$plate1

			# This maps the well on pcrplate to the primerplate (QC check)
			echo $primerplate1 on primer plate 1 maps to $well on the 384 well qPCR plate >> $check

			# get the CT and observed melt curve Tm from qPCR plate
			actemp=`grep -w $well $pcrfile | awk -F"\t" '{print $6}'` 
			CT=`grep -w $well $pcrfile | awk -F"\t" '{print $5}'` 

			# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
			ampname=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F"\t" '{print $1}' | awk -F_  '{print $1"_"$2"_"$3}'`
			caltemp=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F_  '{print $4}' |  sed 's/Tm=//'`
			ampleng=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F_  '{print $5}' |  sed 's/Length=//'`

			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng
			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

			# Inelegent solution to reset counter every 12th positions (ie of a 96 well plate)
			count=`echo $count+1 | bc`
		
			if [[ $count =~ "13" ]];
			then 	{
				unset ampname
				unset caltemp
				unset ampleng
				unset actemp
				unset CT
				count=1
				}
			fi
                done

	        echo $i "Row "$k" done."
        done

###########################################################
# For Primer Plate 2 on top (even) of 384 well qPCR plate # 
###########################################################

	primerplate=$platedir"/"$name"_primer_plate2"
	echo $primerplate
	pname=`echo $primerplate | awk -F"/" '{print $NF}'`

	count=1

	for q in A B C D E F G H 
        	do
                	for r in 2 4 6 8 10 12 14 16 18 20 22 24
                	do
	
			# This is the well on the 384-well qPCR plate	
			well=$q$r

			p2=$count

			        if [[ $p2 -lt "10" ]];
               			 then   {
                        		plate2="0"$p2
                        		}

				else plate2=$p2

        	        	fi
			
			# This gives the mapped well on the 96-well primerplate
	                primerplate2=$q$plate2

			# This maps the well on pcrplate to the primerplate (QC check)
			echo $primerplate2 on primer plate 2 maps to $well on the 384 well qPCR plate >> $check

			# get the CT and observed melt curve Tm from qPCR plate
			actemp=`grep -w $well $pcrfile | awk -F"\t" '{print $6}'` 
			CT=`grep -w $well $pcrfile | awk -F"\t" '{print $5}'` 

			# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
			ampname=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F"\t" '{print $1}' | awk -F_  '{print $1"_"$2"_"$3}'`
			caltemp=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F_  '{print $4}' |  sed 's/Tm=//'`
			ampleng=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F_  '{print $5}' |  sed 's/Length=//'`

			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng
			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

			count=`echo $count+1 | bc`
		
			if [[ $count =~ "13" ]];
				then 	{
					unset ampname
					unset caltemp
					unset ampleng
					unset actemp
					unset CT
					count=1
					}
			fi


		# For incomplete second plates
		# uncomment the following conditions
		#                if [[ $primerplate2 =~ "C09" ]];
		#                then    {
		#                        exit;
		#                        }
		#     		 fi
                
                done

        echo $i "series done."
        done
done

# If only one sample (script can exit here)
# exit;

#############################################################
# For Primer Plate 1 on bottom (odd) of 384 well qPCR plate # 
#############################################################
# For the second half of the 384 well plate
# Rows I-P of 384 well qPCR plate

# For Primer Plate 1
count=1

for m in `echo $bottom`
	do

	pcrfile=$platedir"/"$name"-"$m
	primerplate=$platedir"/"$name"_primer_plate1"
	pname=`echo $primerplate | awk -F"/" '{print $NF}'`

	check=$chk"_"$name"-"$m".txt"
	output=$out$name"-"$m".txt"

	echo "PCR_well,Name_"$m",Actual_Temp,Cal_Temp,CT,Amp_Len" > $output

	echo $pcrfile > $check

	for z in I J K L M N O P 
        	do
                	for y in 1 3 5 7 9 11 13 15 17 19 21 23
                	do
		
			# This is the well in 384-well qPCR plate
			well=$z$y

			p1=$count

		        if [[ $p1 -lt "10" ]];
               			 then   {
        	                	plate1="0"$p1
                        		}

					else plate1=$p1
                	fi

			# Maps the rows qPCR plate to PCR plate
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

			# This gives the mapped well on the 96-well primerplate
        	        primerplate1=$g$plate1

			# This maps the well on pcrplate to the primerplate (QC check)
			echo $primerplate1 on primer plate 1 maps to $well on the 384 well qPCR plate >> $check

			# get the CT and observed melt curve Tm from qPCR plate
			actemp=`grep -w $well $pcrfile | awk -F"\t" '{print $6}'` 
			CT=`grep -w $well $pcrfile | awk -F"\t" '{print $5}'` 

			# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
			ampname=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F"\t" '{print $1}' | awk -F_  '{print $1"_"$2"_"$3}'`
			caltemp=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F_  '{print $4}' |  sed 's/Tm=//'`
			ampleng=`grep $pname $Calculated | grep -m1 $primerplate1 | awk -F_  '{print $5}' |  sed 's/Length=//'`

			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng
			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

			# Inelegent solution to reset counter every 12th positions (ie of a 96 well plate)
			count=`echo $count+1 | bc`
		
			if [[ $count =~ "13" ]];
			then 	{
				unset ampname
				unset caltemp
				unset ampleng
				unset actemp
				unset CT
				count=1
				}
			fi
                done

	        echo $i "Row "$k" done."
        done

#############################################################
# For Primer Plate 2 on bottom (odd) of 384 well qPCR plate # 
#############################################################
	primerplate=$platedir"/"$name"_primer_plate2"
	pname=`echo $primerplate | awk -F"/" '{print $NF}'`
	echo $primerplate

count=1

for z in I J K L M N O P
        do
                for y in 2 4 6 8 10 12 14 16 18 20 22 24
                	do
		
			# This is the well position in the 384-well qPCR plate
			well=$z$y

			p2=$count

			# Maps the rows qPCR plate to PCR plate
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

			# This gives the mapped well on the 96-well primerplate
        	        primerplate2=$g$plate2

			# This maps the well on pcrplate to the primerplate (QC check)
			echo $primerplate2 on primer plate 2 maps to $well on the 384 well qPCR plate >> $check

			# get the CT and observed melt curve Tm from qPCR plate
			actemp=`grep -w $well $pcrfile | awk -F"\t" '{print $6}'` 
			CT=`grep -w $well $pcrfile | awk -F"\t" '{print $5}'` 

			# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
			ampname=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F"\t" '{print $1}' | awk -F_  '{print $1"_"$2"_"$3}'`
			caltemp=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F_  '{print $4}' |  sed 's/Tm=//'`
			ampleng=`grep $pname $Calculated | grep -m1 $primerplate2 | awk -F_  '{print $5}' |  sed 's/Length=//'`

			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng
			echo $well,$ampname,$actemp,$caltemp,$CT,$ampleng >> $output

			# Inelegent solution to reset counter every 12th positions (ie of a 96 well plate)
			count=`echo $count+1 | bc`
		
			if [[ $count =~ "13" ]];
			then 	{
				unset ampname
				unset caltemp
				unset ampleng
				unset actemp
				unset CT
				count=1
				}
			fi

		# For incomplete second plates
		# uncomment the following conditions
	      	#	 if [[ $primerplate2 =~ "C09" ]];
		#                then    {
		#     		         exit;
		#              		}
		#        fi
                done

	        echo $i "Row "$k" done."
        done

done	
exit;

