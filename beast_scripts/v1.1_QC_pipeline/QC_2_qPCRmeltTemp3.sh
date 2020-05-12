#!/bin/sh

# This is a script to format the melt temp into the same file

wd="/home/dyap/Projects/Takeda_T3/CG"

pcrfile=$wd"/siPanel_Plate2"
#Well Position   Sample Name     Target Name     RQ      CT      Ct Mean Ct SD   Delta Ct Mean   Delta Ct SE     Delta Delta Ct  Ct Threshold    Baseline Start  Baseline End    Tm1     Tm2     Tm3     MTP
#A1      siNT-Plate 2    ACTR8@53904009:CHDH@53875026    1.000   35.344  35.344          17.051          0.000   0.040   3       32      83.011                  N

Calculated="/home/dyap/dyap_temp/Temp_cal-output"
#ZFAND1@82615204:SLC10A5@82607158_Tm=77.9_Length=152
#ZNF215@6949038:NLRP14@7059797_Tm=77.9_Length=153

pname=`echo $pcrfile | awk -F"/" '{print $NF}'`

for j in `cat $pcrfile | awk -F"\t" '{print $2}' | tail -n +2 | sort -u | tr " " "_"`
	do
	file=`echo $pname"-"$j"_Temp.csv"`
	output=$wd"/"$file
	echo "#############  "$output
        	echo "Sample,Primer,RQ,actemp,actemp2,actemp3,caltemp,ampleng" > $output

	for i in `cat $pcrfile | awk -F"\t" '{print $3}' | tail -n +2 | sort -u`
		do
		echo $i"-"$j
        	# get the RQ and observed melt curve Tm from qPCR plate
        	actemp1=`grep $i $pcrfile | grep $j | awk -F"\t" '{print $14}'`
        	actemp2=`grep "$i" $pcrfile | grep "$j" | awk -F"\t" '{print $15}'`
        	actemp3=`grep "$i" $pcrfile | grep "$j" | awk -F"\t" '{print $16}'`
        	sample=`grep "$i" $pcrfile | grep "$j" | awk -F"\t" '{print $2}'`
        	RQ=`cat $pcrfile | grep $j | grep $i | awk -F"\t" '{print $4}' | tr -d ","`

		# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
		caltemp=`grep "$i" $Calculated | awk -F_  '{print $2}' |  sed 's/Tm=//'`
		ampleng=`grep "$i" $Calculated | awk -F_  '{print $3}' |  sed 's/Length=//'`

        	echo $sample,$i,$RQ,$actemp1,$actemp2,$actemp3,$caltemp,$ampleng
        	echo $sample,$i,$RQ,$actemp1,$actemp2,$actemp3,$caltemp,$ampleng >> $output
		done
	done

# Command line R script
Rscript QC_3_Temp_Corr.R --no-save --no-restore --args $pname

exit;
