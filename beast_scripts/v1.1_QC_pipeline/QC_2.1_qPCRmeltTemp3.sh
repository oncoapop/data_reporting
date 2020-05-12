#!/bin/sh

# This is a script to format the melt temp into the same file

wd="/home/dyap/Projects/Takeda_T3/CG_factors/siPanel"

#infile=$wd"/siPanel_Plate2_summary.csv"
#pcrfile=$wd"/siPanel_Plate2_summary_unix.csv"
#infile=$wd"/siPanel_Plate3_summary.csv"
#pcrfile=$wd"/siPanel_Plate3_summary_unix.csv"
infile=$wd"/siPanel_Plate4_summary.csv"
pcrfile=$wd"/siPanel_Plate4_summary_unix.csv"

cat $infile | tr '\015' '\012' > $pcrfile

#	Sample.Name	Target.Name	CT	Ct.Mean	GAPDH.Ct	Col6	Col7	Automatic.Ct.Threshold	Ct.Threshold	Baseline.Start	Baseline.End	Tm1	Tm2	Tm3	MTP
#1	siNT-Plate 2	ACTR8@53904009:CHDH@53875026	36.96	36.96	17.1859417	NA	NA	FALSE	0.247344	3	32	83.01065063	NA	NA	N

Calculated="/home/dyap/dyap_temp/Temp_cal-output"
#ZFAND1@82615204:SLC10A5@82607158_Tm=77.9_Length=152
#ZNF215@6949038:NLRP14@7059797_Tm=77.9_Length=153

pname=`echo $pcrfile | awk -F"/" '{print $NF}' | sed 's/_sum*.*$//'`
echo $pname

for j in `cat $pcrfile | awk -F"," '{print $2}' | tail -n +2 | sort -u | tr " " "_"`
	do
	file=`echo $pname"-"$j"_Temp.csv"`
	output=$wd"/"$file
	echo "#############  "$output
        	echo "Sample.Name,Target.Name,Ct.Mean,Tm1,Tm2,Tm3,caltemp,ampleng" > $output

	for i in `cat $pcrfile | awk -F"," '{print $3}' | tail -n +2 | sort -u`
		do
		echo $i"-"$j
        	# get the RQ and observed melt curve Tm from qPCR plate
        	actemp1=`grep $i $pcrfile | grep $j | awk -F"," '{print $10}'`
        	actemp2=`grep "$i" $pcrfile | grep "$j" | awk -F"," '{print $11}'`
        	actemp3=`grep "$i" $pcrfile | grep "$j" | awk -F"," '{print $12}'`
        	sample=`grep "$i" $pcrfile | grep "$j" | awk -F"," '{print $2}'`
        	CT=`cat $pcrfile | grep $j | grep $i | awk -F"," '{print $4}' | tr -d ","`

		# Get the name and calculated Tm and length from $Calculated (Cal_Temp output)
		caltemp=`grep "$i" $Calculated | awk -F_  '{print $2}' |  sed 's/Tm=//'`
		ampleng=`grep "$i" $Calculated | awk -F_  '{print $3}' |  sed 's/Length=//'`

        	echo $sample,$i,$CT,$actemp1,$actemp2,$actemp3,$caltemp,$ampleng
        	echo $sample,$i,$CT,$actemp1,$actemp2,$actemp3,$caltemp,$ampleng >> $output
		done
	done

# Command line R script
###############################################################
# Make sure QC_3.1 is in automated mode to accept command args
##############################################################

echo "Rscript QC_3.1_Temp_Corr.R --no-save --no-restore --args "$pname
Rscript QC_3.1_Temp_Corr.R --no-save --no-restore --args $pname

exit;
