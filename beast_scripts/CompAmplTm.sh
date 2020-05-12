#/bin/sh

# Scipt to compare Amplicon (not primer) Tm with qPCR melting curve

cd /home/dyap/Projects/Tumour_Evol/Primer3_outputs

echo "Well,Actual_Temp,Cal_Temp" > Temp_Plate1.csv

	for i in `cat PE_amplicons_left_amplicontemp.txt | awk -F, '{print $1}'`
	do
	well=`grep $i SA494_primers.txt | awk -F"\t" '{print $1}'`
	acttemp=`grep $i SA494_primers.txt | awk -F"\t" '{print $4}'`
	caltemp=`grep $i PE_amplicons_left_amplicontemp.txt | awk -F"," '{print $2}'`
	echo $well","$acttemp","$caltemp >> Temp_Plate1.csv
	echo $well",     "$acttemp",        "$caltemp 

	done


exit;

