#!/bin/sh

# This script checks the master position file and then calls off the mutations
# from their VCF and makes them into another file for R processing

#######################
#Change this to run ID
run="A4GFU"
# Change this to sample name
name="SA501"

wd="/home/dyap/Projects/MiSeq_Data/Runs/"$run

##############################################
# Change this to where the positions are listed
pos="/home/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_pos.txt"
#pos="/share/lustre/backup/dyap/Projects/Tumour_Evol/SA494/Final_positions_SA494.txt"

##############################################################
# Change this to where the output should go - Project dependent
# For Single Cell Project
# outdir="/home/dyap/Projects/Single_Cell/Analysis/"$name
# For Xeno Project
outdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/"$name

# QC only
cd $wd
#ls $name-*.vcf | more
ls $name*.vcf | more
cat $pos | more


######################################################################### 
# The only unique variable is the position (format depends on input file!)
# Change match pos also

#####################################################################
# For Xeno project where chrxx:xxx-xxx
 for i in `cat $pos | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
# For single cell project
# for i in `cat $pos | awk -F"\t" '{print $3}'`
	do
	{
	echo $i
	if [[ $i > "0" ]];
	             then    	{

#####################################################################
# For Xeno project where chrxx:xxx-xxx
		matchchr=`grep $i $pos | awk -F":" '{print $1}'`
		matchpos=`grep $i $pos | awk -F"-" '{print $2}'`
# For Single Cell Project
#		matchchr=`grep $i $pos | awk -F"\t" '{print "chr"$2}'`
#		matchpos=`grep $i $pos | awk -F"\t" '{print $3}'`

	for j in `ls $name*.vcf`
#	for j in `ls $name-*.vcf`
		do

		{
		grep $matchpos $j | grep $matchchr >> $j.tmp
		}
		done
				}
	fi
	
	}
	done


	for k in `ls *.tmp`
		do
		{
		fname=`echo $k | awk -F"." '{print $1}'`
		grep "#" $fname".vcf" > $outdir"/"$fname".svcf"
		cat $k >> $outdir"/"$fname".svcf"
		cp $name*.vcf $outdir"/"
		}
		done

rm -f *.tmp

echo "Outputs are in "$outdir

exit
