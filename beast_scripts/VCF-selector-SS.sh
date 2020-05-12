#!/bin/sh

# This script checks the master position file and then calls off the mutations
# from their VCF and makes them into another file for R processing
# This script is for single cell analysis

#######################
#Change this to run ID
run="A6MK1"
# Change this to sample name
name="VOA1306-C"
# Change this to the spikein main run 
spike="SA532-1"

########################
# For single cell analysis
# use scwd, else comment out
wd="/home/dyap/Projects/MiSeq_Data/Runs/"$run
scwd1="/home/dyap/Projects/Single_Cell/"$name
scwd="/home/dyap/Projects/Single_Cell/"$name"/"$spike"-spike-in"
mkdir $scwd1
mkdir $scwd2

cd $scwd
source=$wd"/"$name"*.*vcf"
cp $source $scwd/

##############################################
# Change this to where the positions are listed
#pos="/share/lustre/backup/dyap/Projects/Single_Cell/positions/SNV/"$name"_pos.txt"
#pos="/home/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_pos.txt"
#pos="/share/lustre/backup/dyap/Projects/Tumour_Evol/SA494/Final_positions_SA494.txt"

# Change to use AmpliconManifest files
pos=$wd"/ITH-Case2.AmpliconManifest"
##############################################################
# Change this to where the output should go - Project dependent

outdir="/home/dyap/Projects/Single_Cell/Analysis/"$name
# outdir="/home/dyap/Projects/Tumour_Evol/"$name
mkdir $outdir


# QC only
# For single cell $scwd
mkdir $scwd
cd $scwd
ls $name*.vcf | more
cat $pos | more



######################################################################### 
# The only unique variable is the position (format depends on input file!)
# Change match pos also
#####################################################################
# for i in `cat $pos | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
for i in `cat $pos | awk -F"\t" '{print $3}'`
  do
	{
	echo $i
	if [[ $i > "0" ]];
	             then    	{

#####################################################################
#		matchchr=`grep $i $pos | awk -F":" '{print $1}'`
#		matchpos=`grep $i $pos | awk -F"-" '{print $2}'`
		matchchr=`grep $i $pos | awk -F"\t" '{print "chr"$2}'`
		matchpos=`grep $i $pos | awk -F"\t" '{print $3}'`

	for j in `ls $name*.vcf`
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
mv  $outdir $outdir"-"$run

echo "Outputs are in "$outdir"-"$run

exit

