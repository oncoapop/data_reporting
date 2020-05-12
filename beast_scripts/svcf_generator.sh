#!/bin/sh

# This script checks the master position file and then calls off the mutations
# from their VCF and makes them into another file for R processing

#######################
#Change this to run ID
#run="AGAJY"
run="AGALJ"
# Change this to sample name
name="SA535"
# Change this to the alignment run number
no="2"

wd="/home/dyap/Projects/MiSeq_Data/Runs/"$run

if [ ! -d "$wd" ]; then
 	mkdir $wd
fi

##############################################
# Change this to where the AmpliconManifest file is located
pos="/home/dyap/temp/"$name"_Amplicon_Manifest"

# Change this to where the alignments are
indir="/home/dyap/temp/"$run"/Alignment"$no

##############################################################
# Change this to where the output should go - Project dependent
# For Xeno Project
outdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/"$name
if [ ! -d "$outdir" ]; then
 	mkdir $outdir
fi

# Copy to common location
cd $indir
ls
cp *.bam $wd
cp *.bam.bai $wd
cp *.vcf $wd
cp *.vcf.gz $wd
cp SampleSheetUsed.csv $wd

echo "Files copied..."
######################################################################### 
# The only unique variable is the position (format depends on input file!)
# Change match pos also

#####################################################################
# For Xeno project where chrxx_xxx
 for i in `grep "chr" $pos | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
	do
	{
	echo $i
	if [[ $i > "0" ]];
	             then    	{

#####################################################################
# For Xeno project where chrxx_xxx
		matchchr=`grep $i $pos | awk -F":" '{print $1}'`
		matchpos=`grep $i $pos | awk -F"_" '{print $2}'`
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
