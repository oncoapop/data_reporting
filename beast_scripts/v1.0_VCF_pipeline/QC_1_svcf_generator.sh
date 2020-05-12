#!/bin/sh

# This script checks the master position file and then calls off the mutations
# from their VCF and makes them into another file for R processing

#######################
#Change this to run ID
run="AUBU6"
#run="AGAJY"
#run="AGALJ"
#run="AJH03"

# Change this to sample name
name="SA532"
# Change this to the alignment run number
#no="2"
#no=5

wd="/home/dyap/Projects/Tumour_Evol/Xenodrug/"$name"/"$run
if [ ! -d "$wd" ]; then
 	mkdir $wd
fi

############################################################
# Usually this will be the MiSeq run or if only Fastq is generated
# Then this will be the import directory from MiSeq Server

# Change this to where the AmpliconManifest file is located
pos="/home/dyap/Projects/Tumour_Evol/Xenodrug/"$name"/xenodrug."$name"-filtered.AmpliconManifest.txt"
#pos="/home/dyap/temp/"$name"_Amplicon_Manifest"
#cat $pos

# Change this to where the alignments are
#indir="/home/dyap/temp/"$run"/Alignment"$no
indir="/home/dyap/Projects/Tumour_Evol/Xenodrug/"$run

##############################################################
# Change this to where the output should go - Project dependent
# For Xeno Project
outdir=$wd

if [ ! -d "$outdir" ]; then
 	mkdir $outdir
fi

# Copy to common location

cd $indir
ls
#cp *.bam $wd
#cp *.bam.bai $wd
#cp *.vcf $wd
cp SA532*.vcf 131*.vcf $wd
#cp 127*.vcf SA535*.vcf $wd
#cp *.vcf.gz $wd
cp SampleSheetUsed.csv $wd

#echo "Files copied..."
echo "Only Selected Files copied..."
######################################################################### 
# The only unique variable is the position (format depends on input file!)
# Change match pos also

#####################################################################
# For Xeno project where SP_SA535:chrX_73806408:REF_A:ALT_C:ST_Unknown:AN_$ann_gene:REM_
# for i in `grep "chr" $pos | awk -F"\t" '{print $1}' | awk -F":" '{print $2}'`

# For DrugXeno project where SAnnn_chrxx_xxx
 for i in `grep "chr" $pos | awk -F"\t" '{print $1}' | awk -F"_" '{print $2"_"$3}'`
	do
	{
	echo $i
	if [[ $i > "0" ]];
	             then    	{

#####################################################################
# For get from $i where chrxx_xxxx

		matchchr=`echo $i | awk -F"_" '{print $1}'`
		matchpos=`echo $i | awk -F"_" '{print $2}'`
		

	for j in `ls 131*.vcf SA532*.vcf`
		do

		{
		echo $j
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
		cat $k | sort -u >> $outdir"/"$fname".svcf"
		cp $name*.vcf $outdir"/"
		}
		done

# This is essential to remove the *tmp or else second iteration will
# have duplicates
rm -f *.tmp

echo "Outputs are in "$outdir

exit
