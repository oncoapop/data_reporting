#!/bin/sh

# This is a script to take picked positions
# grep them from the ampliconmanifest file and make a new one
# If there is a Phased SNP column then go and find that SNP by searching up and downstream of SNP
# in the original VCFs (by chr and pos)
# and then generating a ampliconmanifest entry for that Phased SNP using the amplicon from the anchor SNP
# Automatically pulls out the Phased SNP if exists in the original VCF and PR>0.9

# This is where the output files from the initial run are stored (by project)
Project="cellmix"
outpath="/home/dyap/public_html/"$Project
dirs="SA036_HCT116-selected.vcf.txt SA040_HTERTL2-selected.vcf.txt"

# This is where the new amplicon manifest is stored as well as the input positions read from
wd="/home/dyap/Projects/PrimerDesign/"$Project"/selected"
newamp=$wd"/Uniq_Mixing.AmpliconManifest"
file=$wd"/Mixing_pos"

# This is the temp file where the vcfs are stored
vcfdir="/home/dyap/dyap_temp/vcf"
shared=$vcfdir"/hTertL2-HCT116-shared.vcf"
sam1="SA036_HCT116-selected.vcf.txt"
sam2="SA040_HTERTL2-selected.vcf.txt"

rm -f $wd"/tmp"
rm $newamp

for f in $sam1 $sam2
	
	do

	outdir=$outpath"/"$f
	amp=$outdir"/"$f".AmpliconManifest.txt"

	echo $amp
	cat $amp | head -n +6 >> $newamp

	vname=`echo $f | awk -F"-" '{print $1}'`
	vcf=$vcfdir"/"$vname".vcf"
	
	sp=`echo $f | awk -F"_" '{print $1}'`

	for i in `grep $sp $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $3}'`
		do
		SA=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $1}'`
                chr=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $2}' | sed 's/chr//'`
                pos=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $3}'`
#		SNP=`grep -m1 $i $file | awk -F"\t" '{print $2}'`
		start=`grep -m1 $i $amp | awk -F"\t" '{print $3}'`
		end=`grep -m1 $i $amp | awk -F"\t" '{print $4}'`

		if [[ $f == $sam1 ]]; 
			then 
				fname=`echo $sam2 | awk -F"-" '{print $1}'`
		fi

		if [[ $f == $sam2 ]]; 
			then 
				fname=`echo $sam1 | awk -F"-" '{print $1}'`
		fi

		testfile=$vcfdir"/"$fname".vcf"				

		TestSchr=`grep -m1 $i $testfile | awk -F"\t" '{print $1}'`
		TestSpos=`grep -m1 $i $testfile | awk -F"\t" '{print $2}'`

		if [[ $TestSchr == $chr ]] && [[ $TestSpos == $pos ]];
			then 
				grep -m1 $pos $amp | sed 's/-selected.vcf.txt//' | sed s'/Unknown/cell-line/' | sed 's/REM_/REM_SharedPos/' >> $newamp
				echo "==================="
				echo "Position is shared!"
				echo "Shared: "$TestSchr"_"$TestSpos
				echo "==================="
				test="FAIL"
			else
				grep -m1 $pos $amp | sed 's/-selected.vcf.txt//' | sed s'/Unknown/cell-line/' | sed 's/+40_L2-HCT/,SA040/' >> $newamp
				echo "chr"$chr"_"$pos
		fi

#		if [[ $SNP =~ "SNP" ]]
#			then

			grep -m1 -A1 -B1 $pos $vcf >> $wd"/tmp"

			testA=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $2}'`
			testB=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $2}'`


			if [[ $testA > $start ]] && [[ $testA < $end ]];
				then	
				pchr=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print "chr"$1}'`
				ppos=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $2}'`
				pref=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $4}'`
				palt=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $5}'`
				ppr=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
				pgt=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
				test="PASS"
				else
				test="FAIL"
			fi

			if [[ $testB > $start ]] && [[ $testB < $end ]];
				then	
				pchr=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print "chr"$1}'`
				ppos=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $2}'`
				pref=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $4}'`
				palt=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $5}'`
				ppr=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
				pgt=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
				test="PASS"
				else
				test="FAIL"
			fi

			if [[ $ppr > 0.9 ]] && [[ $test == "PASS" ]];
				then
				pSA=`grep $pos $newamp | awk -F":" '{print $1}' | awk -F" " '{print $1}'`
				pREF="REF_"$pref
				pALT="ALT_"$palt
				pST=`grep $pos $newamp | awk -F":" '{print $5}' | awk -F" " '{print $1}'`
				pAN=`grep $pos $newamp | awk -F":" '{print $6}' | awk -F" " '{print $1}'`
				
				if [[ $pgt == "0/1" ]]; then pRE="REM_phased_HetSNP"
				fi
				if [[ $pgt == "1/1" ]]; then pRE="REM_phased_HomSNP"
				fi
				if [[ $pgt == "" ]]; then pRE="REM_phased_unknownSNP"
				fi

				pAc=`grep $pos $newamp | awk -F"\t" '{print $2}' | awk -F" " '{print $1}'`
				pAst=`grep $pos $newamp | awk -F"\t" '{print $3}' | awk -F" " '{print $1}'`
				pAen=`grep $pos $newamp | awk -F"\t" '{print $4}' | awk -F" " '{print $1}'`
				pAle=`grep $pos $newamp | awk -F"\t" '{print $5}' | awk -F" " '{print $1}'`
				pAri=`grep $pos $newamp | awk -F"\t" '{print $6}' | awk -F" " '{print $1}'`

				echo -e $pSA":"$pchr"_"$ppos":"$pREF":"$pALT":"$pST":"$pAN":"$pRE"\t"$pAc"\t"$pAst"\t"$pAen"\t"$pAle"\t"$pAri >> $newamp
				echo "Phased SNP: "$pchr"_"$ppos
				test="N"
			fi 
#		fi

		done

 	echo $f" sample done!"

	done

echo $newamp "has to be post-processed..."

cat $newamp

exit;
