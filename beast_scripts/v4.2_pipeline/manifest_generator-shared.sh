#!/bin/sh

# This is a script to take picked SHARED positions
# grep them from the ampliconmanifest file and make a new one
# If there is a Phased SNP column then go and find that SNP by searching up and downstream of SNP
# in the original VCFs (by chr and pos) - picks all Phased SNPs by default (if present)
# and then generating a ampliconmanifest entry for that Phased SNP using the amplicon from the anchor SNP
# Automatically pulls out the Phased SNP if exists in the original VCF and PR>0.9

# This is where the output files from the initial run are stored (by project)
Project="cellmix"
outpath="/home/dyap/public_html/"$Project

# This is where the new amplicon manifest is stored as well as the input positions read from
wd="/home/dyap/Projects/PrimerDesign/"$Project"/selected"
newamp=$wd"/Shared-2.AmpliconManifest"
#file=$wd"/Shared_pos2"
file=$wd"/test_pos"

# This is the temp file where the vcfs are stored
vcfdir="/home/dyap/dyap_temp/vcf"
shared=$vcfdir"/hTertL2-HCT116-shared.vcf"
sam1="HCT116_SA036.vcf"
sam2="hTertL2_SA040.vcf"

rm -f $wd"/tmp"
rm -f $wd"/fail"
rm $newamp

	amp=`echo $outpath"/cat.AmpliconManifest.txt"`

	echo $amp
	cat $amp | head -n +6 >> $newamp

	f=$sam1 # choose one sample as the ref and the other one will be the query

	vname=$f
#	vcf=$vcfdir"/"$vname".vcf"
	vcf=$vcfdir"/"$vname

	for i in `cat $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $3}'`
		do

		# This is from the position file
		SA=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $1}'`
                Schr=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $2}'| sed 's/chr//'`
                Spos=`grep -m1 $i $file | awk -F"\t" '{print $1}' | awk -F"_" '{print $3}' | sed 's/ //'`

		if [[ $SA == "" ]] || [[ $Schr == "" ]] || [[ $Spos == "" ]]
			then 
				echo "FAILED: "$SA"_"$Schr"_"$Spos"--------------------"
				continue
		fi

		# Now we need to get it from the individual vcfs
		echo $vcf

                chr=`grep -w -m1 $i $vcf | awk -F"\t" '{print $1}'`
                pos=`grep -w -m1 $i $vcf | awk -F"\t" '{print $2}'`

		if [[ $chr == "" ]] || [[ $pos == "" ]]
			then 
				echo $SA"_"$Schr"_"$Spos
				echo $vname" : "$chr"_"$pos" - FAILED"

                       # No QC for those deepSeq'ed and found to be Het

                        # This is from the manifest file

                        SAn=`grep -m1 $i $amp | awk -F":" '{print $1}' | awk -F"_" '{print $2}'`
                        Schrn=`grep -m1 $i $amp | awk -F"\t" '{print $1}' | awk -F":" '{print $2}'| awk -F"_" '{print $1}' | sed 's/chr//'`
                        Sposn=`grep -m1 $i $amp | awk -F"\t" '{print $1}' | awk -F":" '{print $2}'| awk -F"_" '{print $2}' | sed 's/ //'`
                        QCn=`grep -m1 $i $file | awk -F"\t" '{print $2}'`

                        echo $SAn"_"$Schrn"_"$Sposn" : "$QCn

                                if [[ $QCn == "NoQC" ]] 
                                        then 
                                        echo "NO_QC: "$SAn"_"$Schrn"_"$Sposn
                                        grep -m1 $Sposn $amp | sed 's/-selected.vcf.txt//' | sed s'/Unknown/cell-line/' | sed 's/REM_/REM_noQC_/' >> $newamp
                                        echo "Position is NOT QC'ed!!!!"
                                        test="oldprimer"
                                fi
                        echo "Cont"
                        continue
		fi

		ref=`grep -m1 $i $vcf | awk -F"\t" '{print $4}'`
		alt=`grep -m1 $i $vcf | awk -F"\t" '{print $5}'`
		pr=`grep -m1 $i $vcf | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
		pgt=`grep -m1 $i $vcf | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
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

		testfile=$vcfdir"/"$fname				


		TestSchr=`grep -m1 -w $i $testfile | awk -F"\t" '{print $1}'`
		TestSpos=`grep -m1 -w $i $testfile | awk -F"\t" '{print $2}'`

		if [[ $TestSchr == "" ]] || [[ $TestSpos == "" ]]
			then 
				echo $SA"_"$Schr"_"$Spos
				echo $vname" : "$chr"_"$pos":REF_"$ref":ALT_"$alt
				echo $fname" : "$TestSchr"_"$TestSpos" - FAILED"
				continue
		fi


		TestSref=`grep -m1 $i $testfile | awk -F"\t" '{print $4}'`
		TestSalt=`grep -m1 $i $testfile | awk -F"\t" '{print $5}'`
		TestSppr=`grep -m1 $i $testfile | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
		TestSpgt=`grep -m1 $i $testfile | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`

	echo $SA"_"$Schr"_"$Spos
	echo $vname" : "$chr"_"$pos":REF_"$ref":ALT_"$alt

		if [[ $Schr != $TestSchr ]]
			then
				for loop in {1..10}
					do
					TestSchr=`grep -m$loop $i $testfile | tail -n1 | awk -F"\t" '{print $1}'`
					echo $loop
					if [[ $Schr == $TestSchr ]];
						then
						TestSref=`grep -m$loop $i $testfile | tail -n1 | awk -F"\t" '{print $4}'`
						TestSalt=`grep -m$loop $i $testfile | tail -n1 | awk -F"\t" '{print $5}'`
						TestSppr=`grep -m$loop $i $testfile | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
						TestSpgt=`grep -m$loop $i $testfile | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
						echo $fname" : "$TestSchr"_"$TestSpos":REF_"$TestSref":ALT_"$TestSalt
						break
					fi

					done
		fi

		if [[ $Schr != $chr ]]
			then
				for loop in {1..10}
					do
			                chr=`grep -m$loop $i $vcf | tail -n1 | awk -F"\t" '{print $1}'`

					if [[ $Schr == $chr ]]; 
						then
						ref=`grep -m$loop $i $vcf | tail -n1 | awk -F"\t" '{print $4}'`
						alt=`grep -m$loop $i $vcf | tail -n1 | awk -F"\t" '{print $5}'`
						pr=`grep -m$loop $i $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
						pgt=`grep -m$loop $i $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
						echo $vname" : "$chr"_"$pos":REF_"$ref":ALT_"$alt
 						break
					fi

					done
		fi

		if [[ $TestSchr == $chr ]] && [[ $TestSpos == $pos ]] && [[ $TestSref == $ref ]] && [[ $TestSalt == $alt ]] && [[ $TestSpgt == $pgt ]];
			then 
				grep -m1 $pos $amp | sed 's/-selected.vcf.txt//' | sed s'/Unknown/cell-line/' | sed 's/REM_/REM_Shared_/' >> $newamp
				echo "Position is shared!!!!"
				test="Shared"
			else
				grep -m1 $Spos $amp | sed 's/-selected.vcf.txt//' | sed s'/Unknown/cell-line/' | sed 's/REM_/REM_failedQC/' >> $newamp
				echo $vname" : "$chr"_"$pos":REF_"$ref":ALT_"$alt":PR="$pr":GT="$pgt >> $wd"/fail"
				echo $fname" : "$TestSchr"_"$TestSpos":REF_"$TestSref":ALT_"$TestSalt":PR="$TestSppr":GT="$TestSpgt >> $wd"/fail"
		fi


# This section is about grepping and annotating any Phased SNPs within the amplicon (FREE INFO)

			grep -m1 -A1 -B1 $pos $vcf >> $wd"/tmp"

			testA=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $2}'`
			testB=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $2}'`


			if [[ $testA > $start ]] && [[ $testA < $end ]];
				then	
				pchr=`grep -m1 -w -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print "chr"$1}'`
				ppos=`grep -m1 -w -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $2}'`
				pref=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $4}'`
				palt=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $5}'`
				ppr=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
				pgt=`grep -m1 -A1 $pos $vcf | tail -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
				test="Y"
			fi

			if [[ $testB > $start ]] && [[ $testB < $end ]];
				then	
				pchr=`grep -m1 -w -B1 $pos $vcf | head -n1 | awk -F"\t" '{print "chr"$1}'`
				ppos=`grep -m1 -w -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $2}'`
				pref=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $4}'`
				palt=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $5}'`
				ppr=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $1}' | sed 's/PR=//'`
				pgt=`grep -m1 -B1 $pos $vcf | head -n1 | awk -F"\t" '{print $8}' | awk -F";" '{print $9}' | sed 's/GT=//'`
				test="Y"
			fi

			if [[ $ppr > 0.9 ]] && [[ $test == "Y" ]];
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

		done


 	echo $f" sample done!"


echo $newamp "has to be post-processed..."

cat $newamp

exit;
