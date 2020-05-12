#!/bin/sh

# This is a pipeline script to Generate the Manifest file given primers and positions
# It generates the exact length of the amplicon and primers for the manifest
# It also does a QC on the primers and positions (from 2 plates of PCR primers (1 plate=forward and reverse pairs))

# Enter the name of sample to run 
name="SA494"

clear 

echo "Creating the Manifest file for :" $name
echo "Press enter to confirm or ctrl-c to exit"
read ans

mkdir /home/dyap/Projects/Tumour_Evol/$name

platedir="/home/dyap/Projects/Tumour_Evol/"$name"/QC"

mkdir $platedir

amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_amplicons.fa"
order="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_order.txt"

# Working Directory
dir="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/"

fname=$name"_p3_output.txt"

# Name processing - Do not change

# Source and Output directories where working files are stored
sourcedir=$dir"SNV/"
outdir="/home/dyap/dyap_temp/"

tmpfile=`echo $fname | awk -F_ '{print $1}'`
# contains the list of positions
infile=$sourcedir$tmpfile"_pos.txt"
# contains the output of this script
outfile=$outdir$tmpfile"_Manifest"
wtpos=$dir$tmpfile"_WT_positions.csv"

manfile=$platedir"/A_"$name".AmpliconManifest"

readfile=$name"_anno.txt"

# This gets the WT position in the middle of a 11bp sequences for matching SNV
echo Processing files
cat $wtpos | awk -F, '{print $2"_"$3,$4}' | tr -d '"' | tr "c" "C" > $outfile"66.tmp"

# Clears the tmp Files
# Headers for all files except $outfile.txt since the header has to be added after sorting for unique lines

echo "" > $outfile.txt
echo ">>>>>>left primer>>>>>...........................>>>>(5bp)(SNV)(5bp)>>>>>........................................................>>>>>right primer>>>>>>>" > $outfile"_check.tmp"
echo "leftpos,leftprimdist,<--------leftdist--------->-----SNV position-----<------------------------rightdist------------------------>,rightprimdist,rightpos" >> $outfile"_check.tmp"
echo "  ">> $outfile"_check.tmp"

echo "SNV_ID,AmpLength,LeftPos,LeftSeqLimit,SNV-coordinates,RightSeqLimit,RightPos" > $outfile"hiseq.tmp"


# No of plates
for no in 1 2

do
 
primerplate=$platedir"/"$name"_primer_plate"$no


for k in A B C D E F G H 
        do
                for l in 01 02 03 04 05 06 07 08 09 10 11 12

                do
		
                primerplate1=$k$l
echo $primerplate1

# This makes sure that both primers that we ordered are found on the same amplicon
		leftprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
# echo $leftprimer
# THe right primer has to be reverse complemented to be able to match the amplicon 
		rightorder=`grep $primerplate1 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'`
		rightprimerrc=`grep $primerplate1 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//' `
		rightprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
                                j = n = split("A C G T", t)
                                for (i = 0; ++i <= n;)
                                map[t[i]] = t[j--]
                                        }
                        {
                                if (/LEFT/) print
                        else {
                                for (i = length; i; i--)
                                printf "%s", map[substr($0, i, 1)]
                                print x
                                }
                        }'`
# echo $rightprimer
		ordname=`grep -B1 -A1 $leftprimer $order | grep -B2 $rightorder | grep "Chr"`
		ampname=`grep -B1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' |  sed 's/.*Chr/Chr/'`
				for m in $ampname
					do
echo $ampname
echo $ordname
				pos=`grep -A1 $m $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $2}' |  sed 's/.*Chr/Chr/'`
echo $pos
				chrom=`grep -A1 $m $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1}'  | sed 's/.*Chr/chr/'`
echo $chrom
#				caltemp=`grep -B1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//'  | sed 's/.*Chr/Chr/'`
				ampleng=`grep -B1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//'  | sed 's/.*Chr/Chr/'`

                        	match=`grep $m $outfile"66.tmp" | awk -F" " '{print $2}'`

# This special case gets the intervening sequence between the left primer and SNP but not the sequences themselves!
echo $match

# For checking this is helpful to check that all is WELL
	                        grep -m1 -B1 --color=always $leftprimer $amplicon | grep --color=always $rightprimer | grep --color=always $match >> $outfile"_check.tmp"


        	                leftdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$leftprimer')[A-Z]*(?='$match')' | wc -c`
              	        	rightdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$match')[A-Z]*(?='$rightprimer')' | wc -c`
				rightprimdist=`echo $rightprimer | wc -c`
				leftprimdist=`echo $leftprimer | wc -c`
				leftpos=`echo $pos - 5 - $leftdist - $leftprimdist | bc `
				rightpos=`echo $pos + 5 + $rightdist + $rightprimdist | bc `
echo ___________________________
# Position checking!
				echo $leftpos,$leftprimdist,$leftdist,---$pos---,$rightdist,$rightprimdist,$rightpos >> $outfile"_check.tmp"

				readfile=$name"_anno.txt"		

				label=`grep $m $sourcedir$readfile | tr " " "_" | tr "," "-"`

					if [ -z "$label" ] 
						then
							label=$m
					fi

# This is for the generation of the manifest file
				echo $m,$chrom,$leftpos,$rightpos,$leftprimdist,$rightprimdist,$label >> $outfile.txt

# This is info for HiSeq 100PE
				leftseqlim=`echo $leftpos + 100 | bc `
				rightseqlim=`echo $rightpos - 100 | bc `

				echo $m,$ampleng,$leftpos,$leftseqlim,$pos,$rightseqlim,$rightpos >> $outfile"hiseq.tmp"
# This is info for Supplemntal Figs
				ucscleftpos=`echo $leftpos +2 | bc `
				ucscrightpos=`echo $rightpos -2 | bc `
				echo $name,$m, $ucscleftpos,$ucscrightpos, $leftprimer,$rightprimerrc >> $outfile"_isPCR_chk.tmp"

# For incomplete  plates
# uncomment the following conditions
#                			if [[ $primerplate1 =~ "C09" ]];
#                				then    {
#                        					exit;
#                        				}
#                			fi
#
	  				done              
                done

        done

echo "Plate "$no" done"

done


echo "[Header]" > $outfile.tmp					
echo $name"  Manifest Version,1" >> $outfile.tmp				
echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $outfile.tmp
echo "	" >> $outfile.tmp
echo "[Regions]" >> $outfile.tmp					
echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" >> $outfile.tmp

cat $outfile.txt | sort -u > $outfile99.tmp

cat $outfile.tmp | tr "," "\t" > $manfile
cat $outfile99.tmp | tr "," "\t" >> $manfile


cat $outfile"_check.tmp" > $platedir"/"$name"_check.txt"
cat $outfile"hiseq.tmp" | sort -u -r  > $platedir"/"$name"_Hiseq.csv"
cat $outfile"_isPCR_chk.tmp" | sort -u > $platedir"/"$name"_isPCR_chk.csv"

rm $outfile*

exit
