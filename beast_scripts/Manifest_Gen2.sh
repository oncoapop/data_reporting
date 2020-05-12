#!/bin/sh

# This is a pipeline script to Generate the Manifest file given primers and positions
# It generates the exact length of the amplicon and primers for the manifest
# It also does a QC on the primers and positions (from 2 plates of PCR primers (1 plate=forward and reverse pairs))

# Enter the name of sample to run 
name="SA501"

platedir="/home/dyap/Projects/Tumour_Evol/"$name"/QC"

amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name"_p3_amplicons.fa"

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

# For Plate 1
count=1

 
primerplate=$platedir"/"$name"_primer_plate1"


for k in A B C D E F G H 
        do
                for l in 1 3 5 7 9 11 13 15 17 19 21 23
                do
		
		p1=$count

		        if [[ $p1 -lt "10" ]];
               		 then   {
                        	plate1="0"$p1
                        	}

			else plate1=$p1

                	fi

                primerplate1=$k$plate1
echo $primerplate1

# This makes sure that both primers that we ordered are found on the same amplicon
		leftprimer=`grep $primerplate1 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
# THe right primer has to be reverse complemented to be able to match the amplicon 
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

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		pos=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $2}' | tr -d ">"`
		chrom=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1}' | tr -d ">" | tr "C" "c"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`

                        match=`grep $ampname $outfile"66.tmp" | awk -F" " '{print $2}'`
# This special case gets the intervening sequence between the left primer and SNP but not the sequences themselves!

# For checking this is helpful to check that all is WELL
                        grep -m1 -B1 --color=always $leftprimer $amplicon | grep --color=always $rightprimer | grep --color=always $match >> $outfile"_check.tmp"


                        leftdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$leftprimer')[A-Z]*(?='$match')' | wc -c`
                        rightdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$match')[A-Z]*(?='$rightprimer')' | wc -c`
			rightprimdist=`echo $rightprimer | wc -c`
			leftprimdist=`echo $leftprimer | wc -c`
			leftpos=`echo $pos - 5 - $leftdist - $leftprimdist | bc `
			rightpos=`echo $pos + 5 + $rightdist + $rightprimdist | bc `

# Position checking!
echo $leftpos,$leftprimdist,$leftdist,---$pos---,$rightdist,$rightprimdist,$rightpos >> $outfile"_check.tmp"


readfile=$name"_anno.txt"		
			label=`grep $ampname $sourcedir$readfile | tr " " "_" | tr "," "-"`

		if [ -z "$label" ] 
			then
				label=$ampname
		fi

# This is for the generation of the manifest file
	echo $ampname,$chrom,$leftpos,$rightpos,$leftprimdist,$rightprimdist,$label >> $outfile.txt

# This is info for HiSeq 100PE
			leftseqlim=`echo $leftpos + 100 | bc `
			rightseqlim=`echo $rightpos - 100 | bc `

	echo $ampname,$ampleng,$leftpos,$leftseqlim,$pos,$rightseqlim,$rightpos >> $outfile"hiseq.tmp"

		count=`echo $count+1 | bc`
		
		if [[ $count =~ "13" ]];
		then 	{
			count=1
			}
		fi

# For incomplete  plates
# uncomment the following conditions
#                if [[ $primerplate1 =~ "C09" ]];
#                then    {
#                        exit;
#                        }
#                fi
#
                
                done

        done

echo "Plate 1 done"

#For Plate 2
primerplate=$platedir"/"$name"_primer_plate2"
count=1

for q in A B C D E F G H 
        do
                for r in 2 4 6 8 10 12 14 16 18 20 22 24
                do
		
		p2=$count

		        if [[ $p2 -lt "10" ]];
               		 then   {
                        	plate2="0"$p2
                        	}

			else plate2=$p2

                	fi

                primerplate2=$q$plate2

echo $primerplate2

# This makes sure that both primers that we ordered are found on the same amplicon
		leftprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $3}' | sed 's/ACACTGACGACATGGTTCTACA//'`
# THe right primer has to be reverse complemented to be able to match the amplicon 
		rightprimer=`grep $primerplate2 $primerplate | awk -F" " '{print $4}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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

		ampname=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1"_"$2}' | tr -d ">"`
		pos=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $2}' | tr -d ">"`
		chrom=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $1}' | tr -d ">" | tr "C" "c"`
		caltemp=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $4}' |  sed 's/Tm=//' | tr -d ">"`
		ampleng=`grep -B1 -m1 $leftprimer $amplicon | grep -A1 "_0_" | grep -B1 $rightprimer | grep ">" | awk -F_  '{print $5}' |  sed 's/Length=//' | tr -d ">"`


                        match=`grep $ampname $outfile"66.tmp" | awk -F" " '{print $2}'`

# For checking this is helpful to check that all is WELL
                        grep -m1 -B1 --color=always $leftprimer $amplicon | grep --color=always $rightprimer | grep --color=always $match >> $outfile"_check.tmp"


                        leftdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$leftprimer')[A-Z]*(?='$match')' | wc -c`
                        rightdist=`grep -m1 -B1 $leftprimer $amplicon | grep -P -o '(?<='$match')[A-Z]*(?='$rightprimer')' | wc -c`
			rightprimdist=`echo $rightprimer | wc -c`
			leftprimdist=`echo $leftprimer | wc -c`
			leftpos=`echo $pos - 5 - $leftdist - $leftprimdist | bc `
			rightpos=`echo $pos + 5 + $rightdist + $rightprimdist | bc `

# Position checking!
echo $leftpos,$leftprimdist,$leftdist,---$pos---,$rightdist,$rightprimdist,$rightpos >> $outfile"_check.tmp"


readfile=$name"_anno.txt"		
			label=`grep $ampname $sourcedir$readfile | tr " " "_" | tr "," "-"`

		if [ -z "$label" ]
		then 	
			label=$ampname
		fi

# This is for the generation of the manifest file
	echo $ampname,$chrom,$leftpos,$rightpos,$leftprimdist,$rightprimdist,$label >> $outfile.txt

# This is info for HiSeq 100PE
			leftseqlim=`echo $leftpos + 100 | bc `
			rightseqlim=`echo $rightpos - 100 | bc `

	echo $ampname,$ampleng,$leftpos,$leftseqlim,$pos,$rightseqlim,$rightpos >> $outfile"hiseq.tmp"

		count=`echo $count+1 | bc`
		
		if [[ $count =~ "13" ]];
		then 	{
			count=1
			}
		fi

# For incomplete second plates
# uncomment the following conditions
#                if [[ $primerplate2 =~ "C09" ]];
#                then    {
#                        exit;
#                        }
#                fi
#
                
                done

        done


cat $outfile.txt | sort -u > $outfile99.tmp


echo "[Header]" > $outfile.tmp					
echo $name"  Manifest Version,1" >> $outfile.tmp				
echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" >> $outfile.tmp
echo "	" >> $outfile.tmp
echo "[Regions]" >> $outfile.tmp					
echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Label" >> $outfile.tmp

cat $outfile.tmp | tr "," "\t" > $manfile
cat $outfile99.tmp | tr "," "\t" >> $manfile

cat $outfile"_check.tmp" > $platedir"/"$name"_check.txt"
cat $outfile"hiseq.tmp" > $platedir"/"$name"_Hiseq.csv"

rm $outfile*

exit
