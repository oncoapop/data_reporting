#!/bin/sh

# Scipt to generate the Manifest file for subset of positions for
# Xenographs samples or any single cell project given primers
# This uses the file that is used to generate the PCR primers

# Enter the name of sample to run
name="SA029"

platedir="/home/dyap/Projects/Single_Cell/positions/SNV"

# 
amplicon="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions/SNV/"$name$
amplicon="/home/dyap/Projects/Single_Cell/positions/"$name"_positions.csv"

# Working Directory
dir="/home/dyap/Projects/Single_Cell/positions/"

out="/home/dyap/dyap_temp/Temp_"
chk="/home/dyap/dyap_temp/check"

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

readfile=$name"_anno.txt"

# This gets the WT position in the middle of a 11bp sequences for matching SNV
echo Processing files
cat $wtpos | awk -F, '{print $2"_"$3,$4}' | tr -d '"' > $outfile"66.tmp"


# There should be 48 positions per single cell run (uses 48x48 FLD 
# system)

primerplate=$platedir"/"$name"_primer_singlecells"

m="SS"
check=$chk"_"$name"-"$m".txt"
output=$out$name"-"$m".txt"


for k in `cat $primerplate | awk '{print $4}'  | sed 's/ACACTGACGACATGGTTCTACA//'`

do
echo "This is k, left primer="$k
       
# This makes sure that both primers that we ordered are found on the same amplicon
		leftprimer=$k
# The right primer has to be reverse complemented to be able to match 
# the amplicon 
		rightprimer=`grep -m1 $k $primerplate | awk -F" " '{print $5}' | sed 's/TACGGTAGCAGAGACTTGGTCT//'  | awk  'BEGIN {
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

echo "Right primer="$rightprimer
echo "This is amplicon="$amplicon
		ampname=`grep $leftprimer $amplicon | awk -F,  '{print $2"_"$3}'| tr -d '"' `
echo "AMPNAME!!!" $ampname
		pos=`grep $leftprimer $amplicon | awk -F,  '{print $3}'| tr -d '"' `
		chrom=`grep $leftprimer $amplicon | awk -F,  '{print $2}'| tr -d '"' `
echo THis is it $ampname,$chrom,$pos

                        match=`grep $ampname $outfile"66.tmp" | awk -F" " '{print $2}'`
echo MATCH $match
# This special case gets the intervening sequence between the left primer and SNP but not the sequences themselves!

                        leftsnv=`grep $leftprimer $amplicon | grep -P -o '(?<='$leftprimer')[A-Z]*(?='$match')' | wc -c`
                        ivsdist=`grep $leftprimer $amplicon | grep -P -o '(?<='$leftprimer')[A-Z]*(?='$rightprimer')' | wc -c`
echo "IVS="$leftsnv,$ivsdist
			rightprimdist=`echo $rightprimer | wc -c`
			leftprimdist=`echo $leftprimer | wc -c`
			leftpos=`echo $pos - 5 - $ivsdist - $leftprimdist | bc `
			rightpos=`echo $leftpos + $leftprimerdist + $ivsdist + $rightprimerdist | bc `

echo left snv,$leftsnv,ivsdist,$ivsdist,rightprimdist,$rightprimdist,leftprimdist,$leftprimdist,leftpos,$leftpos,			
readfile=$name"_anno.txt"		
			label=`grep $ampname $sourcedir$readfile | tr " " "_" | tr "," "-"`

		if [ -z "$label" ] 
			then
				label=$ampname
		fi

#		echo $primerplate1,$primerplate,$leftprimer,$rightprimer,$ampname,$actemp,$caltemp,$ampleng,$CT

# This is for the generation of the manifest file
	echo $label,$chrom,$leftpos,$rightpos,$leftprimdist,$rightprimdist >> $outfile.txt

        done


cat $outfile.txt | sort -u > $platedir"/"$name".Manifest.txt"

echo "done."

rm $outfile*

exit
