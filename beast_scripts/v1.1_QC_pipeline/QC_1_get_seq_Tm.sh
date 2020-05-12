#!/bin/sh

# Script to get the amplicon sequences and then calculate melting temp

projdir="/home/dyap/Projects/Takeda_T3/CG"

cd $projdir
rm -f *.tmp*

pri=$projdir"/CG_primers_ordered.txt"
# Format
# VSIG10@118506182:WSB2@118490283
# WDR77@111985275:OVGP1@111969263

p3design="/home/dyap/Projects/PrimerDesign/Splice/primer3/hct116_htert_p3_design.txt"
# Format
# ACTR8@53904009:CHDH@53875026_160-170,ACTR8@53904009:CHDH@53875026,160,170,WT,CTAAGACAAT,CATCTGACGACACCAAAAAGAAGATGTACAGCTCCATCCTAGTGGTGGGAGGTGGTTTGATGTTTCATAAAGCTCAAGAATTTCTGCAGCACAGAATTCTCAACAAAATGCCACCATCCTTCAGGCGAATTATTGAAAATGTGGATGTGATCACAAGGCCTAAGACAATCACAAAGAGTGTGTAGGCCAGCCCCGGTCACAGAGTGCACCGTATCCTGTCACTTCTGGATGTGAGGGAGAAGTGAGTCATCTCATTCCCCTCCGTGGATCAGAGGACTTGGACTAGATAGAAGCATGTGGTGTCTCCTACGAGGCCTGGGCCGGCCTGGAGCCCTGGCACGGGGAGCCCTGGGGCAGCAGCAATCCCTGGGTGCCCGGGCCCTGGCCAGCGCAGGCTCTGAGAGCCGGGACGAGTACAGCTATGTGGTGGTGGGCGCGGGCTCGGCGGGCTGCGTGCTGGCTGGGAGGCTCACGGAGGACCCCGCCGAGCGCGTGCTGCTGCTGGAGGCCGGGCCCAAGGACGTGCTC

primers="/home/dyap/Projects/PrimerDesign/Splice/primer3/hct116_htert_primer_order.txt"
# Format
#TIMM50@39980603:DLL3@39989832_325-335,chrTIMM50@39980603:DLL3@39989832,298,457,160 bp,TGTCCCGAGAGTCTCCAGATG,21,GGCTTCAGGCAGACTCTGAAG,21,N/A,N/A,TIMM50@39980603:DLL3@39989832_325-335
#TP53RK@45317771:SLC13A3@45242364_301-311,chrTP53RK@45317771:SLC13A3@45242364,195,352,158 bp,TGATCAAGCACCGCTTCCC,19,CACCAGTACACCGCCATGAG,20,N/A,N/A,TP53RK@45317771:SLC13A3@45242364_301-311
#TPD52L2@62505169:DNAJC5@62559688_222-232,chrTPD52L2@62505169:DNAJC5@62559688,129,278,150 bp,CACTGTGGAGAGCTCAAGAGG,21,TCCCCAGAGGTAGACAGTGAG,21,N/A,N/A,TPD52L2@62505169:DNAJC5@62559688_222-232

tempfile="/home/dyap/dyap_temp/Temp_cal-input"
rm -f $tempfile
# The output of the this script is /home/dyap/dyap_temp/Temp_cal-input
# Shared_chr10_78907562 CACATGGCCCCAACTTTTCCTCCTGCCTGATGTCTGGACCCACAGGGGACAGATGGTGGTCCCTGGTCACTGAAGTACCCACAGGTTTGATTATTGAGACTTTTATCTCTGTTGACTACAGGGAGGACCTGAATCAGTTTACATGGTTAAAATAAAAGCGCCATATACAAGGATAAAGCA 
# Format Name<tab>Sequence-in-upper-case

# from v1.0_QC 
# get the primers and match them with the design file
for j in `cat $primers | tail -n +2 | awk -F"," '{print $2}' | sed 's/chr//'`
                                do
					echo $j
                                        seq=`grep -m1 $j $p3design | awk -F"," '{print $7}' | tr -c 'ATCG' 'X' | tr -d "X"`
                                        left=`grep -m1 $j $primers | awk -F"," '{print $6}'`
                                        right=`grep -m1 $j $primers | awk -F"," '{print $8}'`
                                        rightrev=`grep -m1 $j $primers | awk -F"," '{print $8}'  | awk  'BEGIN {
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

                                        if [ ! -z "$seq" ]
                                                then
                                                echo $j" "$seq | grep --color=always $left | grep --colour=always $rightrev
						invseq=`echo $seq | grep -P -o '(?<='$left')[A-Z]*(?='$rightrev')'`
						echo $left,$invseq,$rightrev
                                                echo $j" "$left$invseq$rightrev >> $tempfile
                                        fi
                                done

echo "Calculating Amplicon Melting Temperatures..."

/home/dyap/Scripts/v1.1_QC_pipeline/Temp_cal.sh

exit;

# Expect Final result
# ACTR8@53904009:CHDH@53875026_Tm=81.9_Length=154

