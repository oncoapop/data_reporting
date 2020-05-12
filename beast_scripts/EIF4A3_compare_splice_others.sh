#!/bin/sh
# Script to compare the splice ratios that we obtain with the ratios obtained by Mao et al
# note that Mao et is a mouse E10.5 dataset 

iswd="/share/lustre/abashash/EIF4A3/miso_pipeline/OUTPUT/RUN/TK_18_19_20_miso/outputs/results"


ourwd="/home/dyap/Projects/EIF4A3_paper/toTransfer/EIF4A3project/from_the_other_paper/ourSites"

cd $ourwd

HeLa_SE_T202=`cat all_events | grep -w SE | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_SE_T595=`cat all_events | grep -w SE | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_RI_T202=`cat all_events | grep -w RI | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_RI_T595=`cat all_events | grep -w RI | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_A5_T202=`cat all_events | grep -w A5SS | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_A5_T595=`cat all_events | grep -w A5SS | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_A3_T202=`cat all_events | grep -w A3SS | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_A3_T595=`cat all_events | grep -w A3SS | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_MXE_T202=`cat all_events | grep -w MXE | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_MXE_T595=`cat all_events | grep -w MXE | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_AFE_T202=`cat all_events | grep -w AFE | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_AFE_T595=`cat all_events | grep -w AFE | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

HeLa_ALE_T202=`cat all_events | grep -w ALE | grep HeLa_T_2 | sort -uk 3,3 | grep T_202 | wc -l`
HeLa_ALE_T595=`cat all_events | grep -w ALE | grep HeLa_T_5 | sort -uk 3,3 | grep T_595 | wc -l`

echo "HeLa dataset"
echo "SE :"
echo "T-202 : "$HeLa_SE_T202
echo "T-595 : "$HeLa_SE_T595
comb=`echo "$HeLa_SE_T202 + $HeLa_SE_T595" | bc`
echo "Combined : "$comb

echo "RI :"
echo "T-202 : "$HeLa_RI_T202
echo "T-595 : "$HeLa_RI_T595
comb=`echo "$HeLa_RI_T202 + $HeLa_RI_T595" | bc`
echo "Combined : "$comb

echo "A5SS :"
echo "T-202 : "$HeLa_A5_T202
echo "T-595 : "$HeLa_A5_T595
comb=`echo "$HeLa_A5_T202 + $HeLa_A5_T595" | bc`
echo "Combined : "$comb

echo "A3SS :"
echo "T-202 : "$HeLa_A3_T202
echo "T-595 : "$HeLa_A3_T595
comb=`echo "$HeLa_A3_T202 + $HeLa_A3_T595" | bc`
echo "Combined : "$comb

echo "MXE :"
echo "T-202 : "$HeLa_MXE_T202
echo "T-595 : "$HeLa_MXE_T595
comb=`echo "$HeLa_MXE_T202 + $HeLa_MXE_T595" | bc`
echo "Combined : "$comb

echo "AFE :"
echo "T-202 : "$HeLa_AFE_T202
echo "T-595 : "$HeLa_AFE_T595
comb=`echo "$HeLa_AFE_T202 + $HeLa_AFE_T595" | bc`
echo "Combined : "$comb

echo "ALE :"
echo "T-202 : "$HeLa_ALE_T202
echo "T-595 : "$HeLa_ALE_T595
comb=`echo "$HeLa_ALE_T202 + $HeLa_ALE_T595" | bc`
echo "Combined : "$comb



otherwd="/home/dyap/Projects/EIF4A3_paper/comparisons"

# Using their own criteria for selection of relevant events by MISO

cd $otherwd

Mao_SE=`cat  Mao_Eif4a3_SE.txt | awk -F, '$11 > 20 {print $1,$2,$11}' | wc -l`
Mao_RI=`cat  Mao_Eif4a3_RI.txt | awk -F, '$11 > 20 {print $1,$2,$11}' | wc -l`
Mao_A5=`cat  Mao_Eif4a3_A5.txt | awk -F, '$11 > 20 {print $1,$2,$11}' | wc -l`
Mao_A3=`cat  Mao_Eif4a3_A3.txt | awk -F, '$11 > 20 {print $1,$2,$11}' | wc -l`

echo "Mao et al (mouse embryos):"
echo "SE : "
echo $Mao_SE
echo "RI : "
echo $Mao_RI
echo "A5 : "
echo $Mao_A5
echo "A3 : "
echo $Mao_A3
