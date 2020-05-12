#!/bin/sh

# Script to count the number of bases in various sequences

cusgen="/home/dyap/Projects/Genomes/hg19_rDNA.fa"
output="/home/dyap/Projects/ChIPseq/hg19_rDNA-stats.txt"

A=`grep -o "A" $cusgen | wc -l`
T=`grep -o "T" $cusgen | wc -l`
C=`grep -o "C" $cusgen | wc -l`
G=`grep -o "G" $cusgen | wc -l`

chk=`cat $cusgen | wc -c`

echo "There are "$chk" characters in the file "$cusgen > $output
echo "------------------------------------------------------------" >> $output 

echo $A,$T,$C,$G
tot=`echo "$A+$C+$T+$G" | bc`

ratioA=`echo "(($A*100)/($tot))" | bc`
ratioT=`echo "(($T*100)/($tot))" | bc`
ratioC=`echo "(($C*100)/($tot))" | bc`
ratioG=`echo "(($G*100)/($tot))" | bc`

echo $tot
echo $ratioA
echo $ratioT
echo $ratioC
echo $ratioG

echo "Total ACTG in file =  "$tot >> $output
echo "% A in file        =  "$ratioA >> $output
echo "% T in file        =  "$ratioT >> $output
echo "% C in file        =  "$ratioC >> $output
echo "% G in file        =  "$ratioG >> $output

exit


