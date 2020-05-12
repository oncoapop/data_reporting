#!/bin/sh

# This is a script to compare mouse contamination of MiSeq reads from Peter's Xenograft reads
# This requires 1. Blast outout of amplicons blasted against mouse genome
# blastn -db  Mus_musculus.GRCm38.70.dna.chromosome.1.fa -query SA495_p3_amplicons.fa -out SA495_blast-amp.txt

# This also requires the output file from the bionomal test
# /share/lustre/projects/breast_xeno_evolution/binomial_validation/xeno_SA495_binomvalid_results_table

source="/share/lustre/projects/breast_xeno_evolution/binomial_validation/xeno_SA495_binomvalid_results_table"

query="/home/dyap/Projects/Tumour_Evol/positions/SNV/SA495_blast-amp.txt"

outdir="/home/dyap/Projects/Tumour_Evol/SA495/QC"

outfile=$outdir"/mouse-reads.txt"

tmp="/home/dyap/dyap_temp"

# This section prints the sum, avg of the selected cols per row (section is 8 lines)
cat  $source | awk -F"\t" '{print $5,$6}' | awk '
BEGIN {FS=OFS=" "}
{
sum=0; n=0
for(i=1;i<=NF;i++)
   {sum+=$i; ++n}
   print $0,"sum:"sum,"count:"n,"avg:"sum/n
}' > $tmp/sum.tmp

# This section prints the total sum of the column (numerical value) of the columns
echo "This is the total reads of the NSG sample:"
echo $(( $( cat $tmp"/sum.tmp"  | grep -o sum:[0-9]* | sed 's/sum://' | tr "\n" "+"  | xargs -I{} echo {} 0  ) ))

echo "This is the total number of NSG positions called:"
cat $tmp/sum.tmp | wc -l

exit

