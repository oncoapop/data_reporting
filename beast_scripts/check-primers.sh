#!/bin/sh

primerlist="/home/dyap/Projects/PrimerDesign/Splice/primer3/hct116_htert_primer_order.txt"
#$5 length
#$6 left primer
#$8 right primer

primer_order="/home/dyap/Projects/PrimerDesign/Splice/primer3/hct116_htert-filtered.AmpliconManifest"

temp="/home/dyap/dyap_temp/CG_Panel_primers.csv"
temp2="/home/dyap/dyap_temp/CG_Panel_isPCR"

rm -f $temp

outfile="/home/dyap/Projects/Takeda_T3/CG/CG_Panel_Suppl_Table"
rm -f $outfile

for i in `cat $primer_order | grep  @ | awk -F"\t" '{print $2}'`
	do
	length=`grep -m1 "$i" $primerlist | awk -F"," '{print $5}'`
	left=`grep -m1 "$i" $primerlist | awk -F"," '{print $6}'`
	right=`grep -m1 "$i" $primerlist | awk -F"," '{print $8}'`
	echo $i","$left","$right","$length | tr "," "\t" >> $temp
	done

# isPCR is on beast at 
command="/share/data/apps/isPcr/bin/x86_64/isPcr"

# database (hg19 2bit fa ) at
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"
#database="/home/dyap/Projects/PrimerDesign/manual/gp140.2bit"
#database="/home/dyap/Projects/PrimerDesign/manual/"$name".2bit"

# IF reversecomplement of right primer is NOT required comment this
#flip="-flipReverse"
flip=""

# output format
output=fa       # fasta format (default)
#output=bed     # bed format (tab-delimited; Fields: chrom/start/end/name/score/strand)
#output=psl     # blat format 

# Name of the input file
inputfile=$temp

# Name of the output file
outputfile=$temp2

cat $inputfile
echo $outputfile
$command $database $flip "-out="$output $inputfile $outputfile -maxSize=100000000
echo $command" "$database" " $flip "-out="$output" " $inputfile" " $outputfile" -maxSize=100000000"

grep "@" $outputfile > $outfile

exit
