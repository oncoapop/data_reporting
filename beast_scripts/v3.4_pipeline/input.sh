#!/bin/sh

# I took the NM_ (REF SEQ ID that siRNA is targed against)
# used biomart to get the ESNT and ESNP ID
# script is query.pl (actually used the online tool)

inpath="/home/dyap/Projects/PrimerDesign/manual"
infile=$inpath"/siRNA_CGBP-set-ENSP"

outpath="/home/dyap/Projects/PrimerDesign/manual"
outfile=$outpath"/assoc-factors_seq.txt"

rm -f $outfile

genomepath="/home/dyap/dyap_temp/genomes"
genomefile=$genomepath"/gencode.v19.pc_transcripts.fa"

for i in `cat $infile  | awk '{print $1}'`
	do
	echo $i
	grep -A1 $i $genomefile >> $outfile 
	grep -A1 $i $genomefile 
	done

echo "Records in input file:"
cat $infile | wc -l

echo "Records in output file:" 
cat $outfile | grep -c ">"

