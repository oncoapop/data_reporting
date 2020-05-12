#!/bin/sh

##########################################################
##  Script to get the gene region range from MISO output
##       Dr Damian Yap , Research Scientist
##    oncoapop@sdf.org  Version 1.0 (Jul 2015)
##########################################################

# Format of input
#1 AKT2 SE chr19 chr19:45435716:45435838:-@chr19:45433985:45434132:-@chr19:45433637:45433851:-
#2 APAF1 SE chr12 chr12:99093186:99093347:+@chr12:99097149:99097277:+@chr12:99100263:99100388:+
#3 BCL2L11 RI chr2 chr2:111881310:111881446:+@chr2:111881627:111881716:+
#4 BCL2L11 SE chr2 chr2:111881310:111881446:+@chr2:111881627:111881716:+@chr2:111921710:111926022:+
#5 BCL2L11 SE chr2 chr2:111881310:111881446:+@chr2:111907621:111907724:+@chr2:111921710:111926022:+
#6 BCL2L1 A5SS chr20 chr20:30310151:30309647|30309458:-@chr20:30252261:30253889:-

# Format of output
#(chrom start stop range gene delta) <- no header
#chr19 45433637 45435838 2201 AKT2

pwd="/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15"
cd $pwd
input="eventTable"
outfile="gene_regions.bed"

rm -f $outfile

for i in `cat $input | awk -F" " '{print $1}'`
	do
	echo $i
	chrom=`grep -w $i $input | awk -F" " '{print $4}'`
	test=(`grep -w $i $input | awk -F" " '{print $5}' | awk -F":" '{print $2,$3,$5,$6,$8,$9}' | sed -r ':a;s/([0-9])\|([0-9])/\1 \2/g;ta'`)
	echo ${test[@]}
	end=`printf "%d\n" "${test[@]}" | sort -rn | head -1`
	start=`printf "%d\n" "${test[@]}" | sort -n | head -1`

	range=`echo "$end-$start" | bc`
	gene=`grep -w $i $input | awk -F" " '{print $2}'`
	strand=`grep -w $i $input | awk -F" " '{print $5}' | awk -F":" '{print $NF}'`
	
	echo -e $chrom"\t"$start"\t"$end"\t"$range"\t"$gene"\t"$strand 
	echo -e $chrom"\t"$start"\t"$end"\t"$range"\t"$gene"\t"$strand >> $outfile

	done

exit;

