#!/bin/sh

#################################################
##  Script to convert output to exon bed format
##       Dr Damian Yap , Research Scientist
##    oncoapop@sdf.org  Version 1.0 (Jul 2015)
##################################################

# Format of input
#1 AKT2 SE chr19 chr19:45435716:45435838:-@chr19:45433985:45434132:-@chr19:45433637:45433851:-
#2 APAF1 SE chr12 chr12:99093186:99093347:+@chr12:99097149:99097277:+@chr12:99100263:99100388:+
#3 BCL2L11 RI chr2 chr2:111881310:111881446:+@chr2:111881627:111881716:+
#4 BCL2L11 SE chr2 chr2:111881310:111881446:+@chr2:111881627:111881716:+@chr2:111921710:111926022:+
#5 BCL2L11 SE chr2 chr2:111881310:111881446:+@chr2:111907621:111907724:+@chr2:111921710:111926022:+

# Format of output
#chrom start stop gene score strand type
#1 chr15 73062668 73062770 ENST00000362710 0.000000 -1 exon

pwd="/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15"
cd $pwd
input="hg19_custom.gff3"
ref="eventTable"
outfile1="gene_regions.bed"
outfile2="exon_regions.bed"
outfile3="mRNA_regions.bed"

rm -f $outfile
rm -f $outfile1
rm -f $outfile2
rm -f $outfile3

# For Gene list
for i in `cat $ref | awk -F" " '{print $5}'`
	do
	echo $i
	chrom=`grep -w $i $input | grep "gene" | awk -F"\t" '{print $1}'`
	start=`grep -w $i $input | grep "gene" | awk -F"\t" '{print $4}'`
	end=`grep -w $i $input | grep "gene" | awk -F"\t" '{print $5}'`
	gene=`grep -w $i $ref  | awk -F" " '{print $2}'`
	strand=`grep -w $i $ref | awk -F" " '{print $5}' | awk -F":" '{print $NF}'`
	type=`grep -w $i $input | grep "gene" | awk -F"\t" '{print $3}'`
	
	echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type 
	echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile1
	done

# For exon
for i in `cat $ref | awk -F" " '{print $5}' `
	do
	k=`grep "exon" $input | grep $i | awk -F"\t" '{print $1}'`
	l=`echo ${k[@]} | awk -F" " '{print NF}'`
	m=1
	for j in $(eval echo "{$m..$l}")
		do
		chrom=`grep "exon" $input | grep -w -m$j $i | tail -n1 | awk -F"\t" '{print $1}'`
		start=`grep "exon" $input | grep -w -m$j $i | tail -n1 | awk -F"\t" '{print $4}'`
		end=`grep "exon" $input | grep -w -m$j $i | tail -n1 | awk -F"\t" '{print $5}'`
		name=`grep -w $i $ref  | awk -F" " '{print $2}'`
		strand=`grep -w $i $ref | awk -F" " '{print $5}' | awk -F":" '{print $NF}'`
		type=`grep "exon" $input | grep -w -m$j $i | tail -n1 | awk -F"\t" '{print $3}'`
		event=`grep "exon" $input | grep -w -m$j $i | tail -n1 | awk -F"\t" '{print $2}'`
		gene=$name"_"$event

                if [[ $event == "A5SS" ]]
                        then
                        alt=`echo $end | awk -F "\|" '{print NF-1}'`

                        if [[ $alt == "1" ]]
                                then
                                end1=`echo $end | awk -F"\|" '{print $1}'`
                                end2=`echo $end | awk -F"\|" '{print $2}'`

                                if [[ $start > $end1 ]]
                                        then
                                        begin=$end1
                                        stop=$start

                                        start=$begin
                                        end1=$stop

                                fi

                                if [[ $start > $end2 ]]
                                        then
                                        begin=$end2
                                        stop=$start

                                        start=$begin
                                        end2=$stop


                                fi

				echo -e $chrom"\t"$start"\t"$end1"\t"$gene"\t0\t"$strand"1\t"$type 
				echo -e $chrom"\t"$start"\t"$end1"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile2

				echo -e $chrom"\t"$start"\t"$end2"\t"$gene"\t0\t"$strand"1\t"$type 
				echo -e $chrom"\t"$start"\t"$end2"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile2

                                continue
                        fi
                fi

                if [[ $type == "A3SS" ]]
                        then
                        alt=`echo $start | awk -F "\|" '{print NF-1}'`

                        if [[ $alt == "1" ]]
                                then
                                start1=`echo $start | awk -F"\|" '{print $1}'`
                                start2=`echo $start | awk -F"\|" '{print $2}'`

                                if [[ $start1 > $end ]]
                                        then
                                        begin=$end 
                                        stop=$start1

                                        start1=$begin
                                        end=$stop

                                fi

                                if [[ $start2 > $end ]]
                                        then
                                        begin=$end 
                                        stop=$start2

                                        start2=$begin
                                        end=$stop


                                fi

				echo -e $chrom"\t"$start1"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type 
				echo -e $chrom"\t"$start2"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile2

				echo -e $chrom"\t"$start1"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type 
				echo -e $chrom"\t"$start2"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile2

                                continue
                        fi



                fi

		echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type 
		echo -e $chrom"\t"$start"\t"$end"\t"$gene"\t0\t"$strand"1\t"$type >> $outfile2
		done
	done

# For mRNA

# Copy to Roche folder for R script use
#cat $outfile1 | sort -u > "Roche_150725_HG19_Splice-Sign2_RNA/"$outfile1
#cat $outfile2 | sort -u > "Roche_150725_HG19_Splice-Sign2_RNA/"$outfile2
#cat $outfile1 | sort -u > "Roche_150727_HG19_Splice-Sign2a_RNA/"$outfile1
#cat $outfile2 | sort -u > "Roche_150727_HG19_Splice-Sign2a_RNA/"$outfile2
cat $outfile1 | sort -u > "Roche_150729_custom/"$outfile1
cat $outfile2 | sort -u > "Roche_150729_custom/"$outfile2



exit;

