#!/bin/sh

projdir="/home/dyap/Projects/Takeda_T3/CG"

cd $projdir
rm -f *.tmp*

# Script to take the output of the MiSeq run A92U5
# /share/lustre/projects/takeda_splicing_inhibitor/gene_fusions/data/miseq_validation/miseq_validation_readthrough_filtered_ctrl_norm.tsv
# $1		$2	$3	$4	$5	$6		$7		$8		$9		$10	
# sample_id     chrom   strand  gene1   gene2   breakpt1        breakpt2        ACTB_norm       GAPDH_norm      TUBA1B_norm
# SA464   12      -       ENSG00000176834 ENSG00000176871 118506182       118490283       0.2678053372783339      0.5344913531960581      0.3909628923232267
# SA464   22      -       ENSG00000180957 ENSG00000169184 28254375        28147084        0.00034582796979769066  0.0006902105142068331   0.0005048663506577286
output=$projdir"/miseq_validation_readthrough_filtered_ctrl_norm.tsv"

# and match that to the primer plate and summarize results
# primer plate format
# /home/dyap/Projects/PrimerDesign/Splice/positions/selected_positions.txt
# VSIG10@118506182:WSB2@118490283
# WDR77@111985275:OVGP1@111969263
# XPO4@21357728:N6AMT2@21331756
# cp /home/dyap/Projects/PrimerDesign/Splice/positions/selected_positions.txt CG_primers_ordered.txt

pri=$projdir"/CG_primers_ordered.txt"
bkpt=$projdir"/bkpts.txt"

start=$projdir$"/start.tmp"
end=$projdir$"/end.tmp"

# get the bkpt 1 and bkpt 2 into one file
cat $pri | awk -F"@" '{print $2,$3}' | sed 's/:.*\ /-/' > $bkpt

for i in `cat $bkpt`
	do
	bkpt1=`echo $i | awk -F"-" '{print $1}'`
	bkpt2=`echo $i | awk -F"-" '{print $2}'`

        # Process the annotation file and runs ANNOVAR
        st=`cat $output | awk -F"\t" '{print $6}'`
        en=`cat $output | awk -F"\t" '{print $7}'`

	stnx=`echo "$st + 1" | bc`
	ennx=`echo "$en + 1" | bc`

        cat $output | awk -F"\t" '{print $2" "$6" "$stnx}' >> $start
        cat $output | awk -F"\t" '{print $2" "$7" "$ennx}' >> $end

	# get the bk point
	match=`grep $bkpt1 $output | grep $bkpt2`	
	match2=`grep $bkpt1 $pri | grep $bkpt2`
	name=`echo $match | awk -F" " '{print $2,$3,$4,$5,$6,$7}'`
	
	echo -e $name"\t"$match2
	done


        perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $start ~/bin/ANNOVAR/annovar/humandb/
        cat $start".variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $start"_anno.txt"

        perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $end ~/bin/ANNOVAR/annovar/humandb/
        cat $end".variant_function"  | awk -F" " '{print "Chr"$3 "_" $4, $6, $1, $2}' > $end"_anno.txt"



exit;
