#!/bin/sh

projdir="/home/dyap/Projects/Takeda_T3/CG"

cd $projdir
rm -f *.tmp*

# Script to take the output of the MiSeq run A92U5
# /share/lustre/projects/takeda_splicing_inhibitor/gene_fusions/data/miseq_validation/miseq_validation_readthrough_filtered_ctrl_norm.tsv
#$1		 $2	 $3	 $4	 $5	 $6		 $7		$8	$9	$10	
#sample_id       chrom   strand  gene1   gene2   breakpt1        breakpt2     ACTB_norm        GAPDH_norm      TUBA1B_norm
#SA464   12      -       VSIG10  WSB2    118506182       118490283       0.2678053372783339    0.5344913531960581      0.3909628923232267
#SA464   22      -       PITPNB  MN1     28254375        28147084        0.00034582796979769066        0.0006902105142068331   0.0005048663506577286

input=$projdir"/miseq_validation_readthrough_filtered_ctrl_norm.tsv"

# and match that to the primer plate and summarize results
# primer plate format
# /home/dyap/Projects/PrimerDesign/Splice/positions/selected_positions.txt
# VSIG10@118506182:WSB2@118490283
# WDR77@111985275:OVGP1@111969263
# XPO4@21357728:N6AMT2@21331756
# cp /home/dyap/Projects/PrimerDesign/Splice/positions/selected_positions.txt CG_primers_ordered.txt

pri=$projdir"/CG_primers_ordered.txt"
bkpt=$projdir"/bkpts.txt"
expt=$projdir"/miseq_expt.txt"
output=$projdir"/miseq-results-by-primerset.txt"

rm -f $output

# Roadmap
# which set of primers give which Miseq results

samples="SA464 SA465 SA466 SA467 SA468 SA469 SA470 SA502 SA503 SA504 SA505 SA537 SA538 SA539 SA540"

for i in `cat $pri`
	do
	echo "######"
	echo $i
	# get gene1,2 by annotation
	gene1=`echo $i | awk -F"@" '{print $1}'`
	gene2=`echo $i | awk -F":" '{print $2}' | awk -F"@" '{print $1}'`

	# get the bk point
	stmatch=`grep $gene1 $pri`	
	enmatch=`grep $gene2 $pri`	
	match=`grep $gene1 $pri | grep $gene2`

	echo "Primers = "$match
	echo "Results ="
	echo "Primers = "$match >> $output
	echo "Results =" >> $output

	for j in $samples
		do
		match2=`grep $gene1 $input | grep $gene2 | grep $j`

		if [[ $match2 == "" ]]
			then
			echo $j" N/A N/A N/A N/A N/A N/A N/A N/A N/A" >> $output
			echo $j" N/A N/A N/A N/A N/A N/A N/A N/A N/A"
			else
			echo $match2 >> $output
			echo $match2
		fi
	
		done

	echo "======================="
	echo "=======================" >> $output

	done
exit;
