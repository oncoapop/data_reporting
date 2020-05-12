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

start=$projdir$"/start.tmp"
end=$projdir$"/end.tmp"
anin=$projdir"/annotation_input"

# get the bkpt 1 and bkpt 2 into one file
cat $pri | awk -F"@" '{print $2,$3}' | sed 's/:.*\ /-/' > $bkpt

# positions for annotation (from output)
cat $input | tail -n +2 |  awk -F"\t" '{print "chr"$2":"$6"-"$7}' | sort -u > $anin.tmp

for h in `cat $anin.tmp`
	do
	chr=`echo $h | awk -F":" '{print $1}'`
	st=`echo $h | awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
	en=`echo $h | awk -F"-" '{print $2}'`

	stnx=`echo "$st + 1" | bc`
	ennx=`echo "$en + 1" | bc`

	echo $chr":"$st"-"$stnx >> $start
	echo $chr":"$en"-"$ennx >> $end

	done

echo "Getting WT seq for annotation..." 

twoBitToFa /home/dyap/dyap_temp/genomes/hg19.2bit start.fa -seqList=$start
twoBitToFa /home/dyap/dyap_temp/genomes/hg19.2bit end.fa -seqList=$end

cat start.fa | tr "-" " " | tr ":" " " | sed 's/chr//' | tr ">" "@" | tr "\n" " " | tr "@" "\n" | awk -F" " '{print $0" -"}'  | awk -F" " '{print $1,$2,$2,$4,$5}' > $anin-start
cat end.fa |  tr "-" " " | tr ":" " " | sed 's/chr//' | tr ">" "@" | tr "\n" " " | tr "@" "\n" | awk -F" " '{print $0" -"}'  | awk -F" " '{print $1,$2,$2,$4,$5}' > $anin-end

echo "Annotation of positions..." 

        perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $anin-start ~/bin/ANNOVAR/annovar/humandb/
        cat $anin"-start.variant_function"  | awk -F" " '{print $1, $2, $3, $4}' > $start"_anno.txt"
        perl ~/bin/ANNOVAR/annovar/annotate_variation.pl -buildver hg19 $anin-end ~/bin/ANNOVAR/annovar/humandb/
        cat $anin"-end.variant_function"  | awk -F" " '{print $1, $2, $3, $4}' > $end"_anno.txt"

output1=$start"_anno.txt"
output2=$end"_anno.txt"


for i in `cat $bkpt | sort -u`
	do
	bkpt1=`echo $i | awk -F"-" '{print $1}'`
	bkpt2=`echo $i | awk -F"-" '{print $2}'`

	# get the bk point
	stmatch=`grep $bkpt1 $output1`	
	enmatch=`grep $bkpt2 $output2`	
	match2=`grep $bkpt1 $pri | grep $bkpt2`
	name=`echo $stmatch`
	name2=`echo $enmatch`

	 if [[ $name != "" && $name2 != "NA" ]];
		then	
		{	
		echo -e $bkpt1"="$name"\t"$bkpt2"="$name2"\t"$match2 >> final.tmp
		}
		else
		{
		stagene=`echo $match2 | awk -F":" '{print $1}' | awk -F"@" '{print $1}'`
		endgene=`echo $match2 | awk -F":" '{print $2}' | awk -F"@" '{print $1}'`

		stmatch=`grep $stagene $output1`	
		enmatch=`grep $endgene $output2`	
		name=`echo $stmatch`
		name2=`echo $enmatch`
		echo -e "**************   "$bkpt1"="$name"\t"$bkpt2"="$name2"\t"$match2 >> final.tmp

		}
	
	fi


	# get the bk point from input file (ie analysis file)
	stmatch=`grep -m1 $bkpt1 $input | awk -F"\t" '{print $2,$4,$5,$6,$7}'`	
	enmatch=`grep -m1 $bkpt2 $input | awk -F"\t" '{print $2,$4,$5,$6,$7}'`	
	match2=`grep $bkpt1 $pri | grep $bkpt2`
	name=`echo $stmatch`
	name2=`echo $enmatch`

	 if [[ $name == $name2 ]];
		then	
		{	
		echo -e $name"\t"$match2 >> matching.tmp
		}
		else
		{
		echo -e "**************   "$bkpt1"="$name"\t"$bkpt2"="$name2"\t"$match2 >> matching.tmp
		}
	
	fi

	done


exit;
