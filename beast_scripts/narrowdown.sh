#!/bin/sh

# Script to narrow down the likely positions based on proximity

file="/home/dyap/Projects/Takeda_T3/CG/CG_Panel_Suppl_Table"
primerlist="/home/dyap/Projects/PrimerDesign/Splice/primer3/hct116_htert_primer_order.txt"
#$5 length
#$6 left primer
#$8 right primer
outfile="/home/dyap/Projects/Takeda_T3/CG/CG_Panel_Suppl_Table2.csv"
echo "CG_ID,Left_Primer,Right_Primer,Chr,Amp_Start,Amp_end,Amp_len,UnSplice_Len" > $outfile

for i in `cat $file  | awk -F" " '{print $2}' | sort -u`
	do
	test=`grep "$i" $file | awk -F" " '{print $2}' | wc -l`
	
	for j in `seq 1 $test`
		do
		chr=`grep -m$j "$i" $file | tail -n1 | awk -F" " '{print $1}' | awk -F">" '{print $2}' | awk -F":" '{print $1}'`
	
                plussta=`grep -m$j $i $file | tail -n1 | awk -F" " '{print $1}'| awk -F":" '{print $2}' | awk -F"+" '{print $1}'`
                plusend=`grep -m$j $i $file | tail -n1 | awk -F" " '{print $1}'| awk -F":" '{print $2}' | awk -F"+" '{print $2}'`

                minussta=`grep -m$j $i $file | tail -n1 | awk -F" " '{print $1}'| awk -F":" '{print $2}' | awk -F"-" '{print $1}'`
                minusend=`grep -m$j $i $file | tail -n1 | awk -F" " '{print $1}'| awk -F":" '{print $2}' | awk -F"-" '{print $2}'`


                if [[ $plussta > 0 ]] && [[ $plusend > 0 ]];
                        then
                        sta=$plussta
                        end=$plusend
			strand="+"
                fi

                if [[ $minussta > 0 ]] && [[ $minusend > 0 ]];
                        then
                        sta=$minussta
                        end=$minusend
			strand="-"
                fi

		gene1=`grep -m$j "$i" $file | tail -n1 | awk -F" " '{print $2}' | awk -F"@" '{print $1}'`
		sta1=`grep -m$j "$i" $file | tail -n1 | awk -F"@" '{print $2}' | awk -F":" '{print $1}'`
		gene2=`grep -m$j "$i" $file | tail -n1 | awk -F":" '{print $2}' | awk -F"@" '{print $1}'`
		end2=`grep -m$j "$i" $file | tail -n1 | awk -F"@" '{print $3}' | awk -F" " '{print $1}'`
		
		unsplen=`grep -m$j "$i" $file | tail -n1 | awk -F" " '{print $3}'`
		left=`grep -m$j "$i" $file | tail -n1 | awk -F" " '{print $4}'`
		right=`grep -m$j "$i" $file | tail -n1 | awk -F" " '{print $5}'`

		# check
		chkname=`grep "$i" $primerlist  | awk -F, '{print $2}' | sed 's/chr//'`
		chkleft=`grep "$i" $primerlist  | awk -F, '{print $6}'`
		chkright=`grep "$i" $primerlist  | awk -F, '{print $8}'`
		splen=`grep "$i" $primerlist  | awk -F, '{print $5}'`
		
		if [[ $i != $chkname && $left != $chkleft && $right != $chkright ]]
			then
			continue
		fi

		if [[ $sta == $sta1 && $end == $end2 ]]
			then
			echo $i,$left,$right,$chr,$sta,$end,$strand,$splen,$unsplen >> $outfile
			continue
		fi
 
		diff1=`echo "$sta1-$sta" | bc | tr -d "-"`
		diff2=`echo "$end-$end2" | bc | tr -d "-"`

		echo $diff1,$diff2

		if (( $diff1 < 3000000 && $diff2 < 3000000 ))
			then
				echo $i,$left,$right,$chr,$sta,$end,$strand,$splen,$unsplen >> $outfile
				echo "here"
				continue
		fi
				echo $i,$left,$right,$chr,$sta,$end,$strand,$splen,$unsplen

		done 
	done
