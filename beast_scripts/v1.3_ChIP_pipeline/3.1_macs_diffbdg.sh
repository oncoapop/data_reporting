#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis3/"$Project
command="bdg.sh"
unset counter
clear

for line in `cat /home/dyap/Scripts/v1.3_ChIP_pipeline/table`

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo $sample "started processing." >> $dir/macsjobs.log

	echo $line

	expt=`echo $line | awk -F, '{print $1}'`
	treat=`echo $line | awk -F, '{print $2}'`
	cont=`echo $line | awk -F, '{print $3}'`
	name=`echo $line | awk -F, '{print $4"_"$1}'`


	path=$dir"/EXPT-"$expt
	cd $path
	counter=$(($counter+1))
	if [ $counter = 1 ]
		then
			rm rep.tmp
			rm $command
	fi	
	
	treatname=`ls | grep $name | grep "_treat_pileup.bdg"`
	contname=`ls | grep $name | grep "_control_lambda.bdg"`
	echo $treatname
	echo $contname
	echo $name
	
	echo $name,$treatname,$contname >> rep.tmp


	done

	# Conditions
	cond=`cat rep.tmp | awk -F, '{print $1}' | awk -F# '{print $1}' | sort -u`
	for c in $cond
		do
		reps=`grep $c rep.tmp`
		echo $reps
		echo "++++++++++++++++++++++"

		unset counter
		count1=1

			echo "macs2 bdgdiff \\" >> $command

		for r in $reps
			do
			counter=$(($counter + $count1))
			t=`echo $r | awk -F, '{print $2}'`
			c=`echo $r | awk -F, '{print $3}'`
			n=`echo $r | awk -F, '{print $1}' | awk -F"#" '{print $1}'`

				if [ $counter -lt 3 ]
					then
					echo "--t"$counter"="$t" \\" >> $command
					echo "--c"$counter"="$c" \\" >> $command
				fi

				if [ $counter -eq 3 ]
					then
					echo "--o-prefix "$n >> $command

					echo "macs2 bdgdiff \\" >> $command
					
					echo "--t1="$t" \\" >> $command
					echo "--c1="$c" \\" >> $command
					echo "--o-prefix "$n >> $command

				fi

			done

			echo "--o-prefix "$n >> $command

		done

chmod +x $command
./$command


echo "Renaming bedgraphs..."

# rename .bdg to .bedgraph and add respective track names
for f in `ls *.bdg`
	do 
        out=`echo $f | sed 's/bdg/bedgraph/g'`
	name=`echo $out | sed 's/_treat_pileup.bedgraph//g'`
	echo $name
	header="track type=bedGraph name="$name" description="$name
	echo "$header"  > /home/dyap/dyap_temp/$f.tmp
	cat "$f" >> /home/dyap/dyap_temp/$f.tmp
    	mv /home/dyap/dyap_temp/$f.tmp "$out"
	done

	echo "==========="

echo "done"	

rm *.tmp

exit

##############
