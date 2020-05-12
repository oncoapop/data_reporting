#!/bin/sh

#################################################
##  Script to convert -hg18 files to -hg19 files
##           using chain and liftOver 
##       Dr Damian Yap , Research Scientist
##    oncoapop@sdf.org  Version 1.0 (Sep 2013)
##################################################

# modified Jul 2015 to stand alone script
# $Project (project dir that is)
# $sample group together expts if you have to split them up
# $name = each submission has to have a unique number

Project="Takeda_SpliceSignature"
sample="Sign2_HCT116_May15"
name="eventTable"

basedir="/home/dyap"
projdir=$basedir"/Projects/"$Project"/"$sample
wd=$projdir
p3dir=$projdir

liftover=$basedir"/bin/LiftOver/liftOver"
chain=$basedir"/bin/LiftOver/hg18ToHg19.over.chain"

cd $wd
# eventTable from MISO which is hg18
# Format
#No Gene Typ Chr   chr:start:end:strand@ etc
# 1 AKT2 SE chr19 chr19:45435716:45435838:-@chr19:45433985:45434132:-@chr19:45433637:45433851:-

mpfile=$wd"/"$name

# Output will always be hg19
new=$wd"/"$name"-hg19"
old=$wd"/"$name"_hg18.tmp"

clear
echo "Listing of first 10 lines of: "$mpfile
head $mpfile
rm -f $old
rm -f $old".txt"

for i in `cat $mpfile | awk -F" " '{print $1}'`
	do
	gene=`grep -w $i $mpfile | awk -F" " '{print $2}'`  
	type=`grep -w $i $mpfile | awk -F" " '{print $3}'`  

 	# This gets the # of exons (which is # of "@" + 1)
	count=`grep -w $i $mpfile | awk -F "@" '{print NF-1}'`
	exons=`echo "$count + 1" | bc`

	echo $gene"===="$exons

	for e in `seq 1 $exons`
		do
		var_name="exon${e}"
		chr=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $1}'`
		start=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $2}'`
		end=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $3}'`
		strand=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $4}'`

		if [[ $type == "A5SS" ]]
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

				echo -e $chr"\t"$start"\t"$end1"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name >> $old
				echo -e $chr"\t"$start"\t"$end1"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name 
				echo -e $chr"\t"$start"\t"$end2"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name"a" >> $old
				echo -e $chr"\t"$start"\t"$end2"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name"a" 

				unset var_name
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

				echo -e $chr"\t"$start1"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name >> $old
				echo -e $chr"\t"$start1"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name 
				echo -e $chr"\t"$start2"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name"a" >> $old
				echo -e $chr"\t"$start2"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name"a" 

				unset var_name
				continue
			fi



		fi

		echo -e $chr"\t"$start"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name >> $old
		echo -e $chr"\t"$start"\t"$end"\t"$strand"_"$i"_"$gene"_"$type"_"$var_name 

		unset var_name

		echo "here!"	
		done

		echo  "THERE"
	done

unmap=$wd"/"$name"-unmap.txt"
$liftover $old $chain $new $unmap 

echo "==================================="

echo "These are the hg18 positions:"
cat $old

echo "These are the hg19 positions:"
cat $new

echo "==================================="

echo "These are the positions which cannot be mapped:"

cat $unmap

echo "==================================="

# This section standardizes the output to match the style:

cat $new | tr "\t" "," > $new".csv"

exit;


