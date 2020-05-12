#!/bin/sh

#################################################
##  Script to convert -hg18 files to -hg19 files
##       and get junctions and seq in hg19
##       Dr Damian Yap , Research Scientist
##    dyap@bccrc.ca  Version 1.0 (Aug 2015)
##################################################

# modified Aug 2015 to get junction sequences
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
junc=$wd"/"$name"_junct.tmp"


clear
#echo "Listing of first 10 lines of: "$mpfile
#head $mpfile
rm -f $old
rm -f $old".txt"
rm -f $junc

# Testing if the S/No is unique
	test1=`cat $mpfile | awk -F" " '{print $1}' | wc -l`
	test2=`cat $mpfile | awk -F" " '{print $1}' | sort -u | wc -l`
	if [[  "$test1" != "$test2" ]]
		then
		echo "S/No not unique" 
		exit 1;
		
	fi

for i in `cat $mpfile | awk -F" " '{print $1}'`
	do
	gene=`grep -w $i $mpfile | awk -F" " '{print $2}'`  
	type=`grep -w $i $mpfile | awk -F" " '{print $3}'`  

 	# This gets the # of exons (which is # of "@" + 1)
	count=`grep -w $i $mpfile | awk -F "@" '{print NF-1}'`
	exons=`echo "$count + 1" | bc`

# Getting junction coordinates (or exons coordinates for single exon regions)
# single exon regions 
	if [[ $exons == 1 ]]
		then 
		jchr=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$exons -F"@" '{print $exon}' | awk -F":" '{print $1}'`
		jstart=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$exons -F"@" '{print $exon}' | awk -F":" '{print $2}'`
		jend=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$exons -F"@" '{print $exon}' | awk -F":" '{print $3}'`
		jstrand=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$exons -F"@" '{print $exon}' | awk -F":" '{print $4}'`
		jregion=`echo "$jend - $jstart + 1" | bc`

		echo -e $jchr"\t"$jstart"\t"$jend"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_REFINE" >> $junc
		echo -e $jchr"\t"$jstart"\t"$jend"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_REFINE" 
	fi

# multiple exon regions (junctions separated by @)
	if [[ $exons > 1 ]]
		then 
		for e in `seq 1 $exons`
			do
			var_name="exon${e}"
			jchr1=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $1}'`
			jstart1=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $2}'`
			jend1=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $3}'`
			jstrand1=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$e -F"@" '{print $exon}' | awk -F":" '{print $4}'`

			f=`echo "$e + 1" | bc`
			echo $e,$f

			if [[ $f > $exons ]]
				then
					echo "continue..."
					continue
			fi

			var_name="exon${f}"
			echo "same loop"
			jchr2=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$f -F"@" '{print $exon}' | awk -F":" '{print $1}'`
			jstart2=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$f -F"@" '{print $exon}' | awk -F":" '{print $2}'`
			jend2=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$f -F"@" '{print $exon}' | awk -F":" '{print $3}'`
			jstrand2=`grep -w $i $mpfile | awk -F" " '{print $5}' | awk -v exon=$f -F"@" '{print $exon}' | awk -F":" '{print $4}'`


			if [[ $jchr1 == $jchr2 && $jend1 < $jstart2 ]]
				then 
				echo "plus"
				jchr=$jchr1
				jstart=$jend1
				jend=$jstart2
				jstrand=$jstrand1
							
				jleft=`echo "$jstart - 50" | bc`
				jright=`echo "$jend + 50" | bc`


				if [[ $jleft > $jstart1 && $jright < $jend2 ]]
					then

					jregion=`echo "($jright - $jend) +  ($jstart - $jleft)  + 1" | bc`
					echo "($jright - $jend) +  ($jstart - $jleft)  + 1" | bc
					echo "($jright - $jend)" | bc
					echo "($jstart - $jleft)" | bc

					echo -e $jchr"\t"$jleft"\t"$jstart"\t"$jend"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type 
					echo -e $jchr"\t"$jleft"\t"$jstart"\t"$jend"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type >> $junc

					elif [[ $jleft < $jstart ]]
						then
						jleft=$jstart
					
					jregion=`echo "($jright - $jend) +  ($jstart - $jleft)  + 1" | bc`
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type

					elif [[ $jright > $jend ]]
						then
						jright=$jend
					
					jregion=`echo "($jright - $jend) +  ($jstart - $jleft)  + 1" | bc`
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type
					
					else

					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_OUT-OF-RANGE" >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_OUT-OF-RANGE" 
				fi
			fi

			if [[ $jchr1 == $jchr2 && $jstart2 < $jend1 ]]
				then 
				echo "minus"
				jchr=$jchr1
				jstart=$jstart2
				jend=$jend1
				jstrand=$jstrand1
				unset jregion

				jleft=`echo "$jstart - 50" | bc`
				jright=`echo "$jend + 50" | bc`
echo "any here?"							
				if [[ $jleft > $jstart1 && $jright > $jend2 ]]
					then

					jregion=`echo "($jstart - $jright) + ($jleft - $jend) + 1" | bc`
					echo "($jstart - $jright)" | bc
					echo "($jleft - $jend)" | bc
echo "hello"
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type

					if [[ $jleft < $jstart ]]
						then
						$jleft=$jstart
					fi
					if [[ $jright > $jend ]]
						then
						$jright=$jend
					fi

					jregion=`echo "($jstart - $jright) + ($jleft - $jend) + 1" | bc`

					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type

					else
echo "oh no!"					
					echo "($jstart - $jright)" | bc
					echo "($jleft - $jend)" | bc

					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_OUT-OF-RANGE" >> $junc
					echo -e $jchr"\t"$jleft"\t"$jright"\t"$jstrand"\t"$jregion"\t"$i"_"$gene"_"$type"_OUT-OF-RANGE" 
				fi
			fi

			done
	
	fi

	done
echo $junc
#cat $junc
exit;


##################
	
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


