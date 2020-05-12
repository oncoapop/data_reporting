#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/Projects/ChIPseqAnalysis"

# Type of Analysis
# type="summits.bed"
# type="peaks.narrowPeak"
type="peaks.broadPeak"


# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 216 217 218 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	 echo `date` \n >> $dir/macsjobs.log ; echo $sample "started Summary." >> $dir/macsjobs.log

	cd $dir"/"$Project"/EXPT-"$expt

	rm -f *.sorted
	rm -f diff*
	rm -f drug*

	echo "clear"

	echo "==========================="
	echo $expt

# sorted bed files are required
# sort by chromosome then position
# Use narrow peak rather than summit (which is single base peaks only)
# next to try is broadpeak but need to rerun macs2

	for i in `ls | grep "$type"`
		do
# bed files need to be sorted (narrow peaks is a bed format + extra info about peaks)
		sort -k1,1 -k2,2n $i > $i.sorted
		done

# This section gets all the relevant file names (depending on treatment)
	drugIg=`ls | grep "CX_"$expt"_"$type".sorted" | grep "IgG"`
	druginp=`ls | grep "CX_"$expt"_"$type".sorted" | grep "input"`

	
	Ig=`ls | grep "nodrug_"$expt"_"$type".sorted" | grep "IgG"`
	inp=`ls | grep "nodrug_"$expt"_"$type".sorted" | grep "input"`

	IgG=$Ig" "$drugIg
	GL=`echo "$IgG" | wc -w`
	echo $GL

	Input=$inp" "$druginp
	IL=`echo "$Input" | wc -w`
	echo $IL

# If there are >4 files then add accordingly
# Assigns each file name to a variable

	file1=`echo $IgG | cut -d" " -f1`
	file2=`echo $IgG | cut -d" " -f2`
	file3=`echo $IgG | cut -d" " -f3`
	file4=`echo $IgG | cut -d" " -f4`

	file11=`echo $Input | cut -d" " -f1`
	file12=`echo $Input | cut -d" " -f2`
	file13=`echo $Input | cut -d" " -f3`
	file14=`echo $Input | cut -d" " -f4`
	
# This just prints out all the file names to check that you have all the files
	for t in $file1 $file2 $file3 $file4 $file11 $file12 $file13 $file14
		do
		echo t=$t
		done

# This is the main command to get all the interects from multiple bed files

	/home/dyap/bin/bedtools2/bin/multiIntersectBed -header -i $file1 $file2 $file3 $file4 > diff-IgG.bed

	if [ $GL == 4 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 && $9 == 1 )  print $1,$2,$3 }' diff-IgG.bed > drug-treat_regions_$expt.bed
	fi

	if [ $GL == 3 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 )  print $1,$2,$3 }' diff-IgG.bed > drug-treat_regions_$expt.bed
	fi

	if [ $GL == 2 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 )  print $1,$2,$3 }' diff-IgG.bed > drug-treat_regions_$expt.bed
	fi
	
	if [ $GL == 1 ];
		then
		echo "Nothing to compare!"
		exit;	
	fi

	/home/dyap/bin/bedtools2/bin/multiIntersectBed -header -i $file11 $file12 $file13 $file14 > diff-input.bed

	if [ $IL == 4 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 && $9 == 1 )  print $1,$2,$3 }' diff-input.bed > drug-treat_regions2_$expt.bed
	fi

	if [ $IL == 3 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 && $8 == 1 )  print $1,$2,$3 }' diff-input.bed > drug-treat_regions2_$expt.bed
	fi

	if [ $IL == 2 ];
		then
		awk -F"\t" '{ if ($6 == 0 && $7 == 1 )  print $1,$2,$3 }' diff-input.bed > drug-treat_regions2_$expt.bed
	fi
	
	if [ $IL == 1 ];
		then
		echo "Nothing to compare!"
		exit;	
	fi

	cp  drug-treat_regions2_$expt.bed ~/Projects/ChIPseqAnalysis/CX5461/drug-treat-regions-input-$expt.bed
	cp  drug-treat_regions_$expt.bed ~/Projects/ChIPseqAnalysis/CX5461/drug-treat-regions-IgG-$expt.bed

	done

exit

##############

