#!/bin/sh

# This script goes back and returns the breakpoints of the rearrangements
# ie gets the rearrangement and sample id from the AmpliconManifest
# and then pulls out the original record's chr and pos of each 

# Working directory & files
dir="/home/dyap/Projects/Tumour_Evol/Rearrangements"

source="Destruct.AmpliconManifest.txt"
query="SVs_selected.txt"
tmp="/home/dyap/dyap_temp/destruct-tmp"
tmp2="/home/dyap/dyap_temp/destruct-tmp2"
dest="Destruct.BkPt_summary.csv"

cd $dir
rm -f $dest

cat $source | awk -F"\t" '{print $1}' | awk -F"_" '{print $1"_"$2}' > $tmp
cat $query | awk -F"," '{print $2"_"$1","$3","$4","$5","$6","$11","$12","$13","$14}' > $tmp2

	for i in `cat $tmp | sort -u`
		do
		{
		grep $i $tmp2 >> $dest
		}
	done

exit;

 
