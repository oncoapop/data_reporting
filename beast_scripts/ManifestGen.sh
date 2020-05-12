#/bin/sh
# SCript to generate the Manifest for VBA0038 for reanalysis on Basespace
# This is based on the ______ workflow and is based on the cancer panel manifest

dir="/share/lustre/backup/dyap/Projects/Single_Cell/positions/Annotate/"

annfile=$dir"VBA0038_Annotated_SNVpos.txt"
manfile=$dir"VBA0038SNV.AmpliconManifest"

	for i in A B C D
	do
		for j in 01 02 03 04 05 06 07 08 09 10 11 12
		do
		
		echo $i$j

		line=`grep $i$j $manfile | awk -F"\t" '{print $2, $3, $4}'`
		chr=`grep $i$j $manfile | awk -F"\t" '{print $2}'`
		start=`grep $i$j $manfile | awk -F"\t" '{print $3}'`
		end=`grep $i$j $manfile | awk -F"\t" '{print $4}'`

		ann_start=`grep $start $annfile | awk -F"\t" '{print $6}'`
		ann_end=`grep $end $annfile | awk -F"\t" '{print $6}'`

		echo $chr $start $end "annotated=" $ann_start $ann_end
		done

	echo $i "series done."
	done

