#!/bin/sh

# Name of Project
Project="CX5461"
# Project Directory

dir="/home/dyap/dyap_temp/ChIPseqAnalysis"

# each expt-ID (same as JIRA Ticket ID) (space separated, no quotes, punctuations)
for expt in 216 217 218 219

#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------
	do

	echo `date` \n >> $dir/macsjobs.log ; echo "Started summarizing for EXPT-"$expt >> $dir/macsjobs.log

	cd $dir"/"$Project"/EXPT-"$expt

	rm -f *tmp

# Taken that there is only one no drug control in all cases
# chrom   start   end     num     list    no_drug    drug_conc1   drug_conc2   .....
# column 6 = no drug
# column 7 to ncol = drug treated

# check the number of columes (depends on the drug conc treated)
	ncol=`awk -F' ' 'NR==1{print NF}' diff-IgG.bed`

# print all lines which does not have any peak in the no drug control
	awk -F" " ' $6 == 0 { print $0 }' diff-IgG.bed > IgGnodrug
	awk -F" " ' $6 == 0 { print $0 }' diff-input.bed > Inputnodrug

	for (( i=7; i <= $ncol; i++ ))
		do
		echo $ncol : $i
		awk -F"\t" ' $i == 1 { print $1 FS $2 FS $3 }' IgGnodrug > IgG.$i.tmp
		awk -F"\t" ' $i == 1 { print $1 FS $2 FS $3 }' Inputnodrug > Input.$i.tmp
		done

	echo `date` \n >> $dir/macsjobs.log ; echo "Completed EXPT-"$expt >> $dir/macsjobs.log

	done

exit
##############

