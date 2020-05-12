#!/bin/bash

# Script to split PRISM (BROAD) Data into two parts
# CX-5461 (Hong)
# T-836 (TAKEDA)

wd="/share/lustre/archive/PRISM_BROAD"
cd $wd
pdfin=`find . -name "*.pdf"`
html=`find . -name "*.html"`

# Check to find all filetypes in directory tree
# ls -R | awk -F"." '{print $2}' | sed 's/^\/.*$//g'  | sort -u

#drug="CX-5461"
drug="T-836"

OIFS="$IFS"
IFS=$'\n'

# Copying the html files whole sale
for i in `echo "$html"`
	do
	echo "-------------"
	fullname=`echo $i | sed 's/^.//g'`
	fname=`echo "$i" | rev | awk -v FS='/' '{print $1}' | rev`
	echo $fullname
	echo $fname
	dir=`echo $fullname | sed "s/\/$fname//"`
	dirname=$wd"/"$drug$dir
	cp $i $dirname"/"$fname
	done

for i in `echo "$pdfin"`
	do
	echo "-------------"
	fullname=`echo $i | sed 's/^.//g'`
	fname=`echo "$i" | rev | awk -v FS='/' '{print $1}' | rev`
	echo $fullname
	echo $fname
	dir=`echo $fullname | sed "s/\/$fname//"`
	dirname=$wd"/"$drug$dir

	# If there is the drug name then process the file name
	# If not, it must be a control file common to both, just copy it wholesale
	
	chk=`echo "$fname" | grep "$drug" | wc -l` 
	echo "=============================================>"$chk
	if [ "$chk" -ge "1" ] 
		then
			echo "RENAMING..."
			rename=`echo "$fname" | sed "s/CX5461\ and\ T-836/$drug/"`
			echo "RENAMING TO: "$rename
			cp $i $dirname"/"$rename".todo"
			rm $dirname"/"$rename
			 
		else
			echo "JUST COPYING..."
			cp $i $dirname"/"$fname
	fi

	done

IFS="$OIFS"

echo "done."

