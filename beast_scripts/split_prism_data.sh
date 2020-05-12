#!/bin/bash

# Script to split PRISM (BROAD) Data into two parts
# CX-5461 (Hong)
# T-836 (TAKEDA)

wd="/share/lustre/archive/PRISM_BROAD"
cd $wd
files=`find . -name "*.csv"`

drug="CX-5461"

OIFS="$IFS"
IFS=$'\n'
for i in `echo "$files"`
	do
	echo "-------------"
	fullname=`echo $i | sed 's/^.//g'`
	fname=`echo "$i" | rev | awk -v FS='/' '{print $1}' | rev`
	echo $fullname
	echo $fname
	dir=`echo $fullname | sed "s/\/$fname//"`
	dirname=$wd"/"$drug$dir

	if [ ! -d $dirname ]
		then
		mkdir -p $dirname
	fi


	# If there is the drug name then split file 
	# If not, it must be a control file common to both, just copy it wholesale

	chk=`cat "$i" | grep "$drug" | wc -l` 
	if [ "$chk" -gt "1" ] 
		then
			cat $i | grep $drug > $dirname"/"$fname
		else
			cp $i $dirname"/"$fname
	fi

	done
IFS="$OIFS"
