#!/bin/sh

# Script that ports scripts from one host to another by simply
# taking hard coded paths and make them relative paths

# hard coded home paths 
old="/meta/o/oncoapop/"

new="~"

olddir="/home/dyap/Scripts/v2_pipline"

newdir="/home/dyap/Scripts/v2.2_pipline"


for i in `ls -al $oldir`
	do
	file=`echo $i | awk -F"/" '{print $NF}'`
	cat $i | sed 's/\/meta\/o\/oncoapop\//~/g' > $newdir"/"$file

	done


