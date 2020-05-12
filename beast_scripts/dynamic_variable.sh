#!/bin/bash
# This works as we are using bash 3 (apparently) on beast 
# as of 23 Feb 2017


# This bash script was written to demonstrate the assignment of dynamic variables
# For each file in the directory is assigned to a unique variable 
# which can later be called upon

count=1

for i in $( ls )
	do
	file=`ls | grep $i`

	declare "magic_variable_$count=$file"

	var="magic_variable_$count"

	echo "${!var}"

	count=`echo "$count + 1" | bc`

	done

echo "This is the median file:"

	median=`echo "$count / 2" | bc`

	var="magic_variable_$median"

	echo "${!var}"
