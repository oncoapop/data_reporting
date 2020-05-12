#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)

clear 
cd $dir

# This tests to see if there are any primers returned, if not, the record is skipped
#		test=`grep -A19 Chr19_11943155 ~/Projects/Single_Cell/positions/SNV/SA212_p3_output.txt | grep "PRIMER_PAIR_EXPLAIN="` 
test="something"
#		test=`grep -A19 xxx_11943155 ~/Projects/Single_Cell/positions/SNV/SA212_p3_output.txt | grep "PRIMER_PAIR_EXPLAIN="` 
		echo $test
#			if [[ $test =~ "ok" && $test =~ "0" ]];
			if [[ $test =~ "" ]];
			
				then
					echo $j "blank";
			
				else 
					echo $j "not blank";
					echo $j;
			fi

		if [ -z "$test" ]
then
	echo "Empty!"
 
else
	echo "This is the value of test= "$test
fi		

exit;
