#!/bin/sh

# This is the one that was modified to take two csv files and combine them into one file
# left and right primers

#################################
# Enter the name of sample to run 
name="DAH"
################################

# distingushing file name for match pairs
L1="leftprimers.csv"
R1="rightprimers.csv"
L2="L2.csv"
R2="R2.csv"

platedir="/home/dyap/Projects/CellMix2"

cd $platedir

fname=$name"_primer_plate"*
rm -f $fname

# Get the specific files names

LP1=`ls | grep $name | grep $L1`
RP1=`ls | grep $name | grep $R1`
LP2=`ls | grep $name | grep $L2`
RP2=`ls | grep $name | grep $R2`

echo "These are the file names to be processed:"
echo $LP1
echo $RP1
echo $LP2
echo $RP2

echo "Press enter if correct, or CTRL-C to EXIT..."

read $ans

# The is the Plate reading and compiling subroutine

for j in 1 2
	do
	# This elegant solution means we are generate variable name on the fly 
		var_name="LP${j}"
		left=`echo ${!var_name}`

		var_name="RP${j}"
		right=`echo ${!var_name}`

		plate=`echo $name"_primer_plate"$j`

	for k in A B C D E F G H 
        	do
                	for l in {1..12}
                	do
		
		        	if [[ $l -lt "10" ]];
               		 	then   {
                        		plate1="0"$l
                        		}

				else plate1=$l

                		fi

        	        primerplate=$k$plate1
	
			lpr=`grep $primerplate $left | awk -F"," '{print $3}'`
			rpr=`grep $primerplate $right | awk -F"," '{print $3}'`

			lwell=`grep $primerplate $left | awk -F"," '{print $1}'`
			rwell=`grep $primerplate $right | awk -F"," '{print $1}'`

			# names without the _L and _R suffixes (two different ways)
			lname=`grep $primerplate $left | awk -F"," '{print $2}' | awk -F"_" '{print $1"_"$2"_"$3}'`
			rname=`grep $primerplate $right | awk -F"," '{print $2}'| awk -F"_" '{print $1"_"$2"_"$3}'`

			if [[ $lname -eq $rname && $lwell -eq $rwell ]]
				then
					echo -e $lwell'\t'$lname'\t'$lpr'\t'$rpr'\t'$plate
					echo -e $lwell'\t'$lname'\t'$lpr'\t'$rpr >> $plate
			fi

			done
		done		

        done


exit;

