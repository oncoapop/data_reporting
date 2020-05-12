#!/bin/sh

# This is the one that was modified to take two csv files and combine them into one file
# left and right primers

#################################
# Enter the name of sample to run 
name="SA535"
################################

# distingushing file name for match pairs
L1="leftprimers.csv"
R1="rightprimers.csv"
L2="L2.csv"
R2="R2.csv"
IDTfile="coa.csv"

platedir="/home/dyap/Projects/Tumour_Evol/Xenodrug/"$name

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
	
			# This get the left and right primer sequences

			lpr=`grep $primerplate $left | awk -F"," '{print $3}'`
			rpr=`grep $primerplate $right | awk -F"," '{print $3}'`

			lwell=`grep $primerplate $left | awk -F"," '{print $1}'`
			rwell=`grep $primerplate $right | awk -F"," '{print $1}'`

			# names without the _L and _R suffixes (two different ways)
			lname=`grep $primerplate $left | awk -F"," '{print $2}' | awk -F"_" '{print $1"_"$2"_"$3}'`
			rname=`grep $primerplate $right | awk -F"," '{print $2}'| awk -F"_" '{print $1"_"$2"_"$3}'`

			# names without the _L and _R suffixes (two different ways)
			nameL=`grep $primerplate $left | awk -F"," '{print $2}' | awk -F"_" '{print $1"_"$2"_"$3"_"$4}'`
			nameR=`grep $primerplate $right | awk -F"," '{print $2}'| awk -F"_" '{print $1"_"$2"_"$3"_"$4}'`

			# Primer name (L or R) and Seq from IDT file 
			OrderL=`grep "$nameL"  $IDTfile | awk -F, '{print $6}' | tr -d '"'`
			OrdSeqL=`grep "$nameL" $IDTfile | awk -F, '{print $10}' | tr -d " " | tr -d '"'`
			OrderR=`grep "$nameR"  $IDTfile | awk -F, '{print $6}' | tr -d '"'`
			OrdSeqR=`grep "$nameR" $IDTfile | awk -F, '{print $10}' | tr -d " " | tr -d '"'`


#			echo $nameL,$OrderL,$lpr,$OrdSeqL		
			if [[ $nameL == $OrderL && $lpr == $OrdSeqL ]]
				then 
					LStat="Left_Ordered"
			fi

#			echo $nameR,$OrderR,$rpr,$OrdSeqR		
			if [[ $nameR == $OrderR && $rpr == $OrdSeqR ]]
				then 
					RStat="Right_Ordered"
			fi

			if [[ $lname == $rname && $lwell == $rwell ]]
				then
					echo -e $lwell'\t'$lname'\t'$lpr'\t'$rpr'\t'$LStat'\t'$RStat'\t'$plate
					echo -e $lwell'\t'$lname'\t'$lpr'\t'$rpr'\t'$LStat'\t'$RStat'\t'$plate >> $plate
			fi

			# Clear the Status Flag
			LStat=""
			RStat=""

			done
		done		

        done


exit;

