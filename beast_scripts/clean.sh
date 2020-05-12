#!/bin/sh

project="FL"

clear

echo "This script is to clean all the files when testing primer design pipelines."
echo "---------------------------------------------------------------------------"

echo "All Project "$project" files will be deleted!!!"
echo "Press ENTER to continue. Ctrl-C to exit"
read ans

if [[ $project == "" ]]
	then echo "No Project specified - exiting to avoid mass deletion!!"
	exit 1;
fi

# Deleting files that are left by v4.1 primer design pipeline

rm -f "/home/dyap/dyap_temp/"$project*

rm -fr "/home/dyap/Projects/PrimerDesign/"$project"/primer3"
rm -fr "/home/dyap/Projects/PrimerDesign/"$project"/positions"
rm -fr "/home/dyap/Projects/PrimerDesign/"$project"/Annotate"

