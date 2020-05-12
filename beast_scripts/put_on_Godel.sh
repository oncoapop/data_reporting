#!/bin/sh

# This Script was written by Damian Yap (Jun 2013)

# If script is called within a programme comment these lines
# If not, $expname would be exported with export 

# godelpath="/var/www/html/workflow/primer3check/"

# expname="/home/dyap/Projects/Tumour_Evol/positions/SNV/SA495_p3_check.html"

echo "Copying" $expname

fname=`echo $expname | awk -F"/" '{print $8}'`

echo ".. to "$godelpath " on Godel"

scp $expname dyap@godel.cluster.bccrc.ca:$godelpath$fname

echo "Completed."

exit;

