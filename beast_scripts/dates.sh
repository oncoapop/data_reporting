#!/bin/bash 
# A program that finds the number of days between two dates in the same year. 
# The user inputs one date, presses the return key, then enters another date. 
# The input is of the form "Month day", for example, Jul 4 or Oct 31. 

jan=( $(seq 0 31) ) 
feb=( $(seq 32 60) ) 
mar=( $(seq 61 92) ) 
apr=( $(seq 93 123) ) 
may=( $(seq 124 155) ) 
jun=( $(seq 156 186) ) 
jul=( $(seq 187 218) ) 
aug=( $(seq 219 250) ) 
sep=( $(seq 251 281) ) 
oct=( $(seq 282 313) ) 
nov=( $(seq 314 344) ) 
dec=( $(seq 345 376) ) 

echo    "Enter a month (all small) and day within a (normal, not leap) year: (ie jan 1 & dec 31) " 
read mo_1 da_1 
 # a case for every month 
case $mo_1 in     
         # put into date_1 the value from corresponding array above 
    jan) date_1=${jan[$da_1]} ;; 
    feb) date_1=${feb[$da_1]} ;; 
    mar) date_1=${mar[$da_1]} ;; 
    apr) date_1=${apr[$da_1]} ;; 
    may) date_1=${may[$da_1]} ;; 
    jun) date_1=${jun[$da_1]} ;; 
    jul) date_1=${jul[$da_1]} ;; 
    aug) date_1=${aug[$da_1]} ;; 
    sep) date_1=${sep[$da_1]} ;; 
    oct) date_1=${oct[$da_1]} ;; 
    nov) date_1=${nov[$da_1]} ;; 
    dec) date_1=${dec[$da_1]} ;; 
      *) please enter date format "jan 1"   
esac 
echo '   Enter a second month and day: ' 
read mo_2 da_2 
case $mo_2 in 
    jan) date_2=${jan[$da_2]} ;; 
    feb) date_2=${feb[$da_2]} ;; 
    mar) date_2=${mar[$da_2]} ;; 
    apr) date_2=${apr[$da_2]} ;; 
    may) date_2=${may[$da_2]} ;; 
    jun) date_2=${jun[$da_2]} ;; 
    jul) date_2=${jul[$da_2]} ;; 
    aug) date_2=${aug[$da_2]} ;; 
    sep) date_2=${sep[$da_2]} ;; 
    oct) date_2=${oct[$da_2]} ;; 
    nov) date_2=${nov[$da_2]} ;; 
    dec) date_2=${dec[$da_2]} ;; 
      *) please enter date format "jan 1" 
esac 
if [[ $date_2 -gt $date_1 ]] 
    then number=$(($date_2 - $date_1)) 
    else number=$(($date_1 - $date_2)) 
fi 
echo "The number of days between dates is " $number 

