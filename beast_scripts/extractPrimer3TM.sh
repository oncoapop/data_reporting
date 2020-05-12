#!/bin/sh

cd /home/dyap/dyap_temp

infile="2_primer3_view.txt"
pairnumber="temp_file1.txt"
outfile1="primers_list_temp1.csv"
outfile2="primers_list_temp2.csv"
outfile3="primers_list_temp3.csv"
outfile4="primers_list_temp4.csv"
outfile="primers_list.csv"
rm $outfile
rm $outfile1
rm $outfile3
rm $outfile2
rm $outfile4

rm $pairnumber

grep -e '^PRIMER PICKING' -e '^LEFT PRIMER' -e '^RIGHT PRIMER' -e '^SEQUENCE SIZE' -e '^PRODUCT SIZE' $infile > $pairnumber

# replace phrases in $pairnumber: "PRIMER PICKING RESULTS FOR " with ":", remove extra ,,,
sed -e 's/PRIMER PICKING RESULTS FOR /@/g' -e 's/SEQUENCE SIZE: /SEQUENCE_SIZE,/g' -e 's/LEFT PRIMER/LEFT_PRIMER/g' -e 's/RIGHT PRIMER/RIGHT_PRIMER/g' -e 's/PRODUCT SIZE: /PRODUCT_SIZE,/g' -e 's/COMPL:/COMPL/g' -e 's/ \{1,\}/,/g' $pairnumber > $outfile1

# move all the rows to one row
tr '\n' ',' < $outfile1 > $outfile2

# add a new line where ":" precedes "chr". This will put all ID parameters in separate rows.
tr '@' '\n' < $outfile2 > $outfile3

sed -e 's/LEFT_PRIMER/,,LEFT_PRIMER/g' $outfile3 > $outfile4

# print $1=ID, $2=left primer position(index=0), $3=left primer sequence, $4=right primer position(index=0), $5=right primer sequence, $6=product length, $7=sequence length
awk < $outfile4 -F',' '$5 >0 {print $1","$5","$12","$14","$21","$25","$23}' > $outfile

exit;
