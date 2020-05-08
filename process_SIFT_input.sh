#!/bin/bash

wd="/home/dyap/ind231"
cd $wd

file="expanded_list_grouped_V5b.csv"

cat $file | awk -F"," '{print $1","$3","$7"1,"$5"/"$6}'
