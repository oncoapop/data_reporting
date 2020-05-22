#!/bin/bash

wd="/home/dyap/ind231/SIFT_RESULTS"
cd $wd
outdir="$HOME/ind231/correct_errors"
correct="$HOME/ind231/corrected_errors.csv"

cat *.chrfile.err | awk -F"[" '{print $2}' | awk -F"]" '{print $1}' | tr -d "+" > $correct

cd $outdir
split -d -l 10 $correct correct

