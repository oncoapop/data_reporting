#!/bin/bash

# Script to test inputs to a file

if [ "$1" = "" ] || [ "$2" = "" ] || [[ $1 != *[0-9]"-"[0-9]* ]] || [[ $2 != *[0-9]"-"[0-9]* ]]; then
        echo "usage:  "$0" <optimal size range> <secondary size range>"
        exit 0
fi

# Optimal size range (start-end)
size1=$1

# Seconday size  range (start-end)
size2=$2

echo "Opt size= " $size1
echo "Sec size= " $size2

exit

