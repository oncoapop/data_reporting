#! /bin/bash
histo="======================================================================+"

read datewd value

while [ -n "$datewd" ] ; do
   # Use a default width of 70 for the histogram
   echo -n "$datewd      "
   echo ${histo:0:$value}

   read datewd value
done
