#!/bin/sh

#start
m=1

#end
l=23

for j in  $(eval echo "{$m..$l}")
      do
       var_name="file${j}"
       filename=`echo ${!var_name}`

	echo $var_name,$filename
	
      done
