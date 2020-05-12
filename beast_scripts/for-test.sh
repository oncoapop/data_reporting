#!/bin/sh

projdir="/home/dyap/Projects/Takeda_T3/CG"

cd $projdir
rm -f *.tmp*

samples="SA464 SA465 SA466 SA467 SA468 SA469 SA470 SA502 SA503 SA504 SA505 SA537 SA538 SA539 SA540"

	for j in  $(eval echo "$samples")
		do
		echo $j
		var_name=${j}
        	samplename=`echo ${var_name}`

		echo $var_name, $samplename

		done
exit;
