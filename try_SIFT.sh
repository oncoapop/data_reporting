#!/bin/bash

myperl="/home/dyap/perl-5.14.4/perl"

$myperl ~/SIFTINDEL/bin/SIFT_exome_indels.pl -i ~/test/input37.txt -c hs37 \
	-d /home/dyap/SIFT_INDEL_HG37 -m 1 -o ~/ind231/

