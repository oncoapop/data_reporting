#!/bin/sh

#wkfile="/home/dyap/Projects/PrimerDesign/manual/Splice-sign2_primer3_output.txt"
wkfile="/home/dyap/Projects/PrimerDesign/Tumour_Xenograft/primer3/SA533_primer3_output.txt"

echo "Processing..."
        filelen=` cat $wkfile | wc -l`

for i in `grep "SEQUENCE_ID=" $wkfile`
        do
        echo $i
        n=`grep -A$filelen $i $wkfile | grep -m1 "PRIMER_PAIR_NUM_RETURNED=" | awk -F"=" '{print $2}'`
                if [[ "$n" =~ "ok 0" ]];
                        then continue
                fi

        first=0
        last=`echo "$n - 1" | bc`
        echo $first, $last

        reclen=`grep -A$filelen $i $wkfile | gawk '{ RS = "\n=\n" ; FS = "=" }; END { print $0 }' | wc -l`
        #reclen=`sed -n '0,/^=/p' $wkfile | wc -l`
        echo $reclen

        for j in  $(eval echo "{$first..$last}")
                do
        left=`grep -A$reclen $i $wkfile | grep -m1 "PRIMER_LEFT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
        right=`grep -A$reclen $i $wkfile | grep -m1 "PRIMER_RIGHT_"$j"_SEQUENCE" | awk -F"=" '{print $2}'`
        size=`grep -A$reclen $i $wkfile | grep -m1 "PRIMER_PAIR_"$j"_PRODUCT_SIZE" | awk -F"=" '{print $2}'`
        id=`echo $i | awk -F"=" '{print $2}'`
        snv=`echo $id | awk -F"_" '{print $NF}'`
        chr=`echo $id | awk -F"_" '{print $(NF-1)}'`
                if [ -z "$right" -o -z "$left"  ];
                        then continue
                fi


                echo $id","$chr","$snv","$left","$right","$size
                done

	unset reclen

        echo -ne '#\r'
        done

