#!/bin/sh

# Script to copy vcf into directory to be read
# input directory
#source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/"
source="/share/lustre/archive/MiSeq/MiSeq_Analysis_Files/161118_M02348_0131_000000000-AT094/Data/Intensities/BaseCalls/Alignment/VariantCallingLogs"

# output directory
dest="/home/dyap/Projects/ctDNA/AT094/variants"

cp $source/*[0-9].vcf $dest/
 
# Find somatic mutations
# vcf present on the FFPE but not in the Genome (G -WT, buffy coat)

cd $dest

# Get all the sample IDs

#for i in `ls *[0-9].vcf  | awk -F. '{print $1}' | grep F | sed 's/F//'`
for i in `seq 25 2 52`
	do
	echo $i
	pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $2}'`	

	echo $pos

	next=`echo "$i + 1" | bc`

	# common mutations
	grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_common"


	echo "###################################################"
	done

for i in `seq 26 2 52`
	do
	echo $i
	pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $2}'`	

	echo $pos

	next=`echo "$i - 1" | bc`

	# common mutations
	grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_common"


	echo "###################################################"
	done

### for Germline
for i in `seq 1 2 24`
	do
	echo $i
	pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $2}'`	

	echo $pos

	next=`echo "$i + 1" | bc`

	# common mutations
	grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_common"


	echo "###################################################"
	done

for i in `seq 2 2 24`
	do
	echo $i
	pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $2}'`	

	echo $pos

	next=`echo "$i - 1" | bc`

	# common mutations
	grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "Sample"$i"_common"


	echo "###################################################"
	done
	

#######################
# Common between pools
# first 5 columns

for i in `seq 1 2 24`
        do
        echo $i
        pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'`

        echo $pos

        next=`echo "$i + 1" | bc`

        # common mutations
        grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "Buffy"$i"_common"


        echo "###################################################"
        done

for i in `seq 24 2 48`
        do
        echo $i
        pos=`grep "^chr" *S$i.vcf | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'`

        echo $pos

        next=`echo "$i + 1" | bc`

        # common mutations
        grep "$pos" *S$next.vcf | awk -F"\t" '{print $0}' > "FFPE"$i"_common"


        echo "###################################################"
        done



