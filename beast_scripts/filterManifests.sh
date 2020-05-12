#!/bin/sh

# Script to check generated filtered AmpliconManfests with 
# AmpliconManifests generated using the primer3 display pipeline

basedir="/home/dyap/public_html"

batch="Tumour_Xenograft_Rev"
#set="set2"

sample="SA533"

wd=$basedir"/"$batch"/"$sample

cd $wd

longmanifest=$sample".AmpliconManifest.txt"
#shortmanifest=$sample"-set2.AmpliconManifest"
shortmanifest=$sample"-cal.AmpliconManifest"
manfile=$sample"-filtered.AmpliconManifest"


# Amplicon Manifest Header
                echo "[Header]" > $manfile
#                echo $sample"-"$set"  Manifest Version,1" | tr "," "\t" >> $manfile
                echo $sample"  Manifest Version,1" | tr "," "\t" >> $manfile
                echo "ReferenceGenome,C:\Illumina\MiSeq Reporter\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA" | tr "," "\t" >> $manfile
                echo "  " >> $manfile
                echo "[Regions]" >> $manfile
                echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" | tr "," "\t" >> $manfile

for i in `cat $shortmanifest | tail -n+7 | awk -F"\t" '{print $1}'`
	do
	id=`grep $i $longmanifest | awk -F"\t" '{print $1}'`
	chr=`grep $i $longmanifest | awk -F"\t" '{print $2}'`
	sta=`grep $i $longmanifest | awk -F"\t" '{print $3}'`
	end=`grep $i $longmanifest | awk -F"\t" '{print $4}'`
	llen=`grep $i $longmanifest | awk -F"\t" '{print $5}'`
	rlen=`grep $i $longmanifest | awk -F"\t" '{print $6}'`
	amplen=`grep $i $longmanifest | awk -F"-" '{print $3}' | sed 's/bp//'`

	calid=`grep $i $shortmanifest | awk -F"\t" '{print $1}'`
	calchr=`grep $i $shortmanifest | awk -F"\t" '{print $2}'`
	calsta=`grep $i $shortmanifest | awk -F"\t" '{print $3}'`
	calend=`grep $i $shortmanifest | awk -F"\t" '{print $4}'`
	calllen=`grep $i $shortmanifest | awk -F"\t" '{print $5}'`
	calrlen=`grep $i $shortmanifest | awk -F"\t" '{print $6}'`
	calamplen=`grep $i $longmanifest | awk -F"--" '{print $2}'`

	echo "--------------------------"
	echo $i

	if [ "$id" = "$calid" ] 
		then
			echo "ID match"
	fi

	if [ "$chr" = "$calchr" ] 
		then
			echo "Chr match"
	fi

	if [ "$sta" = "$calsta" ] 
		then
			echo "Start match"
	fi

	if [ "$end" = "$calend" ] 
		then
			echo "End match"
	fi

	if [ "$llen" = "$calllen" ] 
		then
			echo "Left Primer Length match"
	fi

	if [ "$rlen" = "$calrlen" ] 
		then	
			echo "Left Right Length match"
	fi

	if [ "$amplen" = "$calamplen" ] 
		then
			echo "Amplicon Length match"
	fi

	if [ "$id" = "$calid" ] 
		then
			grep $i $longmanifest >> $manfile
	fi

	
	done

echo "done!"

exit;
