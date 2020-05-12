#!/bin/sh

# This Script was written by Damian Yap (Apr 2013)
# To filter and add specific adaptors to primer3 generated primers for MiSeq
# Generate the order file and filtered AmpliconManifest

# Manually select the positions in which have been displayed
# at https://pleione.myseqtools.com/output/eIF4A3/eIF4A3_NMD/eIF4A3_NMD_summary.html

# This script is modified to use the primer supplementary file as the main identifier of
# selected transcripts

##########################################################################
name="eIF4A3_NMD"		# Sample set naming
type="Ill" 			# Illumina adaptors
#type="Fld" 			# Fluidigm adaptors
#type="none"			# NO ADAPTORS
Project="eIF4A3"		# Project name
RefGen="C:\Illumina\Miseq Reporter\Genomes\v19_transcriptome"	# Reference Genome for AmpliconManifest
ver="2"				# AmpliconManifest Version
##########################################################################

htmloutpath="/home/dyap/public_html"
dir=$htmloutpath"/"$Project"/"$name
p3dir="/home/dyap/Projects/PrimerDesign/"$Project"/primer3"
posdir="/home/dyap/Projects/PrimerDesign/"$Project"/positions"

##########################################
# Selected positions should be placed here
sourcefile=$posdir"/"$name"_selected_positions.txt"
# copy $p3dir"/short-listed_positions.txt" to $sourcefile (it is easier to delete the ones not
# wanted then to copy form screen to screen the ones you want to keep)
# Fomat 
# ENST00000207636:525+682 (copied from the HTML output)

cd $dir

suppfile=$p3dir"/"$name"_SupplFig.csv"
newsuppfile=$p3dir"/"$name"_SupplementalTable.csv"
orderleft=$p3dir"/"$name"-"$type"-leftprimers.csv"
orderright=$p3dir"/"$name"-"$type"-rightprimers.csv"
ordertube=$p3dir"/"$name"-"$type"-comb.csv"

						echo "WellPosition,Name,Sequence,Notes" > $orderleft
						echo "WellPosition,Name,Sequence,Notes" > $orderright
						echo "No,Name,Sequence,Notes" > $ordertube
						cat $suppfile | head -n 1 > $newsuppfile
	
oldmanifest=$p3dir"/"$name".AmpliconManifest"
newmanifest=$p3dir"/"$name"-filtered.AmpliconManifest"

rm -fr $newmanifest
rm -fr $newmanifest.tmp

	if [[ $type == "Ill" ]]
		then
		# Illlumina Adaptors (5'->3')
		fa="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
		ra="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
	fi

	if [[ $type == "Fld" ]]
		then	
		# Forward adaptor for Fluidigm
		fa="ACACTGACGACATGGTTCTACA"
		# Reverse adaptor for Fluidigm (5'->3')
		ra="TACGGTAGCAGAGACTTGGTCT"
	fi

	if [[ $type == "none" ]]
		then
		# NO ADAPTORS
		fa=""
		ra=""
	fi


# If file has a header (manual input)
# cat $sourcefile | tail -n +2 | awk -F" " '{print $1}' > pattern
# IF file is generated by a script (no header)
cat $sourcefile | awk -F" " '{print $1}' | sort -u > pattern

# Counters
alpha="A"
count=1
plate=1
# The pattern is what you use to uniquely identify the chose primer set after manual primer selection

		for i in `cat pattern`

		do
		ENST=`echo $i | awk -F":" '{print $1}'`
		start=`echo $i | awk -F":" '{print $2}' | sed 's/[+-]/,/' | awk -F, '{print $1}'`
		end=`echo $i | awk -F":" '{print $2}' | sed 's/[+-]/,/' | awk -F, '{print $2}'`

		echo $ENST,$start,$end

		label=`grep -w "$ENST" $suppfile | grep -w "$start" | grep -m1 -w "$end" | awk -F, '{print $1}'`
		left0=`grep -w "$ENST" $suppfile | grep -w "$start" | grep -m1 -w "$end" | awk -F, '{print $5}'`
		right0=`grep -w "$ENST" $suppfile | grep -w "$start" | grep -m1 -w "$end" | awk -F, '{print $6}'`

		echo "output:"
		echo $label,$left0,$right0

			if [[ $left0 != "" ]];
			
				then

				if [[ $right0 != "" ]];
			
					then

						if [[ $count -lt 10 ]];
							then 
								welpos=$alpha"0"$count
							else
								welpos=$alpha$count
						fi

						echo "Adding adaptors to sequences..."
						echo $welpos","$label"_L,"$fa$left0 >> $orderleft
						echo $welpos","$label"_R,"$ra$right0 >> $orderright
						echo $welpos","$label","$ra$right0","$fa$left0 >> $ordertube
						
						echo "Filtering AmpliconManifest..."
						grep -w "$ENST" $oldmanifest | grep -w "$start" | grep -m1 -w "$end" | awk -v label="$label" -F"\t" '{print label"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  >> $newmanifest.tmp
						grep -w "$ENST" $suppfile | grep -w "$start" | grep -m1 -w "$end" | awk -F"," '{print $1","$2","$3","$4","$5","$6","$7}'  >> $newsuppfile

				fi
			fi

			count=$((count+1))

			if [[ $count == 13 ]];
				then 
				count="1"
				if [[ $alpha == "I" ]]; 
					then 	alpha="A"
						plate=$((plate+1))
						orderleft=$p3dir"/"$name"-"$type"-leftprimers"$plate".csv"
						orderright=$p3dir"/"$name"-"$type"-rightprimers"$plate".csv"
				fi
				if [[ $alpha == "H" ]]; then alpha="I"; fi
				if [[ $alpha == "G" ]]; then alpha="H"; fi
				if [[ $alpha == "F" ]]; then alpha="G"; fi
				if [[ $alpha == "E" ]]; then alpha="F"; fi
				if [[ $alpha == "D" ]]; then alpha="E"; fi
				if [[ $alpha == "C" ]]; then alpha="D"; fi
				if [[ $alpha == "B" ]]; then alpha="C"; fi
				if [[ $alpha == "A" ]]; then alpha="B"; fi

			fi

		
		done

                echo "[Header]" > $newmanifest
                echo $name"  Manifest Version,"$ver | tr "," "\t" >> $newmanifest
                echo "ReferenceGenome,"$RefGen | tr "," "\t" >> $newmanifest
                echo "  " >> $newmanifest
                echo "[Regions]" >> $newmanifest
                echo "Name,Chromosome,Amplicon Start,Amplicon End,Upstream Probe Length,Downstream Probe Length,Comments" | tr "," "\t" >> $newmanifest
		cat $newmanifest.tmp | sort -u >> $newmanifest



		echo "Checking Amplicon Manifest for duplicates in the name column"

		dup=`awk -F"\t" '{print $1}' $newmanifest.tmp | sort | uniq -d `

			if [[ $dup != "" ]];
				then
				echo "DUPLICATES in AMPLICONMANIFEST"
				echo $dup
			fi

rm -f $p3dir/*.tmp

cp $newmanifest $dir"/".
cp $newsuppfile	$dir"/".

echo "Files can be found here:"
echo $htmloutpath"/"$Project"/"$sample

rsync -vr --progress $htmloutpath/$Project/$sample dyap@pleione.myseqtools.com:/var/www/html/output/$Project/
		
echo "done."
exit;