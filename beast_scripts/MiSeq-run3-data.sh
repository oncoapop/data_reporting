#/bin/sh

# script to grep identities from blast aligned files

# written by Damian Yap (Apr 2013)
# written on Ubantu 12, now on Redhat Cluster

# Source directory
indir='/home/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq/Blastfasta' 

# Destination directory
outdir='/home/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq/Blastfasta/summary' 

#position file
posfile=/home/dyap/Projects/Single_Cell/positions/VBA0038_all48.txt

cd $outdir
rm *.check.txt
rm *.sum.txt

cd $indir

	for i in `ls Nuc*-*read.txt`

	do

	echo $i
	name=`echo $i | sed 's/-read.txt//'`
	infile=$i


		for j in {1..48}

		do
				
		WT=`grep -w $j $posfile | awk -F, '{print $5,$4,$6}' | tr -d " "`
		A=`grep -w $j $posfile | awk -F, '{print $5,"A",$6}' | tr -d " "`
		T=`grep -w $j $posfile | awk -F, '{print $5,"T",$6}' | tr -d " "`
		G=`grep -w $j $posfile | awk -F, '{print $5,"G",$6}' | tr -d " "`
		C=`grep -w $j $posfile | awk -F, '{print $5,"C",$6}' | tr -d " "`

		pos=`grep -w $j $posfile | awk -F, '{print $2,$3}' | tr " " ":"`

		countwt=`grep -c $WT $infile`
		counta=`grep -c $A $infile`
		countt=`grep -c $T $infile`
		countg=`grep -c $G $infile`
		countc=`grep -c $C $infile`

		echo $pos,$countwt,$counta,$countt,$countg,$countc >> $outdir/$infile.sum.txt
		echo Pos=$pos,WT=$WT,A=$A,T=$T,G=$G,C=$C >> $outdir/$infile.check.txt

		done

	echo ++++++++++++++++++++++++++++++++++++++++++++++++++

	done


echo All done.

exit;

