#!/bin/sh

# Script to upgrade the primer design pipeline
# Things that it does

# 1. copies makes a copy of the old pipeline into a new directory
# 2. Makes copies of the primer settings also
# 3. changes all the version-old to version-new (if you are using out-of-sync versions
#    they must be upgraded manually)

# QC pipeline
# oldver="v1.0_QC_"
# newver="v1.2_QC_"

# ChIP pipeline
# oldver="v1.2_ChIP_"
# newver="v1.3_ChIP_"

# Targeted Deep Seq pipeline
oldver="v3.3_"
newver="v3.4_"

olddir=$oldver"pipeline"
newdir=$newver"pipeline"

echo "Upgrading "$olddir" to "$newdir
cp -fr $olddir/ $newdir/

cd $newdir
pwd
echo "Before upgrading..."
ls

for file in `ls v*`
	do
	upfile=$(echo $file | sed "s/$oldver/$newver/")

	echo "Upgrading "$file" to "$upfile

	cat $file | sed "s/$oldver/$newver/g" > $upfile

	rm $file

	done

echo "Done."
pwd
ls

exit;

