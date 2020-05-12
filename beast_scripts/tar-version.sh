#!/bin/sh
####################################
#
# Tar version pipeline script.
#
####################################

# Pipeline version
# ver="4.2" # 1,2 Dec 14 for cell mix HCT116/hTert ss
# ver="3.3" # 18 Feb 15 for Splice Signature 2
# ver="4.3" # 2 Mar 15 for cell mix HCT116/hTert & DAH54/DAH55/Shared
# ver="1.1" # 11 Jun 15 for CX5461 ChIPseq analysis pipeline
# ver="1.0" # 23 Jun 15 for QPCR QC pipeline for targeted Deep Seq pipeline
# ver="3" # 18 Aug 15 for upgraded primer design pipeline for targeted Deep Seq pipeline (multiple amplicons per position)
# ver="3.3" # 23 oct 15 for Splice RT using isPCR on transciptome
# ver="3.3" # 15 Mar 16 for Splice RT using isPCR on transciptome (with new factors)
# ver="1.0" # 15 Nov 16 for QC of MiSeq Reporter DNA (Cluster plots of VAF)
# ver="1.2" # 25 Jan 17 for QC of primers orders by qPCR (correlation plots)
# ver="1.2" # 9,10 Feb 17 for QC of primers orders by qPCR (correlation plots)
ver="3.4" # 9,10 Feb 17 for QC of primers orders by qPCR (correlation plots)

# basedir
basedir="/home/dyap"

#pipeline="v"$ver"_VCF_pipeline"
pipeline="v"$ver"_pipeline"
#pipeline="v"$ver"_ChIP_pipeline"
#pipeline="v"$ver"_QC_pipeline"

# What to backup.
backup_files=$basedir"/Scripts/"$pipeline

# Where to backup to.
dest=$basedir"/Projects/PrimerDesign/backup"
remote="/meta/o/oncoapop/backups/beast"
remote2="/backup/beast"

# Create archive filename.
day=$(date +"%y_%m_%d")

archive_file=$pipeline"_"$day".tgz"

# Print start status message.
echo "Backing up $backup_files to $dest/$archive_file"
date
echo

# Backup the files using tar.
tar -czf $dest/$archive_file "$backup_files"

# Print end status message.
echo
echo "Backup finished"
date

# Long listing of files in $dest to check file sizes.
ls -lh $dest
scp $dest"/"$archive_file oncoapop@ma.sdf.org:$remote"/"
scp $dest"/"$archive_file dyap@pleione.myseqdept.org:$remote2"/"

# Restoring
echo "==========================================================================="
echo "To Restore type tar -xzvf $dest/$archive_file from your $HOME"

