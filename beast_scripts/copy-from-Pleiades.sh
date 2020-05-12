#!/bin/bash
clear
# Script to sync the directories of Pleiades and Godel

# uses port 1953 for ssh connection to pleiades
# Syncs the contents of Tumour_Xenograft on Pleiades with v3_pipeline on Godel
# Pipeline (remote directory)
Pipeline="v3_pipeline"

echo "Rsync of files Pleiades to Godel"
echo "Directory on Pleiades : /home/dyap/Scripts/"$Pipeline
echo "Directory on Godel    : /home/dyap/Scripts/"
echo "---------------------------------------------------------------------------"
echo "Press ENTER to continue or Ctrl-C to exit..."
read ans

rsync -aevr --progress --inplace --rsh='ssh -p1953'  dyap@damianeva.dyndns.org:/home/dyap/Scripts/$Pipeline /home/dyap/Scripts/
