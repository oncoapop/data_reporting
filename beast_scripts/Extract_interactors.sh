#!/bin/sh
# Script to extract the interactors from databases

pwd="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
cd "$pwd"

# Uniprot (26 interactions)
# downloaded from http://www.uniprot.org/uniprot/P49760#interaction 
echo "Extracting Uniprot interactors..."
UFile="Uniprot_CLK2_interactors.txt"
cat $UFile | tail -n +2 | awk -F"\t" '{print $1}' > Uniprot_Interactors.txt

# BioGrid (71 unique interactions)
# downloaded from http://thebiogrid.org/107607 4 Apr 2016
echo "Extracting BioGrid interactors..."
Bfile="BIOGRID-GENE-107607-3.4.135.tab2.txt"
cat $Bfile | awk -F"\t" '{print $8"\n"$9}' | sort -u > BioGrid_Interactors.txt

# IntAct (53 unique interactions) 
# downloaded http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=P49760* on 4 Apr 2016
echo "Extracting IntAct interactors..."
IAfile="P49760-.txt"
cat $IAfile | awk -F"\t" '{print $6}' | awk -F"(display_short)" '{print $1}' | awk -F":" '{print $2}' | sed 's/(.*//' | sed 's/_human//' | tr [a-z] [A-Z] | sort -u > IntAct_Interactors.txt

echo "Summarizing all interactors..."

cat Uniprot_Interactors.txt > all_UBIdb_interactors.tmp
cat BioGrid_Interactors.txt >> all_UBIdb_interactors.tmp
cat IntAct_Interactors.txt >> all_UBIdb_interactors.tmp


cat all_UBIdb_interactors.tmp | sort -u | wc -l
cat all_UBIdb_interactors.tmp | sort -u > all_UBIdb_interactors.txt





