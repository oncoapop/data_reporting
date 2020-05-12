#!/bin/sh
# Script to extract the interactors from databases
# Modified for CLK1 interactors

pwd="/home/dyap/Projects/Takeda_T3/MS data/Expt3_CLK1"
cd "$pwd"

# Uniprot (2 interactions)
# downloaded from http://www.uniprot.org/uniprot/P49759#interaction 13 Jun 2016 
echo "Extracting Uniprot interactors..."
UFile="Uniprot_CLK1_interactors.txt"
cat $UFile | tail -n +2 | awk -F"\t" '{print $1}' > Uniprot_Interactors.txt

# BioGrid (168 unique interactions)
# downloaded from http://thebiogrid.org/107606 13 Jun 2016
echo "Extracting BioGrid interactors..."
Bfile="BIOGRID-GENE-107606-3.4.137.tab2.txt"
cat $Bfile | tail -n +2 | awk -F"\t" '{print $8"\n"$9}' | sort -u > BioGrid_Interactors.txt

# IntAct (20 unique interactions) 
# downloaded http://www.ebi.ac.uk/intact/pages/interactions/interactions.xhtml?query=P49759* 14 Jun 2016
echo "Extracting IntAct interactors..."
IAfile="P49759-.txt"
cat $IAfile | awk -F"\t" '{print $6}' | awk -F"(display_short)" '{print $1}' | awk -F":" '{print $2}' | sed 's/(.*//' | sed 's/_human//' | tr [a-z] [A-Z] | sort -u > IntAct_Interactors.txt

echo "Summarizing all interactors..."

cat Uniprot_Interactors.txt > all_UBIdb_interactors.tmp
cat BioGrid_Interactors.txt >> all_UBIdb_interactors.tmp
cat IntAct_Interactors.txt >> all_UBIdb_interactors.tmp


cat all_UBIdb_interactors.tmp | sort -u | wc -l
cat all_UBIdb_interactors.tmp | sort -u | sed 's/^Official*.$//g' | sed 's/^$//' > all_UBIdb_interactors.txt





