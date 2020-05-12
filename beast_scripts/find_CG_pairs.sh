#!/bin/sh

# Script to find matches in validated CG databases

########################
# References:
#Akiva, P., Toporik, A., Edelheit, S., Peretz, Y., Diber, A., Shemesh, R., Novik, A., & Sorek, R. (2006) Transcription-mediated gene fusion in the human genome. Genome Res. 16: 30-36.
#Denoeud, F., Kapranov, P., Ucla, C., Frankish, A., Castelo, R., Drenkow, J., Lagarde, J., Alioto, T., Manzano, C., Chrast, J. et al. (2007) 
#Prominent use of distal 5' transcription start sites and discovery of a large number of additional exons in ENCODE regions.Genome Res. 17: 746-759.
#Kim, N., Shin, S., & Lee, S. (2005a) ECgene: genome-based EST clustering and gene modeling for alternative splicing. Genome Res. 15: 566-576.
#Kim, P., Kim, N., Lee, Y., Kim, B., Shin, Y., & Lee, S. (2005b) ECgene: genome annotation for alternative splicing. Nucleic Acids Res. 33: D75-D79.
#Kim, N., Kim, P., Nam, S., Shin, S., & Lee, S. (2006) ChimerDB--a knowledgebase for fusion sequences. Nucleic Acids Res. 34: D21-D24.
#Parra, G., Reymond, A., Dabbouseh, N., Dermitzakis, E. T., Castelo, R., Thomson, T. M., Antonarakis, S. E., & Guigo, R. (2006) Tandem chimerism as a means to increase protein complexity in the human genome. Genome Res. 16: 37-44.

########################
wd="/home/dyap/Projects/Takeda_T3/CG"
cd $wd
# downloaded from http://metasystems.riken.jp/conjoing/download/ 8 Jan 2016

db1="NCBI_CG_34.txt"
db2="new-CG.txt"
db3="ChimerDB_CGDB.txt"
db4="Akiva_CGDB.txt"
db5="CG_Tissues.txt"
db6="ENCODE_CG.txt"
db7="ParentGeneDist.txt"
db8="Parra_CGDB.txt"
db9="ASA_validated.txt"
db10="Antibody_database.txt"

#############################
#T3_genes="T3_CG.txt"
#genes="all_genes.txt"
#output="matches-all.txt"

CGpanel="miseq_validation_readthrough_filtered_ctrl_norm.tsv"
genes="panel_genes.txt"
output="matches-52.txt"
###############################

echo $db1 > db1.$output
cat $db1 | wc -l >> db1.$output
echo "==============================" >> db1.$output
head -n1 $db1 >> db1.$output

echo $db2 > db2.$output
cat $db2 | wc -l  >> db2.$output
echo "==============================" >> db2.$output
head -n1 $db2 >> db2.$output

echo $db3  > db3.$output
cat $db3 | wc -l >> db3.$output
echo "==============================" >> db3.$output
head -n1 $db3 >> db3.$output

echo $db4 > db4.$output
cat $db4 | wc -l >> db4.$output
echo "==============================" >> db4.$output
head -n1 $db4 >> db4.$output

echo $db5 > db5.$output
cat $db5 | wc -l >> db5.$output
echo "==============================" >> db5.$output
head -n1 $db5 >> db5.$output

echo $db6 > db6.$output
cat $db6 | wc -l >> db6.$output
echo "==============================" >> db6.$output
head -n1 $db6 >> db6.$output

echo $db7 > db7.$output
cat $db7 | wc -l >> db7.$output
echo "==============================" >> db7.$output
head -n1 $db7 >> db7.$output

echo $db8 > db8.$output
cat $db8 | wc -l >> db8.$output
echo "==============================" >> db8.$output
head -n1 $db8 >> db8.$output

echo $db9 > db9.$output
cat $db8 | wc -l >> db9.$output
echo "==============================" >> db9.$output
head -n1 $db9 >> db9.$output

echo $db10 > db10.$output
cat $db10 | wc -l >> db10.$output
echo "==============================" >> db10.$output
head -n1 $db10 >> db10.$output

############################################
#cat $T3_genes | awk -F"@" '{print $1}' | sort -u > $genes
#cat $T3_genes | awk -F":" '{print $2}' | awk -F"@" '{print $1}' | sort -u >> $genes

cat $CGpanel | awk -F"\t" '{print $4}' | sort -u > $genes
cat $CGpanel | awk -F"\t" '{print $5}' | sort -u >> $genes
###########################################

	for i in `cat $genes`
	do
	echo "================"
	grep $i $CGpanel | awk -F"\t" '{print $4"-"$5}' | sort -u
#	grep $i $T3_genes | sort -u
	rdb1=`grep -w "$i" $db1`
	rdb2=`grep -w "$i" $db2`
	rdb3=`grep -w "$i" $db3`
	rdb4=`grep -w "$i" $db4`
	rdb5=`grep -w "$i" $db5`
	rdb6=`grep -w "$i" $db6`
	rdb7=`grep -w "$i" $db7`
	rdb8=`grep -w "$i" $db8`
	rdb9=`grep "$i" $db9`
	rdb10=`grep -w "$i" $db10`
	echo "================"

	if [[ $rdb1 != "" ]]
		then
		echo $rdb1 >> db1.$output
	fi
	if [[ $rdb2 != "" ]]
		then
		echo $rdb2 >> db2.$output
	fi
	if [[ $rdb3 != "" ]]
		then
		echo $rdb3 >> db3.$output
	fi
	if [[ $rdb4 != "" ]]
		then
		echo $rdb4 >> db4.$output
	fi
	if [[ $rdb5 != "" ]]
		then
		echo $rdb5 >> db5.$output
	fi
	if [[ $rdb6 != "" ]]
		then
		echo $rdb6 >> db6.$output
	fi
	if [[ $rdb7 != "" ]]
		then
		echo $rdb7 >> db7.$output
	fi
	if [[ $rdb8 != "" ]]
		then
		echo $rdb8 >> db8.$output
	fi
	if [[ $rdb9 != "" ]]
		then
		echo $rdb9 >> db9.$output
	fi
	if [[ $rdb10 != "" ]]
		then
		echo $rdb10 >> db10.$output
	fi

	done


echo "==============================" >> db1.$output
echo "==============================" >> db2.$output
echo "==============================" >> db3.$output
echo "==============================" >> db4.$output
echo "==============================" >> db5.$output
echo "==============================" >> db6.$output
echo "==============================" >> db7.$output
echo "==============================" >> db8.$output
echo "==============================" >> db9.$output
echo "==============================" >> db10.$output
