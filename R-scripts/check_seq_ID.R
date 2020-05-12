# This script was modified from biostars.org
# https://www.biostars.org/p/22/

# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome");
# biocLite("mygene")

library("mygene")

xli  <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
queryMany(xli, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="human")


