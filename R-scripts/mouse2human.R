# This is an R script that uses BioMart to get the MU <-> MU of ENSG
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# import data
wd="/home/dyap/Projects/eIF4A3_NMD/comparisons"
infile="Mao_2016_Eif4a3.txt"

file=paste(wd,infile,sep="/")

# Read from the input file the list of ENS IDs 
input=read.csv(file=file, StringAsFactors=FALSE)

# Input either the human or the mouse ENS IDs

query_ids<-as.character(input$GeneID)

# Use the BioMart for Human and Mouse comparisons
hu_mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mu_mart = useMart("ensembl", dataset="mmusculus_gene_ensembl") 


# human / mouse
output_ids <- getLDS(attributes=c("ensembl_gene_id"),
           filters="ensembl_gene_id", values=query_ids, mart=mu_mart,
           attributesL=c("ensembl_gene_id"), martL=hu_mart)
