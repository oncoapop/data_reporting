# This R sript is to view and plot certain stats from the output of deFuse pipeline
# The output may change and hence manual inspection is required for confirmation

# The annotation of the output and descriptors can be found here
# https://bitbucket.org/dranew/defuse#markdown-header-output

#install.packages("ggplot2")
#install.packages("plotly")
library(ggplot2)
library(car)

# Directory where Tyler outputs his scripts
#wd="/share/lustre/projects/takeda_splicing_inhibitor/gene_fusions/data/miseq_validation"
wd="/share/lustre/projects/takeda_splicing_inhibitor/gene_fusions/data/SRRM_KD"
outdir="/home/dyap/Projects/Takeda_T3/CG"

#file="miseq_validation_readthrough_filtered.tsv"
file="SRRM_KD_readthrough_filtered.tsv"

input=paste(wd,file,sep="/")

sum<-read.table(file=input, sep="\t", header=TRUE)

pdffile=paste(outdir,"SplitRCount.pdf",sep="/")
pdf(file=pdffile)
#boxplot(splitr_count~sample_id,data=sum, main="Split read Counts by sample (MiSeq run A92U5)", 
#  	xlab="Sample ID", las=2, ylab="No of reads", cex.axis=0.8)
boxplot(splitr_count~sample_id,data=sum, main="Split read Counts by sample (MiSeq run AELDU)", 
  	xlab="Sample ID", las=2, ylab="No of reads", cex.axis=0.8)
dev.off()

boxplot(probability~sample_id,data=sum, main="Probability by sample (MiSeq run A92U5)", 
  	xlab="Sample ID", las=2, ylab="Probability", cex.axis=0.8,labels,id.n=Inf,
	col = "lightgray")
boxplot(probability~sample_id,data=sum, main="Probability by sample (MiSeq run AELDU)", 
  	xlab="Sample ID", las=2, ylab="Probability", cex.axis=0.8,labels,id.n=Inf,
	col = "lightgray")

pdffile=paste(outdir,"SplitRCount_lowrange.pdf",sep="/")
pdf(file=pdffile)
stripchart(splitr_count~sample_id,data=sum,vertical = TRUE, pch=1, method="jitter",ylim=c(0,3000))

dev.off()

# subsetting dataset
# 
