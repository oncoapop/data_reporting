#source("http://bioconductor.org/biocLite.R")
#biocLite("Sushi")

# Modified for the original script from Sam Aparicio (Jun 2015)
# Modified by Damian Yap 9 Jun 2015 and 24 Jul 2015 for RNA-seq data

pwd<-"/home/dyap/Projects/Takeda_SpliceSignature/Sign2_HCT116_May15/Roche_150725_HG19_Splice-Sign2_RNA"
setwd(pwd)

fname="Splice-sign2"

########
# CODE
########

library("Sushi")
pdfname=paste(fname,"pdf",sep=".")
makepdf = TRUE


#FUNCTION FOR READING IN BED FILES VERY QUICKLY
# example of bed file
#chr1    15949184        15949402        -_12_CASP9_SE_exon1;-_13_CASP9_SE_exon1;-_14_CASP9_SE_exon1;-_15_CASP9_SE_exon1
#chr1    15958584        15958707        -_12_CASP9_SE_exon2

read.bed <- function(file) {
 dat <- scan(file=file,
         what=list(character(),integer(),integer(),character()),
         sep="\t", skip=0)
 dat <- data.frame(chr=dat[[1]], chromstart=dat[[2]],
         chromend=dat[[3]], gene=dat[[4]])
 return(dat)
}

read.bed2 <- function(file) {
 dat2 <- scan(file=file,
         what=list(character(),integer(),integer(),character(),integer(),character(),character()),
         sep="\t", skip=0)
 dat2 <- data.frame(chrom=dat2[[1]], chromstart=dat2[[2]],
         chromend=dat2[[3]], gene=dat2[[4]],score=dat2[[5]],strand=dat2[[6]],type=dat2[[7]]
	)
 return(dat2)
}

Cap_reg<-read.bed("Splice-Sign2_capture_targets.bed.txt")
Uncov_reg<-read.bed("Splice-Sign2_predicted_uncovered_targets.bed.txt")
Design_reg<-read.bed("Splice-Sign2_primary_targets.bed.txt")

#READ IN REGIONS -
exons.bed<-read.bed2("exon_regions.bed")

#READ IN REGIONS -
#gene.bed<-data.frame(read.table("gene_regions.bed",col.names=c("chrom","chromstart","chromend","gene","score","strand","type"),sep="\t"))
gene.bed<-read.bed2("gene_regions.bed")

############## OUTPUT ###################
if (makepdf == TRUE)
{
   pdf(pdfname,height=10, width=12)
}


chrom = "chr10"
chromstart = 90770296
chromend = 90771838
range=1
types="exon"
pg = plotGenes(Cap_reg,chrom,chromstart-range,chromend+range,
		colorbycol= SushiColors(4)(4)[1],colorbyrange=c(0,1.0),
		labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
labelgenome( chrom, chromstart,chromend,n=5,scale="bp")

############################## TESTING ################################


#GENERIC FUNCTION FOR PLOTTING FOUR LIBRARIES AT A TIME
bed.region.plot <- function(d1,d2,d3,d4,chrom,chromstart,chromend) 
{try(plotGenes(d1,chrom,chromstart-range,chromend+range,
		colorbycol= SushiColors(4)(4)[1],colorbyrange=c(0,1.0),
		labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box"))
try(plotGenes(d2,chrom,chromstart-range,chromend+range,
		colorbycol= SushiColors(4)(4)[2],colorbyrange=c(0,1.0),
		labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box",overlay=TRUE))
try(plotGenes(d3,chrom,chromstart-range,chromend+range,
		colorbycol= SushiColors(4)(4)[3],colorbyrange=c(0,1.0),
		labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box",overlay=TRUE))
try(plotGenes(d4,chrom,chromstart-range,chromend+range,
		colorbycol= SushiColors(4)(4)[4],colorbyrange=c(0,1.0),
		labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box",overlay=TRUE))
labelgenome( chrom, chromstart,chromend,n=4,scale="bp")
legend("topright",inset=0.025,legend=c("Exon","Capture region","Uncovered region","Design Input"),fill=opaque(SushiColors(4)(4)),border=SushiColors(4)(4),text.font=2,cex=1.0)
}

#FUNCTION FOR EXECUTING PLOTS OVER SEVERAL REGIONS
execute.plot<-function(gene.bed) {
for (n in 1:nrow(gene.bed)) {
bed.region.plot(exons.bed,Cap_reg,Uncov_reg,Design_reg,gene.bed$chrom[n],gene.bed$chromstart[n],gene.bed$chromend[n])
}
}

execute.plot(gene.bed)

if (makepdf == TRUE)
{
   dev.off
}




