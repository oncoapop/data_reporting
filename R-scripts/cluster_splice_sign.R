# use myR R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"

# Script to visualize data
#  source("http://www.bioconductor.org/biocLite.R"); biocLite("Heatplus")
# source("http://www.bioconductor.org/biocLite.R"); biocLite("limma")
# install.packages("gplots")

library(foreign)
library(RColorBrewer)
library(lattice)
library(colorspace)
require(Heatplus)
# library(limma)
library(gplots)

rgb.palette <- colorRampPalette(c("dark green", "red"), space = "rgb")
hmcols <- colorRampPalette(brewer.pal(11,"Spectral"))(1000)
# display.brewer.all()

# added extra col at the end to make it import correctly
counts<-read.table(file="/share/lustre/backup/dyap/Projects/Takeda_SpliceSignature/RocheResults/counts/outputCounts", sep="\t", header=TRUE, strip.white = TRUE, row.names=1, 
stringsAsFactors=F)
ratios<-read.table(file="/share/lustre/backup/dyap/Projects/Takeda_SpliceSignature/RocheResults/ratios/outputRatios", sep="\t", header=TRUE, strip.white = TRUE, row.names=1, 
stringsAsFactors=F)
#ratios<-read.table(file="/share/lustre/backup/dyap/Projects/Takeda_SpliceSignature/RocheResults/ratios/CLK_1_3_outputRatios", sep="\t", header=TRUE, strip.white = TRUE, row.names=1, 
#stringsAsFactors=F)

outdir="/share/lustre/backup/dyap/Projects/Takeda_SpliceSignature/RocheResults/ratios"

rows<-nrow(ratios)
cols<-ncol(ratios)-1

mat1<-as.matrix(ratios[1:cols])
dat1<-as.data.frame(ratios[1:cols])

# subset by cell line or dose
#toMatch <- c("0uM", "Trtd", "Ctrl")
toMatch <- c("HCT116")
mat<-as.matrix(mat1[ , grep(paste(toMatch,collapse="|"), colnames(mat1), value=TRUE)])

# Colside color matrix

drugdose <- data.frame( Sample = rep("", ncol(mat)),
                        Col = rep("", ncol(mat)),
                        stringsAsFactors = FALSE)
count=1
# select colours
# display.brewer.all()

for (a in seq(ncol(mat)))
        {
        drugdose$Sample[a] <- as.character(colnames(mat)[a])

        if (length(grep("Ctrl", drugdose$Sample[a])) > 0)
                {
                count=1
                drugdose$Col[a] <- "#FFFFFF"
                }
        if (length(grep("_0uM", drugdose$Sample[a])) > 0)
                {
                count=1
                drugdose$Col[a] <- "#FFFFFF"
                }
        if (length(grep("CDK12", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"YlOrRd")[count]
                }
        if (length(grep("EIF4A3", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"Oranges")[count]
                }
        if (length(grep("DYRK_.1.2.", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"Greys")[count]
                }
        if (length(grep("SRPK1", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"Purples")[count]
                }
        if (length(grep("CLK_.1.4.", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"Blues")[count]
                }
        if (length(grep("CLK_.1.3.", drugdose$Sample[a])) > 0)
                {
                drugdose$Col[a] <- brewer.pal(9,"Greens")[count]
                }

        count=count+2

        }
palette<-brewer.pal(5,"Paired")

lines <- data.frame( Sample = rep("", ncol(mat)),
                        Col = rep("", ncol(mat)),
                        stringsAsFactors = FALSE)

for (a in seq(ncol(mat)))
        {
        lines$Sample[a] <- as.character(colnames(mat)[a])
        if (length(grep("HeLa", lines$Sample[a])) > 0)
                {
                lines$Col[a] <- palette[1]
                }
        if (length(grep("SK.BR.3", lines$Sample[a])) > 0)
                {
                lines$Col[a] <- palette[2]
                }
        if (length(grep("hTert", lines$Sample[a])) > 0)
                {
                lines$Col[a] <- palette[3]
                }
        if (length(grep("HCT116", lines$Sample[a])) > 0)
                {
                lines$Col[a] <- palette[4]
                }
        }

legend=c("HeLa","SK.BR.3","hTert","HCT116")
fill=palette[1:4]

title="Splice Ratios from RNA-seq data (10uM dose)"
xaxislab2="samples"
pdffile=paste(outdir,"HCT116_Overview.pdf",sep="/")

# Sort by ratios
pdf(pdffile, width=8, height=10)

{
heatmap.2(mat, main=title, xlab=xaxislab2, ylab="Positions", scale="none", Colv=TRUE, key = TRUE
, cexCol=0.8, srtCol=45, cexRow=0.8, col = hmcols, trace="none",margins=c(8,6), ColSideColors=drugdose$Col,dendrogram="none")

#legend("topright",legend=legend, fill=fill, border=TRUE, bty="o", y.intersp = 0.7, cex=1.0)
}

dev.off()

# Sort by cell lines and drugdose
pdf(pdffile, width=8, height=10)

{
heatmap.2(mat, main=title, xlab=xaxislab2, ylab="Positions", scale="none", Colv=TRUE, key = TRUE
, cexCol=0.8, srtCol=45, cexRow=0.8, col = hmcols, trace="none",margins=c(8,6), ColSideColors=lines$Col, colsep=c(1,4,5,6,11,15,16,18,22,26,30,31,32,37,40,41,43,46,48), 
dendrogram="col")

legend("topright",legend=legend, fill=fill, border=TRUE, bty="o", y.intersp = 0.7, cex=1.0)
}

dev.off()

# subset by raw and pulled 
raw<-seq(1,91,2)
pulled<-seq(2,92,2)
nat<-mat[,raw]
pul<-mat[,pulled]

# subset by cell line
HCTr<-as.matrix(nat[ , grep("HCT116", colnames(nat))])
HCTp<-as.matrix(pul[ , grep("HCT116", colnames(pul))])
#mat<-HCTp
drugdose <- data.frame(	Sample = rep("", ncol(mat)),
			Col = rep("", ncol(mat)),
			stringsAsFactors = FALSE)
count=1

for (a in seq(ncol(mat)))
	{
	drugdose$Sample[a] <- as.character(colnames(mat)[a])
	if (length(grep("_0uM", drugdose$Sample[a])) > 0) 
		{
		count=1
		drugdose$Col[a] <- "#FFFFFF"
		count=2
		next
		}
	if (length(grep("_Ctrl", drugdose$Sample[a])) > 0) 
		{
		count=1
		drugdose$Col[a] <- "#FFFFFF"
		count=2
		next
		}
	if (length(grep("EIF4A3", drugdose$Sample[a])) > 0) 
		{
		drugdose$Col[a] <- brewer.pal(9,"Oranges")[count]
		}
	if (length(grep("DYRK_.1.2.", drugdose$Sample[a])) > 0) 
		{
		drugdose$Col[a] <- brewer.pal(9,"Greys")[count]
		}
	if (length(grep("SRPK1", drugdose$Sample[a])) > 0) 
		{
		drugdose$Col[a] <- brewer.pal(9,"Purples")[count]
		}
	if (length(grep("CLK_.1.4.", drugdose$Sample[a])) > 0) 
		{
		drugdose$Col[a] <- brewer.pal(9,"Blues")[count]
		}
	if (length(grep("CLK_.1.3.", drugdose$Sample[a])) > 0) 
		{
		drugdose$Col[a] <- brewer.pal(9,"Greens")[count]
		}
	count=count+1
	print(count)
	}

# Colside color matrix
csc <- data.frame(	Sam = rep("", ncol(mat)),
			Col = rep("", ncol(mat)),
			stringsAsFactors = FALSE)

# Get the colors based on different sample proportions
# This is redundant but left in 
for (j in seq(ncol(mat)))
	{
	match <- colnames(mat)[j]
	csc$Sam[j]<- match
	test <- drugdose[drugdose[,1] %in% c(match),2]
	if (length(test) != 0) 
		{
		csc$Col[j] <- test
		} 
		
	}

title="Splice Ratios from RNA-seq data"
xaxislab2="samples"
pdffile=paste(outdir,"HCT116_drug.pdf",sep="/")


# Sort by ratios
pdf(pdffile, width=8, height=10)
{
par(mar=c(4,2,0,2)+0.1)
heatmap.2(mat, main=title, xlab=xaxislab2, ylab="Positions", scale="none", Colv=TRUE, key = TRUE
, cexCol=0.6, srtCol=45, cexRow=0.5, col = hmcols, trace="none",margins=c(8,6), ColSideColors=drugdose$Col,dendrogram="col")

legend("bottomleft",legend=drugdose$Sample, fill=drugdose$Col, border=TRUE, bty="o", y.intersp = 0.7, cex=1.0)

#colsep=c(4, 8, 12, 16, 24))

#text(x = c(3,5,7,9,11), y = 27, 
#     labels = c("EIF4A3", "DYRK","SRPK1","CLK_1,4","CLK_1,3"),
#     las = 2, col = "black", cex = 0.8, xpd = TRUE)
}
 
dev.off()

full <- heatmap.2(x, trace="none",Colv=FALSE)
nbk=length(full$breaks)-1
hmcols <- colorRampPalette(brewer.pal(11,"Spectral"))(nbk)

###########

pdffile=paste(outdir,"raw.pdf",sep="/")

# Sort by ratios
pdf(pdffile, width=10, height=13)

par(mar = c(12, 3, 8, 10))
title="Splice Ratios from raw RNA-seq"
heatmap.2(x, Colv=full$colDendrogram[[1]], breaks=full$breaks, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.5, cexRow=0.4, col = hmcols, trace="none",srtCol=45)  # raw column subset

dev.off()
###########
pdffile=paste(outdir,"captured.pdf",sep="/")

# Sort by ratios
pdf(pdffile, width=10, height=13)

title="Ratios: in silico capture (RNA-seq)"
par(mar = c(12, 3, 8, 10))
heatmap.2(x, Colv=full$colDendrogram[[2]], breaks=full$breaks, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.5, cexRow=0.4, col = hmcols, trace="none", srtCol=45)  # in-sillico RNA capture column subset

dev.off()

######

amplicons <- data.frame ( Sam= rep("", nrow(primers)),
                         Counts = rep(0, nrow(primers)),
                     stringsAsFactors = FALSE)

for (ri in seq(nrow(counts)))
        {
        amplicons$Sam[ri] <- as.character(primers[ri,1])
        amplicons$Counts[ri] <- as.numeric(gsub("bp", "", strsplit(as.character(primers[ri,5]),
split=":")[[1]][1]))
        }

samples<-colnames(counts)[2:ncol(counts)]
amplicons<-rownames(counts)


#heatmap.2(counts, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
#, cexCol=0.8, cexRow=0.6, col = hmcols, RowSideColors=rsc$Col, ColSideColors=csc$Col, trace="none")


###################

# Sort by mutations
pdf("By_Mutations.pdf", width=10, height=10)

par(mar = c(5, 6, 4, 2))
barplot(sort(table(gene$SampleID),decreasing = TRUE), horiz = TRUE, las = 1,xlab = "Number of Mutations", main = "Mutations per 
Sample")

dev.off()

# Sort by Samples
pdf("By_Samples.pdf", width=10, height=10)

par(mar = c(5, 6, 4, 2))
barplot(sort(table(gene$GENE_NAME),decreasing = TRUE), horiz = TRUE, las = 1,xlab = "Number of Individual Mutations", main = 
"Mutations by Gene (some >1 per sample)")


dev.off()

# two non numeric var
pdf("Combined_anno.pdf", width=10, height=10)

par(mar = c(5, 6, 4, 2))
#genecol<-brewer.pal(n=15, "Spectral")
genecol<-rainbow_hcl(15, start = 1, end = 300)

barplot(table(gene$GENE_NAME,gene$SampleID), horiz = TRUE, las = 1, xlab = "Number of Individual Mutations", main = "Mutations by 
Gene (some >1 per sample)", col=genecol)

legend('topright', unique(sort(gene$GENE_NAME)), lty=1, lwd=10, col=genecol, bty='n', cex=.75)

dev.off()

sum1<-table(gene$GENE_NAME,gene$SampleID)

write.table(sum1,file="/home/dyap/Projects/FFPE/Table.csv",sep=",",row.names=TRUE,col.names=TRUE)

