#source("http://bioconductor.org/biocLite.R")
#iocLite()
#source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
# biocLite(c("GenomicFeatures", "AnnotationDbi"))


#library("VariantAnnotation")
#library("IRanges")
#library("GenomicRanges")
library(foreign)
library(lattice)
library(XLConnect)

#####

# To run this script change the setwd()
# on beast
setwd("/share/lustre/backup/dyap/Projects/Single_Cell/SA212/SA429_spike-in/")

getwd()

dflist <- readWorksheetFromFile(file="/share/lustre/backup/dyap/Projects/Single_Cell/SA212/SA429_spike-in/SA212_MS130605_ResultsSummary.xls", sheet = c("SA212_target_list_130605.txt.de", 
"SA212_target_list_130605.txt.ma","SA212_target_list_130605.txt.pv"), header = TRUE, startRow = c(1,1,1), startCol = c(1,1,1), endCol = c(48,48,48),endRow = c(64,64,64))


sum1 <- dflist$SA212_target_list_130605.txt.ma
sum2 <- dflist$SA212_target_list_130605.txt.de

#####################################

# Drawing heatmap

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])
ff <- data.matrix(sum2[2:ncol(sum2)])

# col headers - unique nuclei
names(sum1)

# Label rownames with ID (these are positions)
rownames(ef) <- sum1$target
rownames(ff) <- sum2$target


########################################################

# pdf("SA429-allfreq.pdf", width=10, height=4)
pdf("SA212-Read_BiNom.pdf", width=10, height=4)

###################
#heatmap(ef, Rowv=NA, Colv=NA, na.rm=TRUE, main="Selected Variant Freq SA494 (Run ID A49AR)", xlab="SA494 Nuclei", ylab="Position", cexCol=0.8, col=rev(heat.colors(1000)))
# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", margins=c(5,10))
###################

rgb.palette <- colorRampPalette(c("dark green", "red"), space = "rgb")

# levelplot(ef, xlab="Position", las=2, ylab="Nuclei", main="BiNom MAF SA212 (Run ID A499H)", col.regions=rgb.palette(10), cuts=10, at=seq(0.2,1.0,0.1), scales=list(x=list(rot=90, 
cex=0.3),y=list(cex=0.5)) )

which.max(ff)
levelplot(ff, xlab="Position", las=2, ylab="Nuclei", main="BiNom Reads SA212 (Run ID A499H)", col.regions=rgb.palette(15), cuts=25, at=seq(20,1500,100), scales=list(x=list(rot=90, 
cex=0.3),y=list(cex=0.5)) )


dev.off()

#####################################################################################################


