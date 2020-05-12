
# Script to visualize data
library(foreign)
library(RColorBrewer)
library(lattice)
library(colorspace)

rgb.palette <- colorRampPalette(c("dark green", "red"), space = "rgb")

gene <- read.table(file = "/home/dyap/Projects/FFPE/gene.txt", stringsAsFactors = FALSE, sep="\t", header=TRUE)

plot(gene$SampleID, gene$GENE_NAME, xlab = "Sample_ID", ylab = "Mutation")


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

