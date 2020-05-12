# Script to visualize data
library(foreign)
library(RColorBrewer)
library(lattice)
library(colorspace)
library(plyr)
library(ggrepel)
library(ggplot2)

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)


input <- read.table(file = "/home/dyap/Projects/Takeda_T3/siRNA/siRNA-test", stringsAsFactors = FALSE, sep="\t", header=TRUE)
gene <-subset(input, siRNA!="siNT"& siRNA!="HCT116")

# Select n colors
n <- length(table(gene$Plate))
# display.brewer.all()
rowcol<-brewer.pal(n,"Set1")

plateinfo <- data.frame(  Plate = rep("", nrow(gene)),
                        siRNA = rep("", nrow(gene)),
                        RQ = rep("", nrow(gene)),
                        geneID = rep("", nrow(gene)),
                        Col = rep("", nrow(gene)),
                        stringsAsFactors = FALSE)

for (a in seq(nrow(gene)))
        {

	rownames(plateinfo)[a] <- paste(gene$siRNA[a],gene$Detector[a],sep="_")
        plateinfo$Plate[a] <- gene$Plate[a]
        plateinfo$siRNA[a] <- gene$siRNA[a]
        plateinfo$RQ[a] <- gene$RQ[a]
        plateinfo$GeneID[a] <- strsplit(gene$siRNA,split="-")[[a]][1]

        for ( b in seq(n) )
                {
                if (gene$Plate[a] == names(table(gene$Plate)[b])) color=rowcol[b]
                }

        plateinfo$Col[a] <- color

        }


legend <-  data.frame(  Plate = rep("", n),
                        Col = rep("", n),
                        stringsAsFactors = FALSE)

for ( c in seq(n) )
	{
	legend$Plate[c] <- names(table(gene$Plate))[c]
	legend$Col[c] <- rowcol[c]
	}

gene<-mutate(gene, siRNAs=ifelse(gene$RQ < 0.4, "RQ < 0.4 (n=103)", "Not selected (n=132)"))

table(gene$siRNAs)

totalgenes<-length(unique(gene$Detector))
totlab<-paste("Total no. of unique genes assayed =",totalgenes, sep=" ")

selgene<-gene[gene$siRNAs == "RQ < 0.4 (n=103)",]
totselgenes<-length(unique(selgene$Detector))
sellab<-paste("No. of unique genes, RQ < 0.4 =",totselgenes, sep=" ")

p = ggplot(gene, aes(x=reorder(siRNA,RQ), y=RQ), pch=20) + 
    geom_point(aes(col=siRNAs)) +
    theme_light()+
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
    scale_color_manual(values=c("grey90", "grey50")) +
    scale_y_continuous(limits=c(0,5)) +
    theme(legend.position = c(0.20, 0.85)) +
    ggtitle("Effect of siRNA on Expression of Domains associated with CG-enriched motifs") +
    theme(plot.title=element_text(face="bold", size=10)) +
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
    labs(x = "siRNA", y="Normalized Expression (Fold change, RQ)") +
    theme(axis.text.x=element_text(size=4, angle = 90, hjust = 1),
        axis.title=element_text(size=12,face="bold")) +
    annotate("text", x = 150, y = 4.5, label = totlab) +
    annotate("text", x = 150, y = 0.2, label = sellab, fontface = "bold") +
    geom_hline(yintercept=0.4, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")

# Plot with entire dataset minus 2 outliers
pdf("siRNA_Summary.pdf", width=15, height=8)
p
dev.off()


# Plot with entire dataset minus 1 outlier
pdf("RQbyplate.pdf", width=15, height=8)

plot(plateinfo$RQ, main="Plate information on RQ", xaxt='n', xlab="siRNA",
	ylab="Fold expression change (normalized)", pch=18, col=plateinfo$Col, cex=1.2, ylim=c(0,5))

legend(150,3,c(legend$Plate), pch=18, col=legend$Col)

axis(1, at=seq(nrow(gene)), labels=gene$siRNA, las=2, cex.axis=0.5)

dev.off()


# only those <1

samples <-subset(plateinfo, RQ < 1)

pdf("RQbyPlate2.pdf", width=15, height=8)

plot(samples$RQ, main="Plate information on RQ", xaxt='n', xlab="siRNA",
        ylab="Fold expression change (normalized)", pch=18, col=samples$Col, cex=1.2, ylim=c(0,1))

legend(50,1,c(legend$Plate), pch=18, col=legend$Col)

axis(1, at=seq(nrow(samples)), labels=samples$siRNA, las=2, cex.axis=0.5)

dev.off()

samples <-subset(plateinfo, RQ < 0.4)

pdf("RQbyPlate3.pdf", width=15, height=8)

plot(samples$RQ, main="Plate information on RQ", xaxt='n', xlab="siRNA",
        ylab="Fold expression change (normalized)", pch=18, col=samples$Col, cex=1.2, ylim=c(0,0.5))

legend(70,0.5,c(legend$Plate), pch=18, col=legend$Col)

axis(1, at=seq(nrow(samples)), labels=samples$siRNA, las=2, cex.axis=0.5)

dev.off()

# get no of genes

geneinfo <- data.frame( siRNA = rep("", nrow(samples)),
                        gene = rep("", nrow(samples)),
                        RQ = rep("", nrow(samples)),
                        Plate = rep("", nrow(samples)),
                        Amp = rep("", nrow(samples)),
                        Col = rep("", nrow(samples)),
                        stringsAsFactors = FALSE)

for ( d in seq(nrow(samples)) )
        {
        geneinfo$gene[d] <- strsplit(samples$siRNA,split="-")[[d]][1]
        geneinfo$siRNA[d] <- samples$siRNA[d]
        geneinfo$Plate[d] <- samples$Plate[d]
        geneinfo$RQ[d] <- samples$RQ[d]
        }

KD<-sub("^si", "", names(table(geneinfo$gene)))
CGenrich<-sub("^si", "", names(table(plateinfo$GeneID)))

all_MS <- read.csv(file = "/home/dyap/Projects/Takeda_T3/siRNA/MS-data_R.csv", stringsAsFactors = FALSE, header=TRUE)

cut=1
p=0.05
pep=2
enrich_MS <-subset(all_MS, (rN1 > cut || rN2 > cut ) & (rC1 > cut || rC2 > cut ) & rU > cut & All_p.adj < p & pepNum > pep)
MS<-enrich_MS$Gene

########################################################################################################
# This marks the genes of interest which are at the interface of the CG factor which have good KD
# as well as the enriched by CLK2 in MS
si<-intersect(CGenrich,KD)
bin<-intersect(CGenrich,MS)

mark<-intersect(si,bin)

# display.brewer.all()
rowcol<-brewer.pal(11,"RdYlGn")


geneinfo <- data.frame( siRNA = rep("", nrow(samples)),
                        gene = rep("", nrow(samples)),
                        RQ = rep("", nrow(samples)),
                        Plate = rep("", nrow(samples)),
                        Col = rep("", nrow(samples)),
                        stringsAsFactors = FALSE)

for (a in seq(nrow(samples)))
        {

        geneinfo$gene[a] <- sub("^si", "", strsplit(samples$siRNA,split="-")[[a]][1])
        geneinfo$siRNA[a] <- samples$siRNA[a]
        geneinfo$Plate[a] <- samples$Plate[a]
        geneinfo$RQ[a] <- samples$RQ[a]

	color<-rowcol[6]

        for ( b in seq(length(mark)) )
                {
                if (geneinfo$gene[a] == names(table(mark)[b])) color=rowcol[1]
                }

        geneinfo$Col[a] <- color


        }

pdf("RQbyMS.pdf", width=15, height=8)

plot(geneinfo$RQ, main="Location of MS enriched genes", xaxt='n', xlab="Gene targeted by siRNA",
        ylab="Fold expression change (normalized)", pch=18, col=geneinfo$Col, cex=1.2, ylim=c(0,0.5))

legend(70,0.5,c("Enriched in CLK2 MS","Not enriched"), pch=18, col=c(rowcol[1],rowcol[6]))

axis(1, at=seq(nrow(geneinfo)), labels=geneinfo$gene, las=2, cex.axis=0.5)

dev.off()

sl<-subset(geneinfo, geneinfo$Col==rowcol[1])

write.table(sl, file="shot-list.csv",sep=",")

########################################################################################################
# This marks the genes of interest which are at the interface of the CG factor which have good KD
si<-intersect(CGenrich,KD)

n<-length(si)
# display.brewer.all()
rowcol<-colorRampPalette(brewer.pal(11,"Spectral"))(n)


geneinfo <- data.frame( siRNA = rep("", nrow(samples)),
                        gene = rep("", nrow(samples)),
                        RQ = rep("", nrow(samples)),
                        Plate = rep("", nrow(samples)),
                        Col = rep("", nrow(samples)),
                        stringsAsFactors = FALSE)

for (a in seq(nrow(samples)))
        {

        geneinfo$gene[a] <- sub("^si", "", strsplit(samples$siRNA,split="-")[[a]][1])
        geneinfo$siRNA[a] <- samples$siRNA[a]
        geneinfo$Plate[a] <- samples$Plate[a]
        geneinfo$RQ[a] <- samples$RQ[a]

        for ( b in seq(length(si)) )
                {
                if (geneinfo$gene[a] == names(table(si)[b])) color=rowcol[b]
                }

        geneinfo$Col[a] <- color


        }


pdf("RQbyGene.pdf", width=15, height=8)

plot(geneinfo$RQ, main="Coloured by genes with effective siRNA", xaxt='n', xlab="Gene targeted by siRNA",
        ylab="Fold expression change (normalized)", pch=18, col=geneinfo$Col, cex=1.2, ylim=c(0,0.5))

legend(5,0.5,names(table(si)), pch=18, col=rowcol, ncol=3, cex=0.8)

axis(1, at=seq(nrow(geneinfo)), labels=geneinfo$gene, las=2, cex.axis=0.5)

dev.off()

fsl<-geneinfo[order(geneinfo$gene),]

write.table(fsl, file="factorKD_shot-list.csv",sep=",")

########################################################################################################
# on command line R CMD javareconf
#install.packages("rJava", dependencies=TRUE)
#install.packages("XLConnect")

# These modules is done on MOMAC14 as it has the java installed
library(XLConnect)

setwd("/Users/dyap/Documents/Collboration\ Projects/Takeda-Splicing/Validation/CG\ factor\ siRNA\ KD/")

# Excel file and sheet names (>1 sheet)
fname="/Users/dyap/Documents/Collboration\ Projects/Takeda-Splicing/Validation/CG\ factor\ siRNA\ KD/20151210/RBP_CG_Assay\ G-CUSTOM-176687.xls"
sname="Sequences_by_Ord_wShippingRack_"
srow=3
scol=1
erow=263
ecol=9

excel <- readWorksheetFromFile(file=fname, sheet = sname, header = TRUE, startRow = srow, startCol = scol, endCol = ecol)
indf<-excel[complete.cases(excel),]

#indfnoc <- excel[!grepl("ON-TARGET", excel$Gene.Symbol), ]
indfnoc<-excel[complete.cases(excel),]

indfnoc$Gene.Symbol
gsrle <- rle(indfnoc$Gene.Symbol)

table(gsrle$values, useNA = "always")
table(table(gsrle$values))
seqnos <- unlist(lapply(gsrle$lengths, function(x) seq(x)))
gnms <- unlist(lapply(seq(along = gsrle$values), function(x) {rep(gsrle$values[x], gsrle$lengths[x]) } ) )
gnmsnos <- paste(gnms, seqnos, sep = "-")

indf$Gene.Symbols.N <- rep(NA_character_, nrow(indf))
indf$Gene.Symbols.N <- gnmsnos

###############################

sl<-read.table(file="factorKD_short-list.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

sl$siRNA.N <-sub("^si", "", sl$siRNA)

indf$Pos<-paste(indf$Plate,indf$Well,sep='-')

sum<-merge(sl,indf, by.x="siRNA.N", by.y="Gene.Symbols.N") 

# Check that this should be NULL
setdiff(sl$siRNA.N,sum$siRNA.N)

# Select the relevant columns
dat<-sum[,c(1,5,12,16,17)]

# make sure the duplicate wells are removed
dd<-dat[!duplicated(dat[, c("Pos")]), ]

# Order according to plate and well position
outdf<-dd[with(dd, order(Pos)), ]


write.csv(outdf,file="CGsiRNA.csv")

