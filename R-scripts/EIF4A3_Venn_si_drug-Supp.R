# Script to generate Venn diagram
# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# install.packages("gplots")
# install.packages("gtools", type="source")
# install.packages('venneuler',,'http://www.rforge.net/')
# If required, update rJava
# R CMD javareconf -e
# install.packages("rJava")
# install.packages('rJava', dependencies=TRUE, repos='http://cran.r-project.org/')

library(gplots)
library(RColorBrewer)
library(colorspace)
library(venneuler)

# For EIF4A3 project directory
# script to compare siRNA genes 
# this the high amd low dose of drugs by MISO only (main Fig 2)

wd="/home/dyap/Projects/EIF4A3_paper/comparisons/compare_genes"
wd2="/home/dyap/Projects/EIF4A3_paper/comparisons"
# Please see README in this directory (wd) for the source of the files 

misodir="/home/dyap/Projects/EIF4A3_paper/comparisons/MISO"
# Ran script to annotate events to gene names and make them unique, so less than events
misosi48file=paste(misodir,"miso_si48_genes",sep="/")
miso48_siRNA<-read.table(file=misosi48file, stringsAsFactors=FALSE, header=FALSE)
# This is the criteria that AliB told me to filter by


#misohighfile=paste(misodir,"miso_Hela_202.highdose",sep="/")
# Ran script to annotate events to gene names and make them unique, so less than events
#misohighfile=paste(misodir,"miso_Hela_202high_genes",sep="/")
#miso_HeLa_T202high<-read.table(file=misohighfile, stringsAsFactors=FALSE, header=FALSE)
# 321 gene unique names

# See README in $wd for source of files
misohighfile=paste(wd,"Hela_202.highdose",sep="/")
# Use old file which I assume to be correct given the similarities with Fig 2 (which uses events rather than genes)
miso_HeLa_T202high<-read.table(file=misohighfile, stringsAsFactors=FALSE, header=FALSE)
# 1100 gene unique names

# See README in	$wd for	source of files
misolowfile=paste(wd,"Hela_202.lowdose",sep="/")
# Use old file which I assume to be correct given the similarities with	Fig 2 (which uses events rather	than genes)
miso_HeLa_T202low<-read.table(file=misolowfile, stringsAsFactors=FALSE, header=FALSE)
# 427 gene unique names

# See README in $wd for source of files
SGfile=paste(wd2,"SG-genes.txt",sep="/")
SG_genes<-read.table(file=SGfile, stringsAsFactors=FALSE, header=FALSE)
### EVENTS NOT GENES - change to GENE (using shell prelim script to annotate genes)

# Fig PRELIM
input<-list(miso48_siRNA,miso_HeLa_T202high,miso_HeLa_T202low)
Reduce(intersect, list(miso48_siRNA,miso_HeLa_T202high,miso_HeLa_T202low))
names(input) <- c("Hela \n(siRNA, \n48 hr)","HeLa \n(T-202, \nhigh dose)","HeLa (T-202, low dose)")

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()

pdffile=paste(wd,"Venn_Fig2c",sep="/")
pdf(file=pdffile)
tmp <- venn(input, show.plot=TRUE)
dev.off()

# The Redraw for better labelling for Gene List for Supplementary Data
names(input) <- c("siRNA","T-202_High","T-202_Low")
tmp <- venn(input, show.plot=FALSE)
overlap<-attr(tmp, "intersections")

outfile=paste(wd,"GeneTable_si_vs_HighLow_Dose.txt",sep="/")
capture.output(overlap, file = outfile)


####################################################
####################################################

# Testing overlap with 63 SG genes only 1 overlap (not interesting);
# use existing table 

#input<-list(overlap$`siRNA:T-202_High:T-202_Low`,SG_genes)
#Reduce(intersect, list(overlap$`siRNA:T-202_High:T-202_Low`,SG_genes))
input<-list(miso48_siRNA,miso_HeLa_T202high,miso_HeLa_T202low,SG_genes)
Reduce(intersect, list(miso48_siRNA,miso_HeLa_T202high,miso_HeLa_T202low,SG_genes))
names(input) <- c("Hela \n(siRNA, \n48 hr)","HeLa \n(T-202, \nhigh dose)","HeLa (T-202, \nlow dose)","SG Genes")

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
#pdffile=paste(wd,"Venn_Supp_new6",sep="/")
pdf(file=pdffile)
#tmp <- venn(input, showSetLogicLabel=TRUE)
tmp <- venn(input, show.plot=TRUE)

overlap_2<-attr(tmp, "intersections")

plot.new()
#tmp <- venn(input, showSetLogicLabel=TRUE)
tmp <- venn(input, show.plot=TRUE, intersections=TRUE)
dev.off()
tmp2 <- venn(input, show.plot=FALSE, intersections=FALSE)
# tmp2 is important for the generation of the colour venn plot

#        mtext("Overlap of siRNA studies", side=3)
#	text(50,350, "Wang et al (MISO)")
#        text(150,380, "Wang et al (DIFF)")
#        text(280,380, "siRNA 48hr (MISO)")
#	text(360,350, "siRNA 48hr (VAST)")


# The interestion of all conditions
overlap_allsi<-attr(tmp, "intersections")

# Overlap of genes in at least 3/4 cases
overlap_allsi_3<-unique(do.call(c,list(overlap_allsi$'A:B:C',overlap_allsi$'B:C:D',overlap_allsi$'A:C:D',overlap_allsi$'A:B:D',overlap_allsi$'A:B:C:D')))
# 64 genes

capture.output(summary(overlap_allsi), file = "Overlaps_ext_siRNA_summary.txt")
capture.output(overlap_allsi, file = "Overlaps_ext_siRNA.txt")

############%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


