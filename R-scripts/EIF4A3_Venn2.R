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
# Convert from python script to R script to plot the update versions

wd="/home/dyap/Projects/EIF4A3_paper/comparisons"
wd3="/home/dyap/Projects/EIF4A3_paper"

fileSG=paste(wd,"SG-genes.txt",sep="/")
SG_genes<-read.table(file=fileSG, stringsAsFactors=FALSE, header=FALSE)

# Own data (MISO)
misodir="/home/dyap/Projects/EIF4A3_paper/comparisons/MISO"
# misosi48file=paste(misodir,"miso_48h_0.1_10",sep="/")
# Ran script to annotate events to gene names and make them unique, so less than events
misosi48file=paste(misodir,"miso_si48_genes",sep="/")
miso48_siRNA<-read.table(file=misosi48file, stringsAsFactors=FALSE, header=FALSE)

#misohighfile=paste(misodir,"miso_Hela_202.highdose",sep="/")
# Ran script to annotate events to gene names and make them unique, so less than events
misohighfile=paste(misodir,"miso_Hela_202high_genes",sep="/")
miso_HeLa_T202high<-read.table(file=misohighfile, stringsAsFactors=FALSE, header=FALSE)

# Own data (VAST)
vastdir="/home/dyap/Projects/EIF4A3_paper/comparisons/VAST"
# Ran script to annotate events to gene names and make them unique, so less than events
#vastsi48file=paste(vastdir,"vast_48h_0.1_10",sep="/")
vastsi48file=paste(vastdir,"vast_si48_genes",sep="/")
vast48_siRNA<-read.table(file=vastsi48file, stringsAsFactors=FALSE, header=FALSE)

# Ran script to annotate events to gene names and make them unique, so less than events
#vasthighfile=paste(vastdir,"vast_Hela_202.highdose",sep="/")
vasthighfile=paste(vastdir,"vast_Hela_202high_genes",sep="/")
vast_HeLa_T202high<-read.table(file=vasthighfile, stringsAsFactors=FALSE, header=FALSE)

### EVENTS NOT GENES - change to GENE (using shell prelim script to annotate genes)
# Or grab from Alborz file (see READ me in directory for source)

# Fig PRELIM
input<-list(miso48_siRNA,miso_HeLa_T202high,vast48_siRNA,vast_HeLa_T202high)
Reduce(intersect, list(miso48_siRNA,miso_HeLa_T202high,vast48_siRNA,vast_HeLa_T202high))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
#tmp <- venn(input, showSetLogicLabel=TRUE)
tmp <- venn(input, show.plot=TRUE, showSetLogicLabel=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

capture.output(summary(overlap), file = "Overlaps.txt")
####################################################
####################################################

# Testing overlap with 63 SG genes 
# use existing table 

input<-list(allcondoverlap,SG_genes)
Reduce(intersect, list(allcondoverlap,SG_genes))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)
overlap_2<-attr(tmp, "intersections")

# Get overlap of siRNA in both 

####################################################
####################################################


# Supp Fig 6f
# comparing with siRNAs:
#./plot_venn.py --inputFiles Miso_EIF4A3.genes DiffSplice_EIF4A3.genes ../Miso_siRNA/48h_genenames --inputNames Miso_EIF4A3_others Diffsplice_EIF4A3_others $wd"/siRNA_new" --inputCols 1 
#--outputpref MisoEIF4A3Others_DiffSpliceOthers_siRNAOur48h

input<-list(WangMiso$V1,WangDiffSpl$V1,miso48_siRNA,vast48_siRNA)
Reduce(intersect, list(WangMiso$V1,WangDiffSpl$V1,miso48_siRNA,vast48_siRNA))
#names(input) <- c("Wang et al (MISO)","Wang et al (DIFF)","siRNA 48hr (MISO)","siRNA 48hr (VAST)")
names(input) <- c("Wang-M","Wang-D","siRNA-M","siRNA-V")

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
#tmp <- venn(input, showSetLogicLabel=TRUE)
pdffile=paste(wd,"Venn_Fig6f",sep="/")
pdf(file=pdffile)
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

# Drawing with colour

# define labels
labels <- c("Wang-M","Wang-D","siRNA-M","siRNA-V")
# define colors
colors <- rainbow_hcl(15)
# name colors with labels
names(colors) <- labels

# Data generated from the tmp2 above in a 4 x 4 matrix
mat <- matrix(c(tmp2[9],tmp2[13],tmp2[11],tmp2[10],tmp2[13],tmp2[5],tmp2[7],tmp2[6],tmp2[11],tmp2[7],tmp2[3],tmp2[4],tmp2[10],tmp2[6],tmp2[4],tmp2[2]), nrow=4,ncol=4)
colnames(mat) <- labels
rownames(mat) <- labels

make_venn <- function(dat, title="",fill_cols=NULL
                      , alpha=.7, title_font=1){
  vv <- venneuler::venneuler(dat)
  if (!is.null(fill_cols)){
    plot(vv, main=title, font.main = title_font
    , col=fill_cols[match(vv$labels, names(colors))], alpha=alpha )
  } else {
    plot(vv, main=title, font.main = title_font, alpha=alpha)
  }
}


make_venn(make_dat(mat),title="Happy Birthday Venn!"
            ,fill_cols=colors, title_font=2)

####################################################
####################################################

# Fig 16 a
input<-list(WangMiso$V1,HeLa_monotonic,HCT116_monotonic)
Reduce(intersect, list(WangMiso$V1,HeLa_monotonic,HCT116_monotonic))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
all3<-attr(tmp, "intersections")$'111'

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

# Fig 16 a PLUS SG
input<-list(WangMiso$V1,HeLa_monotonic,HCT116_monotonic,SG_genes)
Reduce(intersect, list(WangMiso$V1,HeLa_monotonic,HCT116_monotonic,SG_genes))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
#all3<-attr(tmp, "intersections")$'111'
attr(tmp, "intersections")$'0011'
attr(tmp, "intersections")$'0111'
attr(tmp, "intersections")$'0101'
attr(tmp, "intersections")$'1001'


# The interestion of all conditions
overlap<-attr(tmp, "intersections")

####################################
# Fig 16 b
input<-list(WangDiffSpl$V1,HeLa_detected,HCT116_detected)
Reduce(intersect, list(WangDiffSpl$V1,HeLa_detected,HCT116_detected))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

# Fig 16 b PLUS SG
input<-list(WangDiffSpl$V1,HeLa_detected,HCT116_detected,SG_genes)
Reduce(intersect, list(WangDiffSpl$V1,HeLa_detected,HCT116_detected,SG_genes))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

##########################################


####################################
# Fig 16 c
input<-list(WangDiffSpl$V1,HeLa_monotonic,HCT116_monotonic)
Reduce(intersect, list(WangDiffSpl$V1,HeLa_monotonic,HCT116_monotonic))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

# Fig 16 b PLUS	SG
input<-list(WangDiffSpl$V1,HeLa_monotonic,HCT116_monotonic,SG_genes)
Reduce(intersect, list(WangDiffSpl$V1,HeLa_monotonic,HCT116_monotonic,SG_genes))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")

##########################################

# f 
input<-list(WangMiso$V1,WangDiffSpl$V1,siRNA48_new$Gene)
Reduce(intersect, list(WangMiso$V1,WangDiffSpl$V1,siRNA48_new$Gene))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
box("figure",lty="solid", col="green")

plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all conditions
overlap<-attr(tmp, "intersections")


