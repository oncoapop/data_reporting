# Script to generate Venn diagram
# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# install.packages("gplots")
library(gplots)

# For EIF4A3 project directory
wd="/share/lustre/backup/dyap/Projects/eIF4A3_NMD/comparisons"

# Gene expression Up
Upfile="table_up.csv"
up=paste(wd,Upfile,sep="/")

# Read in csv
OurUp<-read.csv(file=up, stringsAsFactors=FALSE)
# any 2 out of 4 conditions (ie Hela/HCT-116 T-202/T-595)
OurUp$sum=rowSums(OurUp[3:6])
OurUpSet<-OurUp[which (OurUp$sum >= 2),]
OurUpGenes<-unique(OurUpSet$Name)


# Mao (mouse) # converted by matching names using shell script human2mouse.sh
mfile="Mao_2016_Eif4a3_human.csv"
mu=paste(wd,mfile,sep="/")

# Read in csv
mut<-read.csv(file=mu, stringsAsFactors=FALSE)
colnames(mut)[3]<-"logFC_Eif4a3_Emx1Cre"

# subsetting the data greater than 2 fold and significant
DeltaMu<-mut[ which( mut$FDR < 0.05), ]

# there are 148 genes significantly altered in E10.5 embryos when Eif4a3 is KD
nrow(DeltaMu)

# of which there are 59 genes with no human counterpart
noHuman<-DeltaMu[DeltaMu$HumanGene %in% "N/A",]
nrow(noHuman)

# 89 relevant genes of which 
MaoUp<-DeltaMu[DeltaMu$logFC_Eif4a3_Emx1Cre > 1,]
nrow(MaoUp) # 63
MaoDown<-DeltaMu[DeltaMu$logFC_Eif4a3_Emx1Cre < 1,]
nrow(MaoDown) # 85


# NO overlap of genes
overlap1<-OurUpSet[OurUpSet$Name %in% MaoUp$GeneName,]
overlap1$GeneName

##################################
# Use KEGG pathway from Mao at al
# and our set of any two
# MISO high dose only
pathMonoUpfile="gene_pathways_monotonically_increasing_info.csv"
pathMonoUp<-read.csv(file=paste(wd,pathMonoUpfile,sep="/"), stringsAsFactors=FALSE)
SigPathUp<-pathMonoUp[pathMonoUp$gene_name %in% OurUpSet$Name,]
SigPathUp$Name<-gsub("GO_","",SigPathUp$GO_process)
pathUp<-unique(gsub("GO_","",SigPathUp$GO_process))


pathMaoMufile="Mao_mouse_GOlt005.txt"
pathMaoMu<-as.list(read.csv(file=paste(wd,pathMaoMufile,sep="/"), stringsAsFactors=FALSE, header=FALSE))
pathMao<-as.list(unique(pathMaoMu))

input<-list(pathUp,pathMao[[1]])
Reduce(intersect, list(pathUp,pathMao))

par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
#box("figure",lty="solid", col="green")

plot.new()
#tmp <- venn(input, showSetLogicLabel=TRUE)
tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
all4<-attr(tmp, "intersections")$'1111'
# The interestion of all conditions
overlap<-attr(tmp, "intersections")


################
msdir="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
ms=paste(msdir,"CLK2_allIP.txt",sep="/")

# For CG factors
cgdir="/home/dyap/Projects/Takeda_T3/CG_factors"
cg=paste(cgdir,"CG-factors",sep="/")

# For TiO2 phosphopeptide expt
tidir="/home/dyap/Projects/Takeda_T3/TiO2"
phos=paste(tidir,"Filtered_top42_splice.txt",sep="/")

msf<-read.table(file=ms)
cgf<-read.table(file=cg)
phosf<-read.table(file=phos)

input<-list(msf,cgf,phosf)

out="/home/dyap/Projects/Takeda_T3"
pdf(file=paste(out,"Venn_Intersect.pdf",sep="/"))
par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
#box("figure",lty="solid", col="green")

#plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
# all4<-attr(tmp, "intersections")$'1111'
# The interestion of all conditions
overlap<-attr(tmp, "intersections")

# This section prints out the genes in the intersection
        createCounter <- function(value) { function(j) { value <<- value+j} }
        createRowCounter <- function(value) { function(j) { value <<- value-r} }
        counter <- createCounter(1)
        rowc <- createRowCounter(30)
	name <- createCounter(1)

#        mtext("IP-MS ectopic expressed FLAG-tagged CLKs", side=3)
#        mtext("Venn diagram", side=3)
#	text(40,350, "CLK1")
	text(200,380, "Decreased phospho-peptide in TiO2")
     	text(55,80, "High Conf interactor of CLK2", srt=-45)
     	text(335,80, "Factors w/motif enrich in CG", srt=45)
#	text(365,350, "CLK3_cs")

dev.off()

outfile=paste(out,"Intersects.txt",sep="/")

# works but no header
#writeLines(unlist(lapply(overlap, paste, collapse=" ")))

sink(outfile) 
lapply(overlap, print) 
sink() 



########################## 
rows=ceiling((length(all4))/10)

	for (r in seq(rows)){
		print(r)
		y <- ((rows-r)*8)-20
              for (j in 1:10) {x <- counter(40)+r; text(x, y, all4[name(1)], srt=35, cex=0.6); print(x)}
	        counter <- createCounter(1)
	        name <- createCounter(r*10)
		print("hello world")
		print(y)
		}


