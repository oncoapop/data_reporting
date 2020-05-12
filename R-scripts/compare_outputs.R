
# source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
library("VariantAnnotation")
library("IRanges")
library("GenomicRanges")
library(foreign)
library(lattice)
require(Heatplus)
library(limma)
library(gplots)
library(RColorBrewer)

library(XLConnect)

#####

# To run this script change the setwd()
inputpath="/home/dyap/Projects/Tumour_Evol/Xenodrug/SA532-A6VL1"

file="xeno_SA532_A6VL1_target_list_stats.xls"

inputfile=paste(inputpath,file,sep="/")

sname="xeno_SA532_A6VL1_target_VarAlle"
srow=1
scol=1
erow=192
ecol=9

excel <- readWorksheetFromFile(file=inputfile, sheet = sname, header = TRUE, startRow = srow, startCol = scol, endCol = ecol)
indf<-excel[complete.cases(excel),]


#####################################################################################################
# To run this script change the setwd()
setwd("/home/dyap/Projects/Tumour_Evol/Xenodrug")
sample="S"
ID="SA532"
#run="AGALJ"
#run="AGAJY"
run="AUBU6"
#check="depth"
check = "freq"
#check = "calls"

# Outputs
dir=paste(paste("/home/dyap/Projects/Tumour_Evol/Xenodrug",ID,sep="/"),run, sep="/")

inname=paste(ID,run,sep="-")
runname=paste(run,"run",sep="-")
exptname=paste(ID,runname,sep="-")

# For output
if (check == "depth") filename=paste(exptname,"reads",sep="_")
if (check == "freq") filename=paste(exptname,"freq",sep="_")
if (check == "calls") filename=paste(exptname,"calls",sep="_")

# For input file name (from QC_2)
if (check == "depth") infilename=paste(inname,"reads",sep="_")
if (check == "freq") infilename=paste(inname,"freq",sep="_")
if (check == "calls") infilename=paste(inname,"calls",sep="_")

# For input file name (from QC_2)
if (check == "depth") sfilename=paste(ID,"reads",sep="-")
if (check == "freq") sfilename=paste(ID,"freq",sep="-")
if (check == "calls") sfilename=paste(ID,"calls",sep="-")


sfile=paste(paste(dir,sfilename,sep="/"),"csv",sep=".")
infile=paste(paste(dir,infilename,sep="/"),"csv",sep=".")
output=paste(dir,filename,sep="/")

csvfile=paste(output, "csv", sep=".")
pdffile=paste(output, "pdf", sep=".")
pdffile2=paste(output, "cluster.pdf", sep="")
pdffile3=paste(output, "cluster2.pdf", sep="")
vennfile=paste(output, "Venn.pdf", sep="-")

title=filename
xaxislab2=paste(ID, "Sample", sep=" ")

prop<-read.table(file=infile ,sep=",", header=TRUE)


####################################################################################################

sum2 <- read.csv(file=sfile)

ID<-strsplit(as.character(sum2$ID),split="_")
int<-do.call(rbind, ID)

cp<-int[,1]
ID2<-strsplit(as.character(cp),split=":")

int<-do.call(rbind, ID2)

sum2$target<-paste(gsub("chr","",int[,1]),int[,2],sep="_")

comb <- 
merge(indf,sum2,by="target")
#merge(sum1,sum2,by="ID",all.y = "TRUE")

# remove col ID
comb$ID<-NULL

#####################################
# Preparing the data for clustering (removing NAs)

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(comb[2:ncol(comb)])

# col headers - sample names
names(comb)

# Filters out all the positions that failed in all samples
indiv=length(comb)-1
filt <-rownames(ef[rowSums(is.na(ef))==indiv,])
filt

# Label rownames with ID
rownames(ef) <- comb$ID

colnames(ef)
ff<-ef[,order(as.numeric(colnames(ef)))]

# If filt=NULL then all primers work so NA = zero
ff[is.na(ff)] <- 0

title="Comparison of SA532 run AUBU6 and A6VL1"
hmcols <- colorRampPalette(brewer.pal(11,"Spectral"))(100)

pdffile="SA532-AUBU6vsA6VL1.pdf"

setwd("/home/dyap/Projects/Tumour_Evol/Xenodrug")
pdf(pdffile, width=7, height=8)

heatmap.2(ff, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.8, cexRow=0.6, col = hmcols, trace="none", dendrogram="col")

dev.off()

