# File to print the CG Panel out in sequence to the samples
# To visualize the output in graphical format

# based on the script my Steve McKinney (Nov 2015)
# modified not to check primer matching
# by Dr Damian Yap Jan 2016

setwd("/home/dyap/Projects/Takeda_T3/CG")

# Typical qPCR file (ignoring first 8 lines as preamble of qPCR output)

# Preprocessed with shell script
inputfile="160114_Takeda_T3_CGPanel_qPCR-processed.csv"
QC<-read.table(inputfile, header=TRUE,blank.lines.skip=TRUE, sep="\t", na.strings="NA", skip=7)

ffile=strsplit(inputfile,split="[.]")[[1]][1]

#Instead get names from file itself (Col 10) from SDS2.4 file
nnames<-nrow(table(QC$Sample))
names(table(QC$Sample))
samples<-names(table(QC$Sample))[1:nnames]
leg<-paste0(samples, collapse=", ")

require("nlme")
QC$sidf <- factor(QC$Sample)
QC$sidn <- as.numeric(QC$sidf)
QC$Primerssp <- gsub(":", "\n", QC$Detector)

trgd <- groupedData(QC$RQ ~ sidn | Primerssp, data = QC, order.groups = FALSE)

# legend checking
len<-length(samples)
for (i in 1:len ) {
	print(i)
	j<-as.character(unique(QC$sidf[grep(pattern=i,QC$sidn)]))
	print (j)
	}


pdffile=paste(ffile,"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)
plot(trgd, main=ffile, sub=leg, aspect = "fill", 
	par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), 
	ylim=c(-100,2000), as.table = TRUE)
dev.off()

pdffile=paste(paste(ffile,"low",sep="_"),"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)
plot(trgd, main=ffile, sub=leg, aspect = "fill", 
	par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), 
	ylim=c(-100,500), as.table = TRUE)
dev.off()

pdffile=paste(paste(ffile,"ultralow",sep="_"),"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)
plot(trgd, main=ffile, sub=leg, aspect = "fill", 
	par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), 
	ylim=c(-10,50), as.table = TRUE)
dev.off()

pdffile=paste(paste(ffile,"log",sep="_"),"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)
logRQ<-log(QC$RQ)
trgd <- groupedData(logRQ ~ sidn | Primerssp, data = QC, order.groups = 
	FALSE)
plot(trgd, main=ffile, sub=leg, aspect = "fill", 
	par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table 
	= TRUE)
dev.off()


# end
