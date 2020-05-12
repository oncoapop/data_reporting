# File to print the CG Panel out in sequence to the samples
# To visualize the output in graphical format

# based on the script my Steve McKinney (Nov 2015)
# modified not to check primer matching
# by Dr Damian Yap Jan 2016

library(RColorBrewer)

setwd("/home/dyap/Projects/Takeda_T3/siRNA")

# Typical qPCR file (ignoring first 8 lines as preamble of qPCR output)

# Preprocessed with shell script
#inputfile="siSRRM1_CG_panel.tsv"
inputfile="20160119_Takeda_SRRM1KD_CGPanel-processed.csv"
QC<-read.table(inputfile, header=TRUE,blank.lines.skip=TRUE, sep="\t", na.strings="NA", skip=0)

ffile=strsplit(inputfile,split="[.]")[[1]][1]

#Instead get names from file itself (Col 10) from SDS2.4 file
#nnames<-nrow(table(QC$Sample))
#names(table(QC$Sample))
#samples<-names(table(QC$Sample))[1:nnames]
#leg<-paste0(samples, collapse=", ")

require("nlme")
QC$Primerssp <- gsub(":", "\n", QC$Detector)

# Subsetting the data for analysis
sub="Calibrated NS-1"
NS1 <- subset(QC, .id=="Day2 Calibrated NS-1" | .id=="Day3A Calibrated NS-1" | .id=="Day3B Calibrated NS-1" | .id=="Day4 Calibrated NS-1", 
	select=c(Primerssp,RQ,Sample))
unsort <- subset(NS1, Sample=="Day2 NS-1 10nM" | Sample=="Day2 siSRRM1 10nM" | Sample=="Day3 NS-1 10nM A" | Sample=="Day3 siSRRM1 10nM A" |
			Sample=="Day2 NS-1 20nM" | Sample=="Day2 siSRRM1 20nM" | Sample=="Day3 NS-1 20nM A" | Sample=="Day3 siSRRM1 20nM A")

#sub="Calibrated NS-2"
#NS2 <- subset(QC, .id=="Day2 Calibrated NS-2" | .id=="Day3A Calibrated NS-2" | .id=="Day3B Calibrated NS-2" | .id=="Day4 Calibrated NS-2", 
#	select=c(Primerssp,RQ,Sample))
#unsort <- subset(NS2, Sample=="Day2 NS-2 10nM" | Sample=="Day2 siSRRM1 10nM" | Sample=="Day3 NS-2 10nM A" | Sample=="Day3 siSRRM1 10nM A" |
#			Sample=="Day2 NS-2 20nM" | Sample=="Day2 siSRRM1 20nM" | Sample=="Day3 NS-2 20nM A" | Sample=="Day3 siSRRM1 20nM A")

#sub="Repeat Calibrated NS-1"
#NS1 <- subset(QC, .id=="Day3 Calibrated NS-1" | .id=="Day3A Calibrated NS-1" | .id=="Day3B Calibrated NS-1" | .id=="Day4 Calibrated NS-1", 
#	select=c(Primerssp,RQ,Sample))
#unsort <- subset(NS1, Sample=="Day4 NS-1 10nM" | Sample=="Day4 siSRRM1 10nM" | Sample=="Day3 NS-1 10nM B" | Sample=="Day3 siSRRM1 10nM B" |
#			Sample=="Day4 NS-1 20nM" | Sample=="Day4 siSRRM1 20nM" | Sample=="Day3 NS-1 20nM B" | Sample=="Day3 siSRRM1 20nM B")

#sub="Repeat Calibrated NS-2"
#NS2 <- subset(QC, .id=="Day2 Calibrated NS-2" | .id=="Day3A Calibrated NS-2" | .id=="Day3B Calibrated NS-2" | .id=="Day4 Calibrated NS-2", 
#	select=c(Primerssp,RQ,Sample))
#unsort <- subset(NS2, Sample=="Day4 NS-2 10nM" | Sample=="Day4 siSRRM1 10nM" | Sample=="Day3 NS-2 10nM B" | Sample=="Day3 siSRRM1 10nM B" |
#			Sample=="Day4 NS-2 20nM" | Sample=="Day4 siSRRM1 20nM" | Sample=="Day3 NS-2 20nM B" | Sample=="Day3 siSRRM1 20nM B")

unsort$sidf <- factor(unsort$Sample)
unsort$sidn <- as.numeric(unsort$sidf)

SRRM1_1 <- unsort[with(unsort, order(sidn)), ]
ncol<-unique(SRRM1_1$sidn)

leginfo <- data.frame(  sidn = rep("", length(ncol)),
			Lab = rep("", length(ncol)),
                        Col = rep("", length(ncol)),
                        stringsAsFactors = FALSE)

# display.brewer.all()
col=brewer.pal(length(ncol), "Paired")

# legend checking
for (i in (unique(SRRM1_1$sidn)) ) {
	j<-as.character(unique(SRRM1_1$Sample[grep(pattern=i,SRRM1_1$sidn)]))
	leginfo$sidn[i]=i
	leginfo$Lab[i]=j
	leginfo$Col[i]=col[i]
	}

for (k in seq(ncol))
	{
	print(k)
	j<-as.character(unique(SRRM1_1$Sample[grep(pattern=k,SRRM1_1$sidn)]))
	print(j)
	l<-as.character(leginfo$Lab[grep(pattern=k,leginfo$sidn)])
	print(l)
	}
leginfo
legend=leginfo$Lab

# Write to PDF file (with cover page)

pdffile=paste(paste(paste(ffile,"log",sep="_"),sub,sep="_"),"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)

logRQ<-log(SRRM1_1$RQ)
trgd <- groupedData(logRQ ~ sidn | Primerssp, data = SRRM1_1, order.groups = FALSE)

#plot.new()
# This is the title page
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")


text(5, 8, ffile)
text(5, 7, sub)

text(5, 6, "QC of library by qPCR")

# This is the plot
plot(trgd, main=ffile, sub=sub, aspect = "fill", col = col,
        par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table
        = TRUE, pch=19)

legend("bottomright",legend=legend, fill=col, border=TRUE, bty="o", y.intersp = 1.0, cex=0.7)

dev.off()




#############

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
