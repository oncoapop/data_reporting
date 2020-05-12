# File to print the CG Panel out in sequence to the samples
# To visualize the output in graphical format

# based on the script my Steve McKinney (Nov 2015)
# modified not to check primer matching
# by Dr Damian Yap Jan 2016

library(RColorBrewer)

#setwd("/home/dyap/Projects/Takeda_T3/siRNA")
setwd("/share/lustre/backup/dyap/Projects/Takeda_T3/CG")

# Typical qPCR file (ignoring first 8 lines as preamble of qPCR output)

# Preprocessed with shell script
#inputfile="siSRRM1_CG_panel.tsv"
#inputfile="20160119_Takeda_SRRM1KD_CGPanel-processed.csv"
inputfile="siPanel_Plate1_CGAssay.sdm-Result\ Data.txt"

QC<-read.table(inputfile, header=TRUE, blank.lines.skip=TRUE, sep="\t", na.strings="NA", skip=0)

ffile=strsplit(inputfile,split="[.]")[[1]][1]

require("nlme")
QC$Primerssp <- gsub(":", "\n", QC$Detector)

# Subsetting the data for analysis
sub="siCG-factors_Plate1"

unsort<-QC

unsort$sidf <- factor(unsort$Sample)
unsort$sidn <- as.numeric(unsort$sidf)

siCG <- unsort[with(unsort, order(sidn)), ]
ncol<-unique(siCG$sidn)

leginfo <- data.frame(  sidn = rep("", length(ncol)),
			Lab = rep("", length(ncol)),
                        Col = rep("", length(ncol)),
                        stringsAsFactors = FALSE)

# display.brewer.all()
col=brewer.pal(length(ncol), "Paired")

# legend checking
for (i in (unique(siCG$sidn)) ) {
	j<-as.character(unique(siCG$Sample[grep(pattern=i,siCG$sidn)]))
	leginfo$sidn[i]=i
	leginfo$Lab[i]=j
	leginfo$Col[i]=col[i]
	}

for (k in seq(ncol))
	{
	print(k)
	j<-as.character(unique(siCG$Sample[grep(pattern=k,siCG$sidn)]))
	print(j)
	l<-as.character(leginfo$Lab[grep(pattern=k,leginfo$sidn)])
	print(l)
	}
leginfo
legend=leginfo$Lab

# Write to PDF file (with cover page)

# For ranging 
temp<-sub

# repeat plots start here
sub=paste(temp,"summary",sep="_")
#range=c(0,20)


pdffile=paste(paste(paste(ffile,"log",sep="_"),sub,sep="_"),"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)

# log data
logRQ<-log(siCG$RQ)
trgd <- groupedData(logRQ ~ sidn | Primerssp, data = siCG, order.groups = FALSE)

# linear data
#lin<-siCG$RQ
#trgd <- groupedData(lin ~ sidn | Primerssp, data = siCG, order.groups = FALSE)

#plot.new()
# This is the title page
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, ffile)
text(5, 7, sub)
text(5, 6, "QC of library by qPCR")

# This is the plot
plot(trgd, main=ffile, sub=sub, aspect = "fill", col = col,
        par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table= TRUE, 
	pch=19)
#	, ylim=range )


legend("bottomright",legend=legend, fill=col, border=TRUE, bty="o", y.intersp = 1.0, cex=0.7)


gain<-table(siCG$sidf[siCG$RQ>2])
loss<-table(siCG$sidf[siCG$RQ<0.5])
unch<-table(siCG$sidf[siCG$RQ<2 | siCG$RQ > 0.5])

# space out the positions of the table numbers on the x axis
sam<-length(gain)+2

text(y=3,x=6,"Summary Statistics from this run", col='black',cex=1.0)

text(y=2.5,x=6,"CG factor siRNAs in this expt", col='black',cex=1.0)

text(y=2,x=c(3:sam),labels=names(table(siCG$sidf)), col='black',cex=0.5)

text(y=1.5,x=1,labels="CGs UP (RQ>2)", col='dark green',cex=0.8)
text(y=1.5,x=c(3:sam),labels=as.numeric(gain), col='dark green',cex=0.8)

text(y=1,x=1,labels="CGs UNCHANGED", col='black',cex=0.8)
text(y=1,x=c(3:sam),labels=as.numeric(unch), col='black',cex=0.8)

text(y=0.5,x=1,labels="CGs DOWN (RQ<0.5)", col='dark red',cex=0.8)
text(y=0.5,x=c(3:sam),labels=as.numeric(loss), col='dark red',cex=0.8)


dev.off()




#############


# Summary data


# end
