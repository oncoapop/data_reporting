# File to print the CG Panel out in sequence to the samples
# To visualize the output in graphical format

# based on the script my Steve McKinney (Nov 2015)
# modified not to check primer matching
# by Dr Damian Yap Jan 2016
# Modified for MiSeq output Feb 2016

# Load libraries
library(foreign)
library(RColorBrewer)
library(lattice)
library(colorspace)
require(Heatplus)
# library(limma)
library(ggplot2)
require("nlme")

projdir <- "/home/dyap/Projects/Takeda_T3/CG"
wd="/home/dyap/Projects/Takeda_T3/CG_factors/screen"
file="KD_CG_factors_readthrough_filtered_ctrl_norm.tsv"

# This is the list of ordered primers for checking
pri <- paste(projdir, "/CG_primers_ordered.txt", sep = "")
pridf <- read.table(file = pri, header = FALSE, stringsAsFactors = FALSE)
names(pridf)[1] <- "Primers"
head(pridf)

# Typical output from Tyler's scripts
inputfile=paste(wd,file,sep="/")
QC<-read.table(inputfile, header=TRUE, blank.lines.skip=TRUE, sep="\t", na.strings="NA")

ffile=strsplit(file,split="[.]")[[1]][1]

#Instead get names from the dataset itself (MiSeq analysed data)
nnames<-nrow(table(QC$sample_id))
names(table(QC$sample_id))
samples<-names(table(QC$sample_id))[1:nnames]
leg<-paste0(samples, collapse=", ")

# This gets the samples from the data itself
# samples <- unique(ppdf$sample_id)
# sdf <- data.frame(Samples = samples, stringsAsFactors = FALSE)

### Ordering of primer gene names in output is arbitrary - need to make Primers variable same for either order of gene names.
QC$bk1 <- paste(QC$gene1, "@", QC$breakpt1, sep = "")
QC$bk2 <- paste(QC$gene2, "@", QC$breakpt2, sep = "")

g1priml <- lapply( QC$bk1, function(x) which( grepl(x, pridf$Primers) ) )
g2priml <- lapply( QC$bk2, function(x) which( grepl(x, pridf$Primers) ) )
intprim <- sapply( seq(nrow(QC) ), function(x) c(intersect(g1priml[[x]], g2priml[[x]] ), NA_integer_)[1] )
if (length(intprim) == nrow(QC)) { QC$PrimNo <- intprim } else { stop("Primer match error") }
QC$PrimersN <- sapply( seq(nrow(QC) ), function(x) c(intersect(g1priml[[x]], g2priml[[x]] ), NA_integer_)[1] )
QC$Primers <-  pridf$Primers[QC$PrimersN]

QC$g1prim <- sapply(QC$bk1, function(x) which(grepl(x, pridf$Primers)))
QC$g2prim <- sapply(QC$bk2, function(x) which(grepl(x, pridf$Primers)))

### ?? Need to make dataframe for all samples and all primer pairs and merge with pipeline data

prsadf <- expand.grid(as.character(samples), as.character(pridf$Primers), stringsAsFactors = FALSE)
names(prsadf) <- c("sample_id", "Primers")

QC$g1nm <- NULL
QC$g2nm <- NULL
QC$g1prim <- NULL
QC$g2prim <- NULL
QC$PrimersN <- NULL

QC$sidf <- factor(QC$sample_id)
QC$sidn <- as.numeric(QC$sidf)
QC$Primerssp <- gsub(":", "\n", QC$Primers)

#trgd <- groupedData(QC$ACTB_norm ~ sidn | Primerssp, data = QC, order.groups = FALSE)

# legend checking
len<-length(samples)
for (i in 1:len ) {
	print(i)
	# \\b marks the boundary so '\\b1\\b' returns 1 and not 10 or 31
	j<-as.character(unique(QC$sidf[grep(pattern=(paste(paste('\\b',i,sep=""),'\\b',sep="")),QC$sidn)]))
	print (j)
	}

# Subsetting the data for analysis
sub="siCG-factors_Plate2,3-AGBDG"

unsort<-QC

unsort$sidf <- factor(unsort$sample_id)
unsort$sidn <- as.numeric(unsort$sidf)

siCG <- unsort[with(unsort, order(sidn)), ]
ncol<-unique(siCG$sidn)

# Try to test to sort of significant deviations from siNT (which is already computated)
# RQ type analysis which does not work as if the CG is not detected in siNT then it is zero!



# Per sample correction of data
siCG$sidf <- as.character(siCG$sidf)
siCG$sidf[ siCG$sidf== "siNT-Plate-3"] <- "siNT-Plate3"

newdf<-cbind(siCG, read.table(text=as.character(siCG$sidf),sep="-"))
names(newdf)[18]<-"Gene"


# compute index of ordered 'median expression' and reassign          
#oind <- order(as.numeric(by(newdf$ACTB_norm, newdf$Gene, median)))    
#newdf$Gene <- ordered(newdf$Gene, levels=levels(newdf$Gene)[oind])
#boxplot(newdf$ACTB_norm ~ newdf$Gene, outline=FALSE, las=2, cex=0.8)
#siNT <- subset(newdf, Gene=="siNT", select=c(ACTB_norm,TUBA1B_norm,GAPDH_norm, Gene))
#simean <-mean(siNT$ACTB_norm)
#simed <-median(siNT$ACTB_norm)
#nsam <- length(ncol)-1
#lines(x=c(0:nsam),y=c(rep(simed,length(ncol))))

# Assign colour by gene
nclr=length(unique(newdf$Gene))
col.spectral <- colorRampPalette(brewer.pal(8,"Dark2"))(nclr)
col.scheme <- data.frame( Gene = rep("", nclr),
                        Col = rep("", nclr),
                        stringsAsFactors = FALSE)
for (c in seq(length(unique(newdf$Gene))))
	{
	col.scheme$Gene[c]<-as.character(unique(newdf$Gene)[c])
	col.scheme$Col[c]<-col.spectral[c]
	}

col.scheme

for (j in seq(nrow(newdf)))
        {
        match <- as.character(newdf$Gene[j])
        test <- col.scheme[col.scheme[,1] %in% c(match),2]
        if (length(test) != 0)
                {
                newdf$Col[j] <- test
                }

        }

# without ordering
plot(newdf$ACTB_norm ~ newdf$sample_id, outline=FALSE, las=2, cex.axis=0.6,
         ylab="Normalized expr (ACTB)",main=sub, col=col.scheme$Col)
legend("top",horiz = FALSE, legend=col.scheme$Gene, fill=col.scheme$Col, border=TRUE,
	bty="o", x.intersp = 0.2, cex=0.8, ncol=5)
########################################
pdffile=paste(sub,"pdf",sep=".")
pdf(file = pdffile, width = 8, height = 10)


# compute index of ordered 'media expression' and reassign          
oind <- order(as.numeric(by(newdf$ACTB_norm, newdf$sample_id, median)))    
newdf$sample_id <- ordered(newdf$sample_id, levels=levels(newdf$sample_id)[oind])
boxplot(newdf$ACTB_norm ~ newdf$sample_id, outline=FALSE, las=2, cex.axis=0.6,
         ylab="Normalized expr (ACTB)",main=sub, boxfill=newdf$Col, notch=FALSE, range = 0)

siNT2 <- subset(newdf, sample_id=="siNT-Plate2", select=c(ACTB_norm,TUBA1B_norm,GAPDH_norm, Gene))
siNT3 <- subset(newdf, sample_id=="siNT-Plate-3", select=c(ACTB_norm,TUBA1B_norm,GAPDH_norm, Gene))
simean2 <-mean(siNT2$ACTB_norm)
simean3 <-mean(siNT3$ACTB_norm)
simed2 <-median(siNT2$ACTB_norm)
simed3 <-median(siNT3$ACTB_norm)
nsam <- length(ncol)
nsi <-length(ncol)+1
lines(x=c(0:nsam),y=c(rep(simed2,nsi)))
lines(x=c(0:nsam),y=c(rep(simed3,nsi)))

dev.off()



# Try plotting a different way
trgd <- groupedData(newdf$ACTB_norm ~ sidn | Primerssp , data = newdf, order.groups = FALSE)
plot(trgd, main=ffile, sub=leg, aspect = "fill",
        par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5),
        as.table = TRUE)





## check ####
siSRSF9 <- subset(newdf, sample_id=="siSRSF9-2", select=c(ACTB_norm,TUBA1B_norm,GAPDH_norm, Gene))
##> siSRSF9
##[1] ACTB_norm   TUBA1B_norm GAPDH_norm  Gene
##<0 rows> (or 0-length row.names)
# QC qPCR failed

#########
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
logRQ<-log(QC$ACTB_norm)
trgd <- groupedData(logRQ ~ sidn | Primerssp, data = QC, order.groups = 
	FALSE)
plot(trgd, main=ffile, sub=leg, aspect = "fill", 
	par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table 
	= TRUE)
dev.off()


# end
