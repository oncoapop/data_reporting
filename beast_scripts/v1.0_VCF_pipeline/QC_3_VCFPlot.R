
###########################################################
#  QC for HeatPlots from MiSeq reporter processed ouputs  #
#       Modified by Dr Damian Yap, 15 Nov 2016	          #
###########################################################


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

#####################################################################################################
# To run this script change the setwd()
setwd("/home/dyap/Projects/Tumour_Evol/Xenodrug")
sample="S"
ID="SA535"
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
xaxislab=paste(ID, "Sample", sep=" ")

# Read in all positions
#htert<-read.table(file="htert_pos.txt")
#htert$V3 <-"htert"
#hct116<-read.table(file="hct116_pos.txt")
#hct116$V3<-"hct116"
#shared<-read.table(file="shared_pos.txt")
#shared$V3 <- "shared"

# Read in Sample Names and Mixing proportions
#prop<-read.table(file="MiSeq samples 18 Jul.txt" ,sep="\t", header=TRUE)
prop<-read.table(file=infile ,sep=",", header=TRUE)


####################################################################################################

sum1 <- read.csv(file=sfile)

###################################
## XXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
## Processing position information

# combine into one file
#hh<-rbind(hct116,htert)
#all<-rbind(shared,hh)


# Select n colors
#n <- length(table(all$V3))
#rowcol<-brewer.pal(n+1, "Accent")

#posinfo <- data.frame(	ID = rep("", nrow(all)),
#			Type = rep("", nrow(all)),
#			Col = rep("", nrow(all)),
#			stringsAsFactors = FALSE)

#for (a in seq(nrow(all)))
#	{
#	chr <- all$V1[a]
#	pos <- all$V2[a]
#	type <- all$V3[a]

#	posinfo$ID[a] <- paste(chr, pos, sep="_")
#	posinfo$Type[a] <- type
	
#	for ( b in seq(n) )
#		{
#		if (type == names(table(all$V3)[b])) color=rowcol[b]
#		}

#	posinfo$Col[a] <- color

#	}


# Legend info posinfo$Types and posinfo$Col
#legend=rownames(table(posinfo[2:3]))
#legend[n+1]="unclassified"

#fill=colnames(table(posinfo[2:3]))
#fill[n+1]=rowcol[n+1]

## XXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
#####################################
# Preparing the data for clustering (removing NAs)

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])

# col headers - sample names
names(sum1)

# Filters out all the positions that failed in all samples
indiv=length(sum1)-1
filt <-rownames(ef[rowSums(is.na(ef))==indiv,])
filt

# Label rownames with ID
rownames(ef) <- sum1$ID

colnames(ef)
ff<-ef[,order(as.numeric(colnames(ef)))]

# If filt=NULL then all primers work so NA = zero
ff[is.na(ff)] <- 0

##############################################
## XXXXXXXXXXXXXXXXXXXXXXXXXXXXX #

# Label according to samplesheet
# Select n colors (first 2 colors from row Col ie cell lines)

#col2<-rowcol[length(prop)]

#colcols <- rev(colorRampPalette(brewer.pal(11,"PRGn"))(1000))

#propinfo <- data.frame(	Sample = rep("", nrow(prop)),
#			Tertp = rep("", nrow(prop)),
#			HCTp = rep("", nrow(prop)),
#			Col = rep("", nrow(prop)),
#			stringsAsFactors = FALSE)
			

#for (a in seq(nrow(prop)))
#	{
#	propinfo$Sample[a] <- prop$Sample[a]
#	propinfo$Tertp[a] <- prop$hTert[a]
#	propinfo$HCTp[a] <- prop$HCT116[a]

#	if (prop$hTert[a]*1000 < 1) 
#		{
#		propinfo$Col[a] <- colcols[1]
#		} else {
#			propinfo$Col[a] <- colcols[prop$hTert[a]*1000]
#			}
#	}



# Colside color matrix
#csc <- data.frame(	Sam = rep("", ncol(ff)),
#			Col = rep("", ncol(ff)),
#			stringsAsFactors = FALSE)

# Get the colors for different sample proportions
#for (j in seq(ncol(ff)))
#	{
#	match <- colnames(ff)[j]
#	csc$Sam[j]<- match
#	test <- propinfo[propinfo[,1] %in% c(match),4]
#	if (length(test) != 0) 
#		{
#		csc$Col[j] <- test
#		} 
#		
#	}
	
# Rowside color matrix
#rsc <- data.frame(	ID = rep("", nrow(ff)),
#			Col = rep("", nrow(ff)),
#			stringsAsFactors = FALSE)
#
# Get the colors for different types of positions
#for (j in seq(nrow(ff)))
#	{
#	match <- rownames(ff)[j]
#	rsc$ID[j]<- match
#	test <- posinfo[posinfo[,1] %in% c(match),3]
#	if (length(test) != 0) 
#		{
#		rsc$Col[j] <- test
#		} else {
#		# Unclassified positions
#		rsc$Col[j] <- rowcol[n+1]
#		}
#		
#	}

## XXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
##################################

# Plot Colours
# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", margins=c(5,10))
# hmcols<-colorRampPalette(c("dark green","red"))(100)
hmcols <- colorRampPalette(brewer.pal(11,"Spectral"))(100)
# display.brewer.all()
################################################### HEatMaps2

title=paste("Heatmap of",sep=" ",filename)          

xaxislab2=""    

pdf(pdffile3, width=7, height=8)

#heatmap.2(ff, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
#, cexCol=0.8, cexRow=0.6, col = hmcols, RowSideColors=rsc$Col, ColSideColors=csc$Col, 
#trace="none")

heatmap.2(ff, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.8, cexRow=0.6, col = hmcols, trace="none")

legend("topright",legend=legend, fill=fill, border=TRUE, bty="o", y.intersp = 0.7, cex=0.7)

dev.off()



###########################UNDER DEV#################################



##############################PLOT SELECTED VARIANTS ONLY ##########################

# New matrix gf (selected positions)

gf <- ff[rownames(ff) %in% posinfo$ID,]

# Rowside color matrix
rsc <- data.frame(	ID = rep("", nrow(gf)),
			Col = rep("", nrow(gf)),
			stringsAsFactors = FALSE)

# Get the colors for different types of positions
for (j in seq(nrow(gf)))
	{
	match <- rownames(gf)[j]
	rsc$ID[j]<- match
	test <- posinfo[posinfo[,1] %in% c(match),3]
	if (length(test) != 0) 
		{
		rsc$Col[j] <- test
		} else {
		# Unclassified positions
		rsc$Col[j] <- rowcol[n+1]
		}
		
	}

# Plot Colours
# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", margins=c(5,10))
# hmcols<-colorRampPalette(c("dark green","red"))(100)
hmcols <- colorRampPalette(brewer.pal(11,"Spectral"))(100)
# display.brewer.all()
################################################### HEatMaps2

title=paste("Heatmap of",sep=" ",filename)

xaxislab2=""

pdf(pdffile2, width=7, height=8)

heatmap.2(gf, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.8, cexRow=0.6, col = hmcols, ColSideColors=csc$Col, RowSideColors=rsc$Col, 
trace="none")

legend("topright",legend=legend, fill=fill, border=TRUE, bty="o", y.intersp = 0.7, cex=0.7)

dev.off()

#########################END SELECTED VARIANTS ####################################
########################################################## levelplot

title="Level Plot QC of HCT116-hTert mixing expt"

pdf(pdffile, width=6, height=6)
levelplot(ff, main=title, xlab=xaxislab, ylab="Position", aspect="fill", cexCol=0.8, 
col.regions = hmcols)
dev.off()


pdf(pdffile, width=6, height=6)
reg2 = regHeatmap(ff, legend=2,breaks=-4:4)
plot(reg2)
dev.off()


########################################################


