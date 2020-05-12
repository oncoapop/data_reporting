# QC for MiSeq VCF files 
# QC For Single Cell Project data from VCF

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

#####################################################################################################
# To run this script change the sample and directory names

basedir="/share/lustre/backup/dyap/Projects"
sample="HCT116"
run="AMU63"
#check="depth"
check = "freq"
#check = "calls"

workdir=paste(paste(basedir, sample, sep="/"),run,sep="/")
setwd(workdir)

# Outputs
dir=paste(basedir, sample, sep="/")

exptname=paste(sample,run,sep="_")
if (check == "depth") filename=paste(exptname,"reads",sep="_")
if (check == "freq") filename=paste(exptname,"freq",sep="_")
if (check == "calls") filename=paste(exptname,"calls",sep="_")

output=paste(dir,filename,sep="/")

csvfile=paste(output, "csv", sep=".")
pdffile=paste(output, "pdf", sep=".")
pdffile2=paste(output, "cluster.pdf", sep="")
pdffile3=paste(output, "cluster2.pdf", sep="")
vennfile=paste(output, "Venn.pdf", sep="-")

title=filename
xaxislab=paste(sample, "Sample", sep=" ")

# Read in all positions and whether they are germline or somatic
# For the bed file that Jas produced position is start+1 = SNV position
#posfile = paste(workdir, "SA029_SNV_typing.txt", sep="/")
#postype <- read.table(file=posfile)
#postype$V5 <- paste(gsub("chr","",postype$V1), postype$V2+1, sep="_")


########################################################################################

vcfdir=paste(workdir, run, sep="/")
setwd(vcfdir)

# Read Samplesheet in
# Change the skip lines to ensure it is read correctly in

line<-system("grep -n 'Data' SampleSheet.csv | awk -F':' '{print $1}'", intern=TRUE)
samplesheet<-read.csv("SampleSheet.csv", skip=line, header=TRUE)


####################################################################################################


# read txt files with names of data files .vcf
# pat=paste(sample,"*.*vcf",sep="")
pat=paste("","*.*vcf",sep="")
file_names = list.files(pattern = pat);
file_names

# Extract all the VCFs into a concatenated VCF list
vcf_list = lapply(file_names, readVcf, "hg19", sep = "\t")

# Change the number of samples
# check with "list(vcf_list)"
samples <- length(vcf_list)

createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)


#################################################
# THIS IS THE CURRENT MODULE THAT WORKS !!!!!   #
#################################################



sumdf <- data.frame(	Sample_ID = rep("", samples),
			Variants = rep(0, samples),
			Rows = rep("", samples),
			stringsAsFactors = FALSE)

for (rj in seq(samples)) 	{
			sid <- colnames(vcf_list[[rj]])
			varn <- nrow(vcf_list[[rj]])

sumdf$Samples_ID[rj] <- sid
sumdf$Variants[rj] <- varn


# For each of the list of samples
len <- nrow(vcf_list[[rj]])

	if ( len > 0 ){


d.frame <- data.frame(	     ID = rep("", len),
 			     framename = rep (0, len),
			     stringsAsFactors = FALSE)

if (check == "freq" ) names(d.frame)[2] <- "Varalfreq"
if (check == "depth" ) names(d.frame)[2] <- "Seqdepth"

for (ri in seq(len) ) 	{
# This extracts the postion information instead of the rs name for dbSNP positions (chr_position)
d.frame$ID[ri] <- paste(gsub("chr","",seqnames(vcf_list[[rj]][ri])), as.character(start(vcf_list[[rj]][ri])),sep="_")
#	d.frame$ID[ri] <- rownames(vcf_list[[rj]][ri])
if (check == "calls")	d.frame$Calls[ri] <- geno(vcf_list[[rj]][ri])$GT
if (check == "freq")	d.frame$Varalfreq[ri] <- geno(vcf_list[[rj]][ri])$VF
if (check == "depth")	d.frame$Seqdepth[ri] <- geno(vcf_list[[rj]][ri])$DP

			} 

# unique freq val or seq depth for summary plot
# must be the same length as names(d.frame)
names(d.frame)[2] <- paste(sid)


			} else { d.frame <- "NULL" };

assign(paste("Nuclei", rj, sep=""), d.frame)

# for first value 
	if ( rj == 1 && d.frame != "NULL"  ) {
		first <- d.frame }
 
# combining successive data.frames
	if ( rj == 2 && d.frame != "NULL"  ) { 
		sum1 <- merge(first, d.frame, by="ID", all=TRUE) }

# combining successive data.frames
	if ( rj > 2 && d.frame != "NULL"  ) { 
		a <- count(1)
		sum1 <- merge(sum1, d.frame, by="ID", all=TRUE) 
		} else { print("skip") }

}

write.table(sum1,file=csvfile,sep=",",row.names=FALSE,col.names=TRUE)

###################################
## Processing position information

# Colors rows in terms of whether pos is germline or somatic
# Select n colors
# n <- length(table(postype$V4))
# rowcol<-brewer.pal(n+1, "Accent")
# types = names(table(postype$V4))

#posinfo <- data.frame(	ID = rep("", nrow(postype)),
#			Type = rep("", nrow(postype)),
#			Col = rep("", nrow(postype)),
#			stringsAsFactors = FALSE)

#for (a in seq(nrow(postype)))
#	{
#	posinfo$ID[a] <- as.character(postype$V5[a])
#	posinfo$Type[a] <- as.character(postype$V4[a])
	
#	for ( b in seq(n) )
#		{
#		if (posinfo$Type[a] == types[b]) color=rowcol[b]
#		}

#	posinfo$Col[a] <- color

#	}


# Legend info posinfo$Types and posinfo$Col
# legend=rownames(table(posinfo[2:3]))
# legend[n+1]="unclassified"

#fill=colnames(table(posinfo[2:3]))
#fill[n+1]=rowcol[n+1]

#####################################
# Preparing the data for clustering (removing NAs)

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])

# col headers - unique nuclei
names(sum1)

# Filters out all the positions that failed in 80% samples
indiv=floor(length(sum1)-(0.2*ncol(sum1)))
filt <-rownames(ef[rowSums(is.na(ef))==indiv,])
filt

# Label rownames with ID
rownames(ef) <- sum1$ID

# Label colnames (sample ID with Sample Names from Samplesheet)

for (j in seq(ncol(ef)))
	{
	match <- colnames(ef)[j]
	test <- as.character(samplesheet[samplesheet[,1] %in% c(match),2])
	if (length(test) != 0) 
		{
		colnames(ef)[j] <- test
		} 
		
	}

colnames(ef)
ff<-ef[,order(as.numeric(colnames(ef)))]

#If filt=NULL then all primers work so NA = zero
ff[is.na(ff)] <- 0
#This condition is not true for single cell experiments!

##############################################

# Label according to samplesheet
# Select n colors (first 2 colors from row Col ie cell lines)

#col2<-rowcol[1:2]

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
rsc <- data.frame(	ID = rep("", nrow(ff)),
			Col = rep("", nrow(ff)),
			stringsAsFactors = FALSE)

# Get the colors for different types of positions
for (j in seq(nrow(ef)))
	{
	match <- rownames(ef)[j]
	rsc$ID[j]<- match
	test <- postype[postype[,5] %in% c(match),4]
	if (length(test) != 0) 
		{
		rsc$Col[j] <- rowcol[test]
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

title=paste(sample, "Identification Experiments - all variants", sep=" ")

xaxislab2=paste("Samples from Run", run, sep=" ")

pdf(pdffile3, width=7, height=8)

heatmap.2(ff, main=title, xlab=xaxislab2, ylab="Positions", scale="none", key = TRUE
, cexCol=0.8, cexRow=0.6, col = hmcols, srtCol=20, trace="none")

legend("topright",legend=legend, fill=fill, border=TRUE, bty="o", y.intersp = 0.7, cex=0.7)

dev.off()



###########################UNDER DEV#################################



##############################PLOT SELECTED VARIANTS ONLY ##########################

# New matrix gf (selected positions)

gf <- ef[rownames(ef) %in% postype$V5,]

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

title="HCT116-hTert Cell Mixing Expt (Sel Var Freq)"

xaxislab2=paste("Samples from Run", run, sep=" ")

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


