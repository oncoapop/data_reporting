# QC for TruSeq custom and Cancer Panel run based on the 
# QC For Standard Targeted Sequencing Output on Miseq
# standard vcf output from miSeq reporter
# selected vcfs processed from AmpliconManifest using the script scvf_generator.sh

# Updated for QC Nov 2015 and Nov 2016
# Damian Yap, PhD

################################################################################
####        WSOP2013-001 Target Validation DeepSeq on MiSeq                 ####
####           SOP2013-002:  Post-Miseq Run Procedure SOP                   ####
####        	   Modified for custom amplicon manifests	    	    ####
####                written by Dr D Yap, 22 May 2013                        ####
################################################################################


# This script requires R  3.0.1
# reads VCF 4.1 formats as produced by MiSeq Reporter
# https://samtools.github.io/hts-specs/VCFv4.1.pdf

#source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
#biocLite(c("GenomicFeatures", "AnnotationDbi"))

library("VariantAnnotation")
library("IRanges")
library("GenomicRanges")
library(foreign)
library(lattice)
library("RColorBrewer")

#####


#############################
# For QC of Run files

# For Pipeline QC
# Name of Single Cell Sample
name="Expt1"
RUNID="AT094"

# Variants from the TruSeq pipeline (pre-processed)
direc=paste(paste("/home/dyap/Projects/ctDNA",RUNID, sep="/"),"variants",sep="/")

setwd(direc)
getwd()

##############################


# Select type of variants
expt="all"
#expt="selected"
# mut="WTonly"
  mut="TUMOURonly"
# Do not think that "selected" is useful here

# Modification to plot all relevant values
append="summary"
type=append

csvfile=paste(paste(name,append,sep="-"),"csv",sep=".")
pdffile=paste(paste(name,append,sep="-"),"pdf",sep=".")
namefile=paste(paste(paste(name,RUNID, sep="-"), append, sep="_"),"csv",sep=".")

graphtitle=paste(expt, paste(paste(paste(type, name, sep=" "), "Run-ID", sep=" "), RUNID, sep="-"), sep=" ")

#######################################

if ( mut == "WTonly") pat ="WTonly*"
if ( mut == "TUMOURonly") pat ="TUMOURonly*"

#notemptys = list.files(pattern = pat)
file_names = list.files(pattern = pat)


info = file.info(file_names)
empty = rownames(info[info$size == 0, ])
notempty = rownames(info[info$size > 0, ])


##########
# Extract all the postions into a concatenated list
# all_list = lapply(file_names, read.csv, sep = ",")
ne_list = lapply(notempty, read.csv, sep = ",", header=FALSE)

# Get the number of samples
# check with "list(ne_list)"
samples <- length(ne_list)


createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)



#############################################################################
#     This is automatically changed when different Expt type is selected    #
#############################################################################


name.frame <- data.frame(	Sample_ID = rep("", samples),
			filename = rep("", samples),
			Description = rep("", samples),
			stringsAsFactors = FALSE)

sumdf <- data.frame(	Sample_ID = rep("", samples),
			Variants = rep(0, samples),
			Rows = rep("", samples),
			stringsAsFactors = FALSE)


##########################################################

for (rj in seq(samples)) 	{
			sid <- notempty[[rj]]
			varn <- nrow(ne_list[[rj]])

sumdf$Samples_ID[rj] <- sid
sumdf$Variants[rj] <- varn

			
# For each of the list of samples
len <- nrow(ne_list[[rj]])

	if ( len > 0 ){
		

d.frame <- data.frame(	     ID = rep("chr_position", len),
			     MUT = rep ("REF>ALT", len),
			     stringsAsFactors = FALSE)

for (ri in seq(len) ) 	{

	d.frame$ID[ri] <- paste(as.character(ne_list[[rj]]$V1)[ri],as.character(ne_list[[rj]]$V2)[ri],sep="_")
	d.frame$MUT[ri] <- paste(as.character(ne_list[[rj]]$V4)[ri],as.character(ne_list[[rj]]$V5)[ri],sep=">")

			} 

# unique freq val or seq depth for summary plot
# must be the same length as names(d.frame)

names(d.frame)[2] <- paste(sid)


			} else { d.frame <- data.frame (ID ="pos") };

assign(paste("Variant", rj, sep=""), d.frame)

		# if NULL frame (just add the header as placeholder)
	if ( d.frame == "NULL"  ) {

					d.frame$ID[ri] <- "X"
 					}

# for first value 
	if ( rj == 1 ) {
		first <- d.frame } 
 
# combining successive data.frames
	if ( rj == 2 ) { 
		sum1 <- merge(first, d.frame, by="ID", all=TRUE) }

# combining successive data.frames
	if ( rj > 2 ) { 
		a <- count(1)
		sum1 <- merge(sum1, d.frame, by="ID", all=TRUE) 
		} else { print("skip") }

}


################################### STOP HERE AND CHECK COLUMN NAMES FOR LABELLING OF CONTROLS

name.frame$filename[]<-notempty[]
name.frame$Sample_ID<-names(sum1)[2:(samples+1)]
sumdf

# Check list(ne_list)
#names(sum1)

# Labels for Main Run files
for (des in seq(samples)) 
	{
	point<-des+1
	names(sum1)[point] <- strsplit(name.frame$filename[des],split="_S")[[1]][1]
	name.frame$Description[des] <- strsplit(name.frame$filename[des],split="_S")[[1]][1]
	}

# Check that that autolabelling is correct
name.frame

# Removing duplicates (not required)
#sum2<-sum1
#sum1<-sum2[!duplicated(lapply(sum2, summary))]

########################################
# STOP HERE TO CHECK BEFORE PROCEEDING #
########################################

# Write checkfiles
write.table(name.frame,file=namefile,sep=",",row.names=FALSE,col.names=TRUE)

# Write csv
write.table(sum1,file=csvfile,sep=",",row.names=FALSE,col.names=TRUE)

##########################################################

# Drawing heatmap

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])

# col headers - unique nuclei
names(sum1)

# Label rownames with ID (these are positions)
rownames(ef) <- sum1$ID



########################################################

pdf(pdffile, width=10, height=10)

###################
#heatmap(ef, Rowv=NA, Colv=NA, na.rm=TRUE, main="Selected Variant Freq SA494 (Run ID A49AR)", xlab="SA494 Nuclei", ylab="Position", aspect="fill", cexCol=0.8, col=rev(heat.colors(1000)))
# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", aspect="fill", margins=c(5,10))
###################

rgb.palette <- colorRampPalette(c("dark green", "red"), space = "rgb")
# rgb.palette <- brewer.pal(10, "Spectral")


if (type == "freq") levelplot(ef, xlab="Variant", las=2, ylab="Samples", main=graphtitle, col.regions=rgb.palette(100), 
cuts=20, at=seq(0.005,1.0,0.05), aspect="fill", scales=list(x=list(rot=75, cex=0.4),y=list(cex=0.8)) )

if (type == "reads") levelplot(ef, xlab="Variant", las=2, ylab="Samples", main=graphtitle, col.regions=rgb.palette(200), 
cuts=200, at=seq(20,5020,200), aspect="fill", scales=list(x=list(rot=90, cex=0.8),y=list(cex=0.8)) )

if (type == "strandbias") levelplot(ef, xlab="Variant", las=2, ylab="Samples", main=graphtitle, col.regions=rgb.palette(100), 
cuts=10, at=seq(-99,0,10), aspect="fill", scales=list(x=list(rot=75, cex=0.3),y=list(cex=0.8)) )


dev.off()

