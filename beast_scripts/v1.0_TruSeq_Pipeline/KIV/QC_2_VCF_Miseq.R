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

# Tumour Evol / XenoDrug
direc=paste("/home/dyap/Projects/ctDNA",RUNID, sep="/")

setwd(direc)

##############################


# Select type of variants
expt="all"
#expt="selected"
mut="somatic"
#mut="germline"
# Do not think that "selected" is useful here

# Select type of QC - only one type
# type="freq"
# type="reads"
# type="calls"
type="strandbias"

if ( expt == "all" && type == "reads") append="allreads"
if ( expt == "all" && type == "freq") append="allfreq"
if ( mut == "germline" && type == "strandbias") append="gl_strandbias"
if ( expt == "selected" && type == "reads") append="reads"
if ( expt == "selected" && type == "freq") append="freq"
if ( mut == "somatic" && type == "strandbias") append="som_strandbias"

csvfile=paste(paste(name,append,sep="-"),"csv",sep=".")
pdffile=paste(paste(name,append,sep="-"),"pdf",sep=".")
namefile=paste(paste(paste(name,RUNID, sep="-"), append, sep="_"),"csv",sep=".")

graphtitle=paste(expt, paste(paste(paste(type, name, sep=" "), "Run-ID", sep=" "), RUNID, sep="-"), sep=" ")

#######################################


# change to either all or selected vcf or svcf
# read txt files with names of the form *.vcf (all variants called my BWA/GATK - MSR)
# file_names = list.files(pattern = '*[0-9].vcf');
# read txt files with names of the form *.svcf (only variants that we look for)
# file_names = list.files(pattern = '*[0-9].svcf');
# if ( expt == "all") pat =paste(name, "*.*[0-9].vcf", sep="")
# if ( expt == "selected") pat =paste(name, "*.*[0-9].svcf", sep="")
if ( mut == "germline") pat ="*germline"
if ( mut == "somatic") pat ="*somatic"

#file_names = list.files(pattern = pat)
file_names = list.files(pattern = pat)

# Extract all the VCFs into a concatenated VCF list
vcf_list = lapply(file_names, readVcf, "hg19", sep = "\t")

# Get the number of samples
# check with "list(vcf_list)"
samples <- length(vcf_list)


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
			sid <- colnames(vcf_list[[rj]])
			varn <- nrow(vcf_list[[rj]])

sumdf$Samples_ID[rj] <- sid
sumdf$Variants[rj] <- varn

			
# For each of the list of samples
len <- nrow(vcf_list[[rj]])

	if ( len > 0 ){
		

d.frame <- data.frame(	     ID = rep("position", len),
			     framename = rep (0, len),
			     stringsAsFactors = FALSE)

if (type == "calls") names(d.frame)[2] <- "Calls"
if (type == "freq") names(d.frame)[2] <- "Varalfreq"
if (type == "reads") names(d.frame)[2] <- "Seqdepth"
if (type == "strandbias") names(d.frame)[2] <- "StrandBias"

for (ri in seq(len) ) 	{

	d.frame$ID[ri] <- rownames(vcf_list[[rj]][ri])
if (type == "calls")  d.frame$Calls[ri] <- geno(vcf_list[[rj]][ri])$GT
if (type == "freq")  d.frame$Varalfreq[ri] <- geno(vcf_list[[rj]][ri])$VF
if (type == "reads") d.frame$Seqdepth[ri] <- geno(vcf_list[[rj]][ri])$DP
if (type == "strandbias") d.frame$StrandBias[ri] <- geno(vcf_list[[rj]][ri])$SB

			} 

# unique freq val or seq depth for summary plot
# must be the same length as names(d.frame)

names(d.frame)[2] <- paste(sid)


			} else { d.frame <- data.frame (ID ="pos") };

assign(paste("Nuclei", rj, sep=""), d.frame)

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

name.frame$filename[]<-file_names[]
name.frame$Sample_ID<-names(sum1)[2:(samples+1)]
sumdf

# Check list(vcf_list)
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

# Removing duplicates
sum2<-sum1
sum1<-sum2[!duplicated(lapply(sum2, summary))]

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

