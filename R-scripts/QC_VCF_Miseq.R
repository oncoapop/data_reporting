
# QC For Standard Targeted Sequencing Output on Miseq
# standard vcf output from miSeq reporter
# selected vcfs processef from AmpliconManifest using the script scvf_generator.sh
# Updated for QC Nov 2015
# Damian Yap

################################################################################
####        WSOP2013-001 Target Validation DeepSeq on MiSeq                 ####
####           SOP2013-002:  Post-Miseq Run Procedure SOP                   ####
####        							    ####
####                written by Dr D Yap, 22 May 2013                        ####
################################################################################


# This script requires R  3.0.1
# 
#source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
#biocLite(c("GenomicFeatures", "AnnotationDbi"))


library("VariantAnnotation")
library("IRanges")
library("GenomicRanges")
library(foreign)
library(lattice)
#####


#############################
# For QC of Run files

# For Pipeline QC
# Name of Single Cell Sample
name="SA535"
RUNID="AGALJ&AGAJY"

# MOMAC 14
# RunID="130603_M00897_0033_000000000-A49AR"
# Alignment="Alignment"
# Aligndir=paste(paste("/Volumes/Monco/MiSeq Run Files/MiSeq Analysis Files",RunID,sep="/"), 
# direc=paste(Aligndir,Alignment,sep="/")

# Single Cell
#basedir="/home/dyap/Projects/Single_Cell/Analysis"
#run=paste(name,RUNID,sep="-")
#direc=paste(basedir, run, sep="/")

# Tumour Evol / Xeno
direc=paste("/home/dyap/Projects/Tumour_Evol", name, sep="/")

setwd(direc)

##############################


# Select type of variants
#expt="all"
expt="selected"
# Select type of QC - only one type
type="freq"
# type="reads"
# type="calls"

if ( expt == "all" && type == "reads") append="allreads"
if ( expt == "all" && type == "freq") append="allfreq"
if ( expt == "selected" && type == "reads") append="reads"
if ( expt == "selected" && type == "freq") append="freq"

csvfile=paste(paste(name,append,sep="-"),"csv",sep=".")
pdffile=paste(paste(name,append,sep="-"),"pdf",sep=".")
namefile=paste(paste(paste("Name", name, sep="-"), append, sep="-"),"csv",sep=".")

graphtitle=paste(expt, paste(paste(paste(type, name, sep=" "), "Run-ID", sep=" "), RUNID, sep="-"), sep=" ")

#######################################


# change to either all or selected vcf or svcf
# read txt files with names of the form *.vcf (all variants called my BWA/GATK - MSR)
# file_names = list.files(pattern = '*[0-9].vcf');
# read txt files with names of the form *.svcf (only variants that we look for)
# file_names = list.files(pattern = '*[0-9].svcf');
if ( expt == "all") pat =paste(name, "*.*[0-9].vcf", sep="")
if ( expt == "selected") pat =paste(name, "*.*[0-9].svcf", sep="")

#file_names = list.files(pattern = pat)
file_names = list.files(pattern = pat)

# Extract all the VCFs into a concatenated VCF list
vcf_list = lapply(file_names, readVcf, "hg19", sep = "\t")

# Get the number of samples
# check with "list(vcf_list)"
samples <- length(vcf_list)


createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)



#########################################################
#     Change this to output Var Freq or Seq Depth (no of reads)      #
#########################################################


name.frame <- data.frame(	Sample_ID = rep("", samples),
			filename = rep("", samples),
			Description = rep("", samples),
			stringsAsFactors = FALSE)

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
		

d.frame <- data.frame(	     ID = rep("position", len),
			     framename = rep (0, len),
#			     Seqdepth = rep(0, len),
			     stringsAsFactors = FALSE)

if (type == "calls") names(d.frame)[2] <- "Calls"
if (type == "freq") names(d.frame)[2] <- "Varalfreq"
if (type == "reads") names(d.frame)[2] <- "Seqdepth"

for (ri in seq(len) ) 	{

	d.frame$ID[ri] <- rownames(vcf_list[[rj]][ri])
if (type == "calls")  d.frame$Calls[ri] <- geno(vcf_list[[rj]][ri])$GT
if (type == "freq")  d.frame$Varalfreq[ri] <- geno(vcf_list[[rj]][ri])$VF
if (type == "reads") d.frame$Seqdepth[ri] <- geno(vcf_list[[rj]][ri])$DP

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
# 
write.table(obj_name,file="C:\\Users\\dyap_000\\Documents\\R\\SA494\\SA494.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)

##########################################################

# Drawing heatmap

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])

# col headers - unique nuclei
names(sum1)

# Label rownames with ID (these are positions)
rownames(ef) <- sum1$ID



########################################################

pdf(pdffile, width=10, height=4)

###################
#heatmap(ef, Rowv=NA, Colv=NA, na.rm=TRUE, main="Selected Variant Freq SA494 (Run ID A49AR)", xlab="SA494 Nuclei", ylab="Position", aspect="fill", cexCol=0.8, col=rev(heat.colors(1000)))
# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", aspect="fill", margins=c(5,10))
###################

rgb.palette <- colorRampPalette(c("dark green", "red"), space = "rgb")

if (type == "freq") levelplot(ef, xlab="Variant", las=2, ylab="Samples", main=graphtitle, col.regions=rgb.palette(10), 
cuts=10, at=seq(0.2,1.0,0.1), aspect="fill", scales=list(x=list(rot=90, cex=0.3),y=list(cex=0.5)) )

if (type == "reads") levelplot(ef, xlab="Variant", las=2, ylab="Samples", main=graphtitle, col.regions=rgb.palette(200), 
cuts=200, at=seq(20,5020,200), aspect="fill", scales=list(x=list(rot=90, cex=0.3),y=list(cex=0.5)) )


dev.off()

