##################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels for ddPCR probe design			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jan 2017) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); 
#source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); 
#biocLite("BSgenome.Hsapiens.UCSC.hg19"); 

library('BSgenome.Hsapiens.UCSC.hg19')
library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)

#################################################
# Directory structure - uncomment for first running of script
Project="ctDNA"
runID="AT094"
dir="validation"
homebase="/home/dyap/Projects/"
wdir=paste(homebase,Project,runID,dir,sep="/")

setwd(wdir)
getwd()

###############################################
# Save input files under $homebase/validation #
###############################################

# Check if the hg19 library is installed, if not install it and load it
genome="BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP.20120608)

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
SNP <-   installed.SNPs()
# If dbSNP is not currently installed load the latest version from source
if (dbSNP == SNP) print("Latest SNPS loaded") 
#	source("http://www.bioconductor.org/biocLite.R"); biocLite(dbSNP)

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

##############################################
######            User defined variables               ######
# Directory and file references
basedir=paste(homebase,Project,sep="/")
sourcedir=wdir
outpath=wdir

######################
# These are the input files
snvfile="expt1.csv"
#Sample,Chr,Pos,REF,ALT,Gene,Mut
#BOB047,17,7579455,C,T,TP53,A78T
#BOB049,3,178936094,C,A,PIK3CA,Q546K

#######################################
# This is the name of the design file - change it 
p3file=paste("SNV_design","txt",sep=".")

###############################################
file1 = paste(outpath,"SNV_Annotate.csv",sep="/")

############################
outfile=paste(wdir,p3file,sep="/")
input=paste(wdir,snvfile,sep="/")

# offsets (sequences on either side of SNV,indel for design space)
snvoffset=100
#indeloffset=200
WToffset=1

# Select the appropriate Genome (mask) - for reference only
#BSg="Hsapiens" # normal-reference
BSg="SNP_Hsapiens" # SNP-hard masked genome

##############################################

indf <- read.csv(file=input,  stringsAsFactors = FALSE, header=FALSE)

# This is how the data is found in the dataset
names(indf)<-c("SampleID","Chr","Pos","REF","ALT","Gene","Mut")

                    
outdf <- data.frame( ID = rep("", nrow(indf)),
                     Seq = rep("", nrow(indf)),
#                     Check = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)
                     
# Change this                     
offset <- snvoffset
                     
for (ri in seq(nrow(indf))) {
	
	Sample <- indf$SampleID[ri]
    	Chr <- paste("chr", indf$Chr[ri], sep="")
  	Pos <- as.numeric(indf$Pos[ri])
  	Ref <- indf$REF[ri]
  	Alt <- indf$ALT[ri]
  	Gene <- indf$Gene[ri]
  	Mut <- indf$Mut[ri]
	
# The indel or SNV has to match exactly so we get the reference sequences from hg19
 idseq <- as.character(getSeq(Hsapiens,Chr,Pos,Pos))
 
if (idseq == Ref) print("REF match")
if (!idseq == Ref) print("REF DOES NOT match")

 # This design space comprises of upstream and downstream annotated sequence and the exact reference seq
 # for hg19 for the SNV position and/or indel (This is important for matching)
  dseq <- 
paste(paste(as.character(getSeq(Hsapiens,Chr,Pos-offset,Pos-1)),"[",idseq,"/",Alt,"]",sep=""),
              as.character(getSeq(Hsapiens,Chr,Pos+1,Pos+offset)),sep="")

#check <-  as.character(getSeq(Hsapiens,Chr,Pos-5,Pos+5))
              
  outdf$ID[ri] <- paste(Sample,Gene,Mut,Chr,Pos,sep="_")          
  outdf$Seq[ri] <- dseq
#  outdf$Check[ri] <- check
  }

# Output file ID, chr, start, end, indel sequence, context for matching, design space 
# seq for probe design
write.csv(outdf, file = outfile)





