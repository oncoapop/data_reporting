##################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels to design primers for 			
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jul 2013) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); 
source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); 
biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)

#################################################
# Directory structure - uncomment for first running of script
Project="TNBC"
homebase="/Users/dyap/Documents/Breast Cancer"
setwd(homebase)
# system('mkdir TNBC')
setwd(paste(homebase,Project,sep="/"))
# system('mkdir primer3')
# system('mkdir positions')
# system('mkdir Annotate')
getwd()

#######################################
# Save input files under $homebase/positions#
#######################################

# Check if the hg19 library is install, if not install it and load it
genome="BSgenome.Hsapiens.UCSC.hg19"
library('BSgenome.Hsapiens.UCSC.hg19')
library(SNPlocs.Hsapiens.dbSNP.20120608)

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
SNP <-   installed.SNPs()
# If dbSNP is not currently installed load the latest version from source
# if (!dbSNP %in%  SNP) source("http://www.bioconductor.org/biocLite.R"); 
biocLite(dbSNP)

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

##############################################
######            User defined variables               ######
# Directory and file references
basedir=paste(homebase,Project,sep="/")
sourcedir=paste(basedir,"positions", sep="/")
# outdir=paste(basedir,"positions", sep="/")
p3dir=paste(basedir,"primer3", sep="/")
outpath=paste(basedir,"Annotate", sep="/")

######################
# These are the input files
snvfile="primerIn-TNBC-SNV-fix-List.txt"
indelfile="primerIn-TNBC-indel-fix-List.txt"

#######################################
# This is the name of the design file - change it 
p3file=paste("SNV_design","txt",sep=".")

###############################################
file1 = paste(outpath,"SNV_Annotate.csv",sep="/")

############################
outfile=paste(p3dir,p3file,sep="/")
input=paste(sourcedir,snvfile,sep="/")

# offsets (sequences on either side of SNV,indel for design space)
snvoffset=200
indeloffset=200
WToffset=5

# Select the appropriate Genome (mask) - for reference only
#BSg="Hsapiens" # normal-reference
#BSg="SNP_Hsapiens" # SNP-hard masked genome

##############################################

indf <- read.table(file=input,  stringsAsFactors = FALSE, header=TRUE)

# This is how the data is found in the dataset
table(indf$samp)

#SA <- names(table(indf$samp))
# dataset <- subset(indf, indf$samp == SA)

                    
outdf <- data.frame(ID = rep("", nrow(indf)),
		     Chr = rep("", nrow(indf)),
                     Start = rep(0, nrow(indf)),
                     End = rep(0, nrow(indf)),
                     Indel = rep("", nrow(indf)),
                     Context = rep("", nrow(indf)),
                     Design = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)
                     
# Change this                     
offset <- snvoffset
                     
for (ri in seq(nrow(indf))) {
	
	id <- indf$ID[ri]
    	chr <- paste("chr", indf$chr[ri], sep="")
  	start <- as.numeric(indf$startPos[ri])
    	end <- as.numeric(indf$endPos[ri])

# The indel or SNV has to match exactly so we get the reference sequences from hg19
 idseq <- as.character(getSeq(Hsapiens,chr,start,end))
 
 # If the index <5 or it is an SNV then we need 5 bp on either side to match (for 
visualization only)
 if (nchar(idseq) < 5) cxt <- 
as.character(getSeq(Hsapiens,chr,start-WToffset,end+WToffset)) else cxt <- idseq
 
 # This design space comprises of upstream and downstream annotated sequence and the 
exact reference seq
 # for hg19 for the SNV position and/or indel (This is important for matching)
  dseq <- 
paste(paste(as.character(getSeq(SNP_Hsapiens,chr,start-offset,start-1)),idseq,sep=""),
              as.character(getSeq(SNP_Hsapiens,chr,end+1,end+offset)),sep="")
              
  outdf$ID[ri] <- id          
  outdf$Chr[ri] <- chr
  outdf$Start[ri] <- start
  outdf$End[ri] <- end
  outdf$Indel[ri] <- idseq
  outdf$Context[ri] <- cxt
  outdf$Design[ri] <- dseq
  }

# Output file ID, chr, start, end, indel sequence, context for matching, design space 
seq for primer3 design
write.csv(outdf, file = outfile)


######################

# For annotation files

andf1 <- data.frame(Chr1 = rep("", nrow(indf)),
                     Pos1 = rep(0, nrow(indf)),
                     Pos2 = rep(0, nrow(indf)),
                     WT1 = rep("", nrow(indf)),
                     SNV1 = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)
        
andf2 <- data.frame(Chr2 = rep("", nrow(indf)),
                     Pos3 = rep(0, nrow(indf)),
                     Pos4 = rep(0, nrow(indf)),
                     WT2 = rep("", nrow(indf)),
                     SNV2 = rep("", nrow(indf)),
                     stringsAsFactors = FALSE)       
                     
for (ri in seq(nrow(indf))) {
# For GetSeq we need the Chr prefix for chromosome but not for ANNOVAR  
# Assume indels are <40bp on the same chromosome
   id 		<- indf$ID[ri]
   chrom 	<- indf$chr[ri]
   chr 		<- paste("chr", indf$chr[ri], sep="")
   pos1 	<- as.numeric(indf$startPos[ri])
   pos2 	<- as.numeric(indf$endPos[ri])
  
    wt1 <- as.character(getSeq(SNP_Hsapiens,chr,pos1,pos2))
    
# Fake the SNV to be just the complement of WT position (as SNV allele is not known)

if (wt1=="A") snv1 <- "T"
if (wt1=="C") snv1 <- "G"
if (wt1=="G") snv1 <- "C"
if (wt1=="T") snv1 <- "A"
  
  andf1$Chr1[ri] <- chrom
  andf1$Pos1[ri] <- pos1
  andf1$Pos2[ri] <- pos2
  andf1$WT1[ri] <- wt1
  andf1$SNV1[ri] <-snv1


  }


# Format for ANNOVAR  <15 43762161 43762161 T C>
write.csv(andf1, file = file1 )




