#################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels to check primers already made for   		
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jul 2013) 			
##################################################

###############################################################
# Run as source("/home/dyap/R/R-scripts/check_primers_isPCR.R")
# If fails, correct and rerun
# IF for the first time, use firstrow=1
firstrow=1 
################################################################

# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); 
# biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')
# install.packages("XLConnect")
# install.packages("XML", repos = "http://www.omegahat.org/R")
# CRAN (http://cran.r-project.org/)
# install.packages("audio")

library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)
library('BSgenome.Hsapiens.UCSC.hg19')
library(XLConnect)
library(XML)
# For debugging
library(audio)
require(stats)
library(seqinr)

#################################################
# Directory structure and file names

basedir="/share/lustre/backup/dyap/Projects/MiSeq_Data"

primerdir=paste(basedir,"Primer_Order_Files", sep="/")

manifestdir=paste(basedir,"AmpliconManifests", sep="/")

######################################################

# Input #1 is the primer excel sheet
setwd(primerdir)

file_names = list.files(pattern="*Primers*")
file_list = paste(primerdir,file_names, sep="/")
maxfiles<-length(file_list)

# Removing SA494
#files<-file_list[c(1:2,4:maxfiles)]

# Just checking SA494
fname<-file_list[3]

#################

	# Input #2 which is the AmpliconManfest file
	# set input directory
	indir=manifestdir
	# This strips the path and characters after sample name (check when changed)
	sample=substring(fname,66,70)
	manfile=paste(sample,"AmpliconManifest",sep=".")
	infile2=paste(indir, manfile, sep="/")

# workbook <- system.file(file, package = "XLConnect")
# Load workbook
wb <- loadWorkbook(fname, create = FALSE)
# Query available worksheets
sheets <- getSheets(wb)

# This command reads all the sheets in a workbook into a data_list
primers <- readWorksheet(wb, sheets)

sumdir="/home/dyap/Projects/Tumour_Evol/Summary"
PCRfile <- paste(sumdir, paste(sample,"isPCR_chk.csv", sep="_"), sep="/")
bedchkfile <- paste(sumdir, paste(sample,"isPCR_QC.tsv", sep="_"), sep="/")
manfile <- paste(sumdir, paste(sample,"Primers.csv", sep="_"), sep="/")
outfile <- paste(sumdir, paste(sample,"SNP_UCSCman_chk.csv", sep="_"), sep="/")
myfile <- paste(sumdir, paste(sample,"SNP_myman_chk.csv", sep="_"), sep="/")
isinfile <- paste(sumdir, paste(sample,"isPCR_input.tsv", sep="_"), sep="/")
isoutfile <- paste(sumdir, paste(sample,"isPCR_output.tsv", sep="_"), sep="/")

############### Matching the forward and reverse primers ########################
# calculating the number of plates
plates <- length(names(primers))
pairs <- plates/2
if (!isTRUE(all(pairs == floor(pairs)))) stop("Odd number of primer plates cannot be matched")

# If the script continues that means we have even number of plates
# Prepare arrays (data.frames) for primer sets (F and R)

for (pr in seq(plates)) {

	# Fwd and Rev_RC are the primer sequences that were ordered

   d.frame <- data.frame(Plate = rep("", nrow(primers[[plates]])),
		     Well = rep("", nrow(primers[[plates]])),
                     Name = rep("", nrow(primers[[plates]])),
                     Fwd = rep("", nrow(primers[[plates]])),
                     Rev = rep("", nrow(primers[[plates]])),
                     Rev_RC = rep("", nrow(primers[[plates]])),
                     stringsAsFactors = FALSE)

			}      

# This takes the forward primers and puts them into a dataframe.
# This function takes the last n characters of a string

	substrRight <- function(x, n)	{
  					substr(x, nchar(x)-n+1, nchar(x))
					}


for (rw in seq(names(primers)))  {

			platename <- names(primers)[rw]
			plate <- substr(platename,1,nchar(platename)-1)
			plateid <- substr(platename,nchar(platename)-1,nchar(platename))
			test <- substrRight(platename,1)
			if (test == "R")
				next
				else { 
 			for (rx in seq(nrow(primers[[plates]]))) {

					d.frame$Plate[rx] <- plate
					d.frame$Well[rx] <- primers[[rw]]$Well[rx]
					d.frame$Name[rx] <- substr(primers[[rw]]$Name[rx],1,nchar(primers[[rw]]$Name[rx])-2)


						left <-  primers[[rw]]$Sequence[rx]

				# Get the matching REV plate by primer name (cuts F and R suffixies)

						 matchF <- substr(primers[[rw]]$Name[rx],1,nchar(primers[[rw]]$Name[rx])-2) 
						 matchR <- substr(primers[[rw+1]]$Name[rx],1,nchar(primers[[rw+1]]$Name[rx])-2)
						if ( matchF == matchR ) right <- primers[[rw+1]]$Sequence[rx]

						# Removal of adaptors (Fluidigm)
						leftadapt="ACACTGACGACATGGTTCTACA"
						# (5'->3' of reverse adaptor)
						rightadapt="TACGGTAGCAGAGACTTGGTCT"

				# removal of the adaptor sequences from ordered primers
					lpriseq <- gsub(leftadapt, "", left)
					rpriseq <- gsub(rightadapt, "", right)
		
				# Reverse Complement the Right primer without adaptor
					x <- DNAString(rpriseq)
					rpriseqr <- as.character(reverseComplement(x))
		
					d.frame$Fwd[rx] <- lpriseq
					d.frame$Rev[rx] <- rpriseqr
					d.frame$Rev_RC[rx] <- rpriseq
		
				# for first value 
				if ( rw == 1 ) {
						first <- d.frame }
 
 				# Reverse plates are skipped
				# combining successive data.frames
				if ( rw == 3 ) { 
						sum1 <- rbind(first,d.frame) }

				# combining successive data.frames
				if ( rw > 3 ) { 
						sum1 <- rbind(sum1,d.frame) 
							} 

                     					}
					}
				}

print(sum1)

# remove duplictes in primer names
sum2<-sum1[!duplicated(sum1[,3]),]

sum2
################### READ IN PRIMER ORDER FILE AND REMOVE BITS #####################
# Only got the amplicons and position info
# Need to match up with actual order files

orderdir="/share/lustre/backup/dyap/Projects/Tumour_Evol/positions"
orderfile=paste(paste(orderdir,sample,sep="/"),"positions.csv",sep="_")

sum1<-read.csv(file=orderfile, header=TRUE, stringsAsFactors = FALSE)



################## Get hg19 positions from UCSC inSilico PCR web #################################

	manifest=as.data.frame(read.table(infile2, header=TRUE, skip=5, sep="\t", stringsAsFactors = FALSE))

PCR <- data.frame(myID = rep("", nrow(sum1)),
                     UCSCID = rep("", nrow(sum1)),
                     UCSCchr = rep("", nrow(sum1)),
                     UCSCStart = rep(0, nrow(sum1)),
                     UCSCEnd = rep(0, nrow(sum1)),
                     UCSCAmplen = rep(0, nrow(sum1)),
                     UCSCLpri = rep("", nrow(sum1)),
                     UCSCRpri = rep("", nrow(sum1)),
                     UCSCAmp = rep("NA", nrow(sum1)),
                     QCFlag = rep("", nrow(sum1)),                      
                     stringsAsFactors = FALSE)

man <- data.frame(ID = rep("", nrow(sum1)),
                     Chr = rep("", nrow(sum1)),
		     Pos = rep(0, nrow(sum1)),
                     Start = rep(0, nrow(sum1)),
                     End = rep(0, nrow(sum1)),
                     Lprilen = rep(0, nrow(sum1)),
                     Rprilen = rep(0, nrow(sum1)),
                     Description = rep("UCSC-isPCR", nrow(sum1)),
                     stringsAsFactors = FALSE)

isin <- data.frame(ID = rep("", nrow(sum1)),
                     Lpriseq = rep("", nrow(sum1)),
                     Rpriseq = rep("", nrow(sum1)),
                     stringsAsFactors = FALSE)

outdf <- data.frame(ID = rep("", nrow(sum1)),
                     AmpSNP = rep("", nrow(sum1)),
                     LpriSNP = rep("", nrow(sum1)),
                     RpriSNP = rep("", nrow(sum1)),
                     Lpriseq = rep("None", nrow(sum1)),
                     Rpriseq = rep("None", nrow(sum1)),
                     Ampliseq = rep("None", nrow(sum1)),
                     SNPLpriseq = rep("NA", nrow(sum1)),
                     SNPRpriseq = rep("NA", nrow(sum1)),
                     SNPAmpliseq = rep("NA", nrow(sum1)),                      
                     stringsAsFactors = FALSE)

	print("Getting positions from UCSC")

########################### COMMAND LINE SI SILICO PCR #####################################

lastrow<-nrow(sum1)
for (no in (firstrow:lastrow))

	{
	ampid <- paste(sum1$Chr[no],sum1$Pos[no],sep="_")
	fwd <- sum2$Fwd[no]
	rev <- sum2$Rev[no]

	print(ampid)

	PCR$myID[no] <- paste(sample,ampid,sep="_")
	man$ID[no] <- paste(sample,ampid,sep="_")

	
	# This is the most important SNV position!!
	pos <- as.numeric(strsplit(ampid, split="_")[[1]][2])
	man$Pos[no] <- as.numeric(pos)


	isin$ID[no]		<- ampid
	isin$Lpriseq[no]	<- fwd
	isin$Rpriseq[no]	<- rev

	}

write.table(isin, file=isinfile, append = FALSE, quote = FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)

#    isPcr options database query output

isPcr="/share/data/apps/isPcr/bin/x86_64/isPcr"
options="-flipReverse -out=bed"
database="/share/data/apps/isPcr/isPcrSrc/isPcr/data/genomes/twoBit/hg19.2bit"
command=paste(paste(paste(paste(isPcr,options,sep=" "),database, sep=" "), isinfile, sep=" "), isoutfile, sep=" ")

system(command)

bedfile <- read.table(isoutfile, stringsAsFactors = FALSE)

# Label bedfile
names(bedfile)[1] = "chr"
names(bedfile)[2] = "start"
names(bedfile)[3] = "end"
names(bedfile)[4] = "ID"
names(bedfile)[5] = "Score"
names(bedfile)[6] = "Strand"

# Add extra column for QC
bedfile["QCFlag"] <- "Pass"

######################################## QC  #############################################
# While waiting perform the QC
# Enter MiSeq reads
reads <- 151
	
for (ri in seq(nrow(bedfile)))
	{

	chr <- gsub("C", "c", strsplit(bedfile$ID[ri],split="_")[[1]][1])
	pos <- as.numeric(strsplit(bedfile$ID[ri],split="_")[[1]][2])
	
	fiveprime <- pos - as.numeric(bedfile$start[ri])
	threeprime <- as.numeric(bedfile$end[ri]) - pos
	chrom <- bedfile$chr[ri]
	

# Checking overlapping for PE reads

	if ( is.na(fiveprime) | is.na(threeprime) ) {
				print("no UCSC!")
	  			flag4="No UCSC"
	  			bedfile$QCFlag[ri] <- paste(bedfile$QCFlag[ri], flag4, sep=" ") 
				} else {
	if ( fiveprime > reads ) 
		{ print("5' dist insufficient")
		  flag1="5'"
  		  bedfile$QCFlag[ri] <- paste(bedfile$QCFlag[ri], flag1, sep=" ") }
	if ( threeprime > reads ) 
		{ print("3' dist insufficient")
		  flag2="3'"
		  bedfile$QCFlag[ri] <- paste(bedfile$QCFlag[ri], flag2, sep=" ") }
	if ( fiveprime < 0 | threeprime < 0 | chrom != chr ) 
		{ print("FAIL")
		  flag3="SNV offside"	
		  bedfile$QCFlag[ri] <- paste(bedfile$QCFlag[ri], flag3, sep=" ") } 
					}
	
write.table(bedfile, file=bedchkfile, append = FALSE, quote = FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)

	print("Checking SNPs")
	
	#####################################################################################################
	# CHECKING THE PRESENCE OF SNPs (NON SAMPLE SPECIFIC) IN THE AMPLICON AND PRIMERS
	#####################################################################################################


	}
#}


	
print("Done!")


#######################################################################################



