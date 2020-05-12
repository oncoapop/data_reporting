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

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
print("Latest SNP database")
print(dbSNP)

SNP <-   installed.SNPs()
print("Installed SNP database")
print(SNP)

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

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
#for (fname in files)
#{ 

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

# This only assigns the data frame to the primer pairs
#pair <- pr/2
#print(pair)
#if (isTRUE(all(pair == floor(pair)))) assign(paste("Plate", pair, sep=""), d.frame)
}      

# Sections work till here
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

sum1<-sum2
################### READ IN PRIMER ORDER FILE AND REMOVE BITS #####################
# Need to relabel primers with amplicon IDs 
# Read in the primer order file and remove primer= etc and secondary primers

# This file is all in one column

orderdir="/home/dyap/Projects/Tumour_Evol/positions/SNV"
orderfile=paste(paste(orderdir,sample,sep="/"),"p3_order.txt",sep="_")

file <- read.table(orderfile, sep="\n", skip=2)

# Removes the second pair of primer sequences
pat="PRIMER_LEFT_1_SEQUENCE="
file1 <- subset(file, !grepl(pat, file$V1))
pat="PRIMER_RIGHT_1_SEQUENCE="
file2 <- subset(file1, !grepl(pat, file1$V1))

# This removes the primer_0_ headings
clean <- function(ttt){
	gsub("PRIMER_LEFT_0_SEQUENCE=", "", ttt)	
			}
file2[] <- sapply(file2, clean)

clean <- function(ttt){
	gsub("PRIMER_RIGHT_0_SEQUENCE=", "", ttt)
			}
file2[] <- sapply(file2, clean)

####### Inject the chr-pos amplicon ID into sum1 ###########

for (ri in seq(nrow(sum1)))

	{
	saname <-  substr(sum1[ri,3],1,5)
	leftpri <- sum1[ri,4]
	# The right primer is reverse complemented
	rightprir <- sum1[ri,5]
	# Reverse Complement the Right primer 
					x <- DNAString(rightprir)
					rightpri <- as.character(reverseComplement(x))

	#This grep command with $V1 specified returns the COLUMN match
	leftmatch <- grep(leftpri,file2$V1)
	rightmatch <- grep(rightpri,file2$V1)

	if ( length(leftmatch) == 0 | length(rightmatch) == 0 ) 
		{
		if ( length(leftmatch) == 0 ) print("Unmatched Left Primer")
		if ( length(rightmatch) == 0 ) print("Unmatched Right Primer")
		next
		} else {

# The ID chr-pos is 1 before leftmatch and 2 before right match in one column 
file2

        if (file2$V1[leftmatch-1] %in% file2$V1[rightmatch-2])
                      {
                      if (length(file2$V1[rightmatch-2]) != 1 | file2$V1[leftmatch-1] != 1 )
                                 {
				if (length(file2$V1[leftmatch-1]) == 1 ) label <- file2$V1[leftmatch-1]
				if (length(file2$V1[rightmatch-2]) == 1 ) label <- file2$V1[rightmatch-2]
				if (length(file2$V1[rightmatch-2]) != 1 && file2$V1[leftmatch-1] != 1 ) 
						{
						label <- capture.output(cat(file2$V1[leftmatch-1],sep="@"))
						}
                                 } else {
                                         label <- file2$V1[leftmatch-1]
                                        }
                                sum1[ri,3] <- label
                     }
			}
        }

        print("Checking unique amplicons")

######################### Ensure amplicons are unique ####################
# We do this by assigning uniquely assigning the amplicon name to each of the duplicate amplicons

# First delete the duplicated rows in col 3 ie Name
sum2<-unique(sum1)
sum2 <- sum1[!duplicated(sum1[3]),] 
# Alternatively
# sum2 <- sum1[!duplicated(sum1$Name),] 
sum4=NULL

# For this, the ones with multiple names are pulled out (ie with @)
# Then each row is duplicated n+1 times
# Put back the the duplicated row with the multiple IDs removed

for (i in grep("@", sum2$Name))
	{
	dup<-sum2[i,]
	name<-strsplit(sum2[i,3],split="@")
	print(name[[1]])
	shared<-length(name[[1]])

	for (j in seq(shared))
		{ 

		print(name[[1]][j])
		assign(paste("Name",j,sep=""), name[[1]][j])
		sum3<-sum2[rep(i,each=shared),]
		sum3$Name[j]<-get(paste("Name",j,sep=""))
		print(sum3$Name[j])

		# Putting it back into the data.frame
		sum4<-rbind(sum4,sum3)

		}
	sum5<-rbind(sum4,sum2)
	}

sum6<- unique(sum5)
sum7 <- unique(sum6[!duplicated(sum6[3]),])
sum8 <- sum7[! grepl("@", sum7$Name),]

# Checking manually
sum9<-sum8[with(sum8, order(Name)), ]
# edit(sum9)

# For ease , change back to sum1
sum1<-sum9

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
	ampid <- sum1[no,3]
	fwd <- sum1[no,4]
	rev <- sum1[no,5]

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



