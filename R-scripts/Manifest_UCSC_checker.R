#################################################
##  Script to check UCSC positions of the     
##          primers already made for     	
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Aug 2013) 			
##################################################

###############################################################
# Run as source("/home/dyap/R/R-scripts/Manifest_UCSC_checker.R")
#
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

#################################################
# Directory structure and file names

# ASROCK II
# basedir="E:\\DATA\\DOCUMENTS\\WORK\\QC"
# Beast.Cluster
basedir="/home/dyap/Projects/MiSeq_Data"

primerdir=paste(basedir,"Primer_Order_Files", sep="/")
#primerdir=basedir

manifestdir=paste(basedir,"AmpliconManifests", sep="/")
#manifestdir=basedir

######################################################

# Input #1 is the primer excel sheet
setwd(primerdir)

file_names = list.files(pattern="*Primers*.*xlsx")
files = paste(primerdir,file_names, sep="/")
maxfiles<-length(files)

#################
for (fname in files[1:maxfiles])
{

#For testing
# print (fname)}
#fname<-files[2]

	# Input #2 which is the AmpliconManfest file
	# set input directory
	indir=manifestdir
	# This strips the path and characters after sample name (check when changed)
	sample=substring(strsplit(fname,split="/")[[1]][2],1,5)
	manfile=paste(sample,"AmpliconManifest",sep=".")
	infile2=paste(indir, manfile, sep="/")

# workbook <- system.file(file, package = "XLConnect")
# Load workbook
wb <- loadWorkbook(fname, create = FALSE)
# Query available worksheets
sheets <- getSheets(wb)

# This command reads all the sheets in a workbook into a data_list
primers <- readWorksheet(wb, sheets)

sumdir=basedir

outfile <- paste(sumdir, paste(sample,"SupplFig", sep="_"), sep="/")

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

##################################################################################



################# LOAD THE isPCR PRIMER FILE #####################################

#For testing
# print (fname)}
#fname<-files[2]

	# Input #2 which is the AmpliconManfest file
	# set input directory
	indir=manifestdir
	# This strips the path and characters after sample name (check when changed)
	sample=substring(strsplit(fname,split="/")[[1]][2],1,5)
	manfile=paste(sample,"AmpliconManifest",sep=".")
	infile2=paste(indir, manfile, sep="/")

	prim <- read.csv(fname, header=TRUE, stringsAsFactors = FALSE)

	###############LOAD THE RESPECTIVE AmpliconManifest ###################

	# Check that the files are in unix format
	# unix commmand line: awk '{ gsub("\r", "\n"); print $0;}' SA493.AmpliconManifest.txt > SA493.AmpliconManifest

	manifest=as.data.frame(read.table(infile2, header=TRUE, skip=5, sep="\t", stringsAsFactors = FALSE))

	############# GENERATED NEW MANIFEST WITH UCSC CHECKED COORDINATES #################################

	newman <- data.frame(Name = rep("", nrow(manifest)),
		     Chromosome = rep("", nrow(manifest)),
                     Amplicon.Start = rep(0, nrow(manifest)),
                     Amplicon.End = rep(0, nrow(manifest)),
                     Upstream.Probe.Length = rep(0, nrow(manifest)),
		     Downstream.Probe.Length = rep(0, nrow(manifest)),
                     Comments = rep("", nrow(manifest)),
                     stringsAsFactors = FALSE)

	for (man in seq(nrow(manifest)))
		{
		# generate the match from the manifest file
		match=paste(sample,manifest$Name[man],sep="_")

		if (match %in% prim$ID)
			{
			newman$Name[man] <- match
			newman$Chromosome[man] <- manifest$Chromosome[man]
			newman$Comments[man] <- manifest$Comments[man]

			# match the manifest ID with the ID in primer file (col2 ie [,2]) generated by UCSC
			# for left primer
			manleftlen <- as.numeric(manifest$Upstream.Probe.Length[man])
			ucscleftlen <- as.numeric(prim[prim[,2] %in% c(match),7])
			ucscstart <- as.numeric(prim[prim[,2] %in% c(match),5])
			manstart <- as.numeric(manifest$Amplicon.Start[man])
				if (ucscstart == manstart+2) 
					{
					newman$Amplicon.Start[man] <- ucscstart
					} 
				if (ucscleftlen == manleftlen-1) 
					{
					newman$Upstream.Probe.Length[man] <- ucscleftlen
					} 

				if (length(ucscstart) == 0)   newman$Amplicon.Start[man] <- manstart+2
				if (length(ucscleftlen) == 0) newman$Upstream.Probe.Length[man] <- manleftlen-1

			# for rightprimer
			manrightlen <- as.numeric(manifest$Downstream.Probe.Length[man])
			ucscrightlen <- as.numeric(prim[prim[,2] %in% c(match),8])
			ucscend <- as.numeric(prim[prim[,2] %in% c(match),6])
			manend <- as.numeric(manifest$Amplicon.End[man])
				if (ucscend == manend-2) 
					{
					newman$Amplicon.End[man] <- ucscend
					} 
				if (ucscrightlen == manrightlen-1)
					{
					newman$Downstream.Probe.Length[man] <- ucscrightlen
					} 
				if (length(ucscend) == 0)      newman$Amplicon.End[man] <- endstart-2
				if (length(ucscrightlen) == 0) newman$Downstream.Probe.Length[man] <- manrightlen-1

				print(ucscrightlen)
				print(manrightlen)
				print(ucscend)
				print(manend)
				
			} 

 		
		}

}











####################################################################################









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

# The ID chr-pos is 1 before leftmatch and 2 before right match in one column file2

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
edit(sum9)

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

########################### USCS CONNECTION FOR SI SILICO PCR #####################################

lastrow<-nrow(sum1)
for (no in (firstrow:lastrow))

	{
	ampid <- sum1[no,3]
	fwd <- sum1[no,4]
	rev <- sum1[no,5]
	add1="http://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=342829929&org=Human&db=hg19&wp_target=genome&wp_f="
	add2="&Submit=submit&wp_size=400&wp_perfect=17&wp_good=17&wp_flipReverse=on&boolshad.wp_flipReverse=0"
	sep="&wp_r="
	url=paste(paste(paste(add1,fwd,sep=""),rev,sep=sep),add2,sep="")

	doc2<-readHTMLTable(url)
	parse<-as.character(doc2[[1]][2,1])
	print(ampid)

	PCR$myID[no] <- paste(sample,ampid,sep="_")
	man$ID[no] <- paste(sample,ampid,sep="_")

	
	string <- strsplit(parse,split=">")[[1]][2]
	PCR$UCSCchr[no] <- strsplit(string,split=":")[[1]][1]
	man$Chr[no] <- strsplit(string,split=":")[[1]][1]

	# This is the most important SNV position!!
	pos <- as.numeric(strsplit(ampid, split="_")[[1]][2])
	man$Pos[no] <- as.numeric(pos)

	# This is to split the start and end coords in the format "nnnnn+nnnnn"
	coord <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][1]
	coord1 <- as.character(gsub('([[:punct:]])','\n', coord))
	PCR$UCSCStart[no] <- strsplit(coord1, split="\n")[[1]][1]
	PCR$UCSCEnd[no] <- strsplit(coord1, split="\n")[[1]][2]
	man$Start[no] <- as.numeric(strsplit(coord1, split="\n")[[1]][1])
	man$End[no] <- as.numeric(strsplit(coord1, split="\n")[[1]][2])

	amplen <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][2]
	PCR$UCSCAmplen[no] <- substr(amplen,1,nchar(amplen)-2)
	PCR$UCSCLpri[no] <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][3]
	man$Lprilen[no] <- nchar(PCR$UCSCLpri[no])

	str2 <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][4]

	PCR$UCSCRpri[no] <- strsplit(str2,split="\n")[[1]][1]
	man$Rprilen[no] <- nchar(PCR$UCSCRpri[no])

	# This is the string length of the amplicon split into bits by \n
	strlen <- length(strsplit(str2,split="\n")[[1]])

	PCR$UCSCAmp[no] <- paste(strsplit(str2,split="\n")[[1]][2:strlen],collapse="")
	PCR$UCSCID[no] <- paste(paste(sample,PCR$UCSCchr[no],sep="_"), as.character(gsub('([[:punct:]])','-', coord)), sep=":")

	##########################################################################################	

	# If submitted automatedly it is recommended to wait 15 sec between requests
	wait(15)

	print("Wait 15s")



	######################################## QC  #############################################
	# While waiting perform the QC
	# Enter MiSeq reads
	reads <- 151
	
	fiveprime <- pos - as.numeric(PCR$UCSCStart[no])
	threeprime <- as.numeric(PCR$UCSCEnd[no]) - pos
	
	# Checking overlapping for PE reads

	if ( is.na(fiveprime) | is.na(threeprime) ) {
					print("no UCSC!")
		  			flag4="No UCSC"
		  			PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag4, sep=" ") 
				} else {
	if ( fiveprime > reads ) 
		{ print("5' dist insufficient")
		  flag1="5'"
		  PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag1, sep=" ") }
	if ( threeprime > reads ) 
		{ print("3' dist insufficient")
		  flag2="3'"
		  PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag2, sep=" ") }
	if ( fiveprime < 0 | threeprime < 0) 
		{ print("FAIL")
		  flag3="SNV offside"	
		  PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag3, sep=" ") }
					}

	print("Checking SNPs")
	
	#####################################################################################################
	# CHECKING THE PRESENCE OF SNPs (NON SAMPLE SPECIFIC) IN THE AMPLICON AND PRIMERS
	#####################################################################################################

	chr <- man$Chr[no]
	start <- as.numeric(man$Start[no])
	end <- as.numeric(man$End[no])
	leftlen <- as.numeric(man$Lprilen[no])
	rgtlen <- as.numeric(man$Rprilen[no])

	if ( is.na(chr) | is.na(start) | is.na(end)) {
						outdf$ID[no] <- paste(sample,ampid,sep="_")
						outdf$Ampliseq[no] <- "NO AMPLICON"
						next
					} else {

	# Get the amplicons 1. Reference and 2. SNPmasked
 						ampliseq <- as.character(getSeq(Hsapiens,chr,start,end))
 						snpampliseq <- as.character(getSeq(SNP_Hsapiens,chr,start,end))
 
	# Get left and right primers 1. Designed Seq 2. SNPmasked 
						leftend <- start + leftlen - 1
						lpriseq <- PCR$UCSCLpri[no]
						snplpriseq <- as.character(getSeq(SNP_Hsapiens,chr,start,leftend))

						rightstart <- end - rgtlen + 1
						snprpriseq <- as.character(getSeq(SNP_Hsapiens,chr,rightstart,end))
	# Reverse Complement the Right primer 
					x <- DNAString(PCR$UCSCRpri[no])
					rpriseq <- as.character(reverseComplement(x))

		# Testing to see if the sequence are identical
		# if they are not identical they contain SNPs or are wrong

				if (ampliseq == snpampliseq) ampsnp <- "ok" else { ampsnp <- "SNP" }
				if (lpriseq == snplpriseq) lprisnp <- "ok" else { lprisnp <- "SNP" }
				if (rpriseq == snprpriseq) rprisnp <- "ok" else { rprisnp <- "SNP" }


	# writing the output to the dataframe 
			outdf$ID[no] <- paste(sample,ampid,sep="_")
                      	outdf$AmpSNP[no] <- ampsnp
                     	outdf$LpriSNP[no] <- lprisnp
                     	outdf$RpriSNP[no] <- rprisnp
                     	outdf$Lpriseq[no] <- lpriseq
                     	outdf$Rpriseq[no] <- rpriseq
                     	outdf$Ampliseq[no] <- ampliseq
                     	outdf$SNPLpriseq[no] <- snplpriseq
                     	outdf$SNPRpriseq[no] <- snprpriseq
                     	outdf$SNPAmpliseq[no] <- snpampliseq

						}



	########################################################################################################

write.csv(PCR, file=PCRfile)
write.csv(man, file=manfile)
write.csv(outdf, file=outfile)

}

print(sample)
print("Done!")

}

#######################################################################################



