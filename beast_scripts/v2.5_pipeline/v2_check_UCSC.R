#################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels to check primers already made for     	
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Aug 2013) 			
##################################################

###############################################################
# Run as source("/meta/o/oncoapop/Scripts/v2_check_UCSC.R")
# IF for the first time, use firstrow=1

run=1
firstrow=1

# Remember to change the output file to run number
# Change the input / output path

################################################################

# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); 
# source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); 
# biocLite("BSgenome.Hsapiens.UCSC.hg19"); 
# library('BSgenome.Hsapiens.UCSC.hg19')
# install.packages("XLConnect")
# install.packages("XML", repos = "http://www.omegahat.org/R")
# CRAN (http://cran.r-project.org/)
# install.packages("audio")

library(XML)
# For debugging
#library(audio)
require(stats)


################################################################

# if run directly uncomment the sample name
# Command line `Rscript v2_GetSeq.R --no-save --no-restore 
# --args $dir%$sample%$file%$Miseq`

# This takes the 4th argument (see str above) which is sample name
args <- commandArgs(trailingOnly = TRUE)
input <- args[4]

# To test this programme in R using source
# commandArgs <- function() "TEST%123%20130926214630%150"
# source(file="/meta/o/oncoapop/Scripts/v2_pipeline/v2_check_UCSC.R")
# For testing only uncomment for production
# input <- "TEST%Case2%20130927231049%150"

Project <- strsplit(input, split="%")[[1]][1]
sample <- strsplit(input, split="%")[[1]][2]
posfile <- strsplit(input, split="%")[[1]][3]
Miseq <- as.numeric(strsplit(input, split="%")[[1]][4])

print("Directory")
print(Project)
print("Sample_ID")
print(sample)
print("File")
print(posfile)
print("MiSeq Read")
print(Miseq)

sumdir <- paste(paste(paste("/meta/o/oncoapop/Projects/Pipeline", Project, sep="/"), sample, sep="/"), "primer3", sep="/")

################################################################
print("################################################################")

timestamp <- Sys.time()
print(timestamp)

PCRfile <- paste(sumdir, paste(posfile,"isPCR_chk.csv", sep="_"), sep="/")
manfile <- paste(sumdir, paste(posfile,"Primers.csv", sep="_"), sep="/")
#outfile <- paste(sumdir, paste(posfile,"SNP_UCSCman_chk.csv", sep="_"), sep="/")
#myfile <- paste(sumdir, paste(posfile,"SNP_myman_chk.csv", sep="_"), sep="/")
#pipeline <- paste(sumdir, paste(posfile,"isPCR-output.txt", sep="_"), sep="/")

p3dir=sumdir
tmpdir="/meta/o/oncoapop/temp"
infile=paste(paste(sumdir, posfile, sep="/"), "primerlist.txt",sep="_")

# Format of file (from primer3pipeline.sh)
# Sample_chr2_95542476,chr2_95542476,TCTCTTCATCGACCGCCAGAAG,CCAGTTCAGATGGACAGCAAC,189

print("Read input file")
anyfile<-as.data.frame(read.csv(infile, header=FALSE))

# Removes incomplete enteries
primerfile <- anyfile[complete.cases(anyfile),]

primers <- data.frame (ID = rep("", nrow(primerfile)),
		       Chr = rep("", nrow(primerfile)),
		       Pos = rep(0, nrow(primerfile)),
		       SA  = rep("", nrow(primerfile)),
		       Left = rep("", nrow(primerfile)),
		       Right = rep("", nrow(primerfile)),
		       Length = rep(0, nrow(primerfile)),
		 	stringsAsFactors = FALSE)
print("sort primers")
timestamp <- Sys.time()
print(timestamp)

for ( ri in seq(nrow(primerfile)))
	{
	primers$ID[ri] <- as.character(primerfile$V1[ri])
	primers$Chr[ri] <- strsplit(as.character(primerfile$V2[ri]),split="_")[[1]][1]
	primers$Pos[ri] <- as.numeric(strsplit(as.character(primerfile$V2[ri]),split="_")[[1]][2])
	primers$SA[ri] <- strsplit(as.character(primerfile$V1[ri]),split="_")[[1]][1]
	primers$Left[ri] <- as.character(primerfile$V3[ri])
	primers$Right[ri] <- as.character(primerfile$V4[ri])
	if (length(primerfile$V5[ri] >0 ))	primers$Length[ri] <- as.numeric(primerfile$V5[ri])
	} 	

print(primers)
print("done")
timestamp <- Sys.time()
print(timestamp)

################## Get hg19 positions from UCSC inSilico PCR web ###############################


PCR <- data.frame(myID = rep("", nrow(primers)),
                     UCSCID = rep("", nrow(primers)),
                     UCSCchr = rep("", nrow(primers)),
                     UCSCStart = rep(0, nrow(primers)),
                     UCSCEnd = rep(0, nrow(primers)),
                     UCSCAmplen = rep(0, nrow(primers)),
                     UCSCLpri = rep("", nrow(primers)),
                     UCSCRpri = rep("", nrow(primers)),
                     UCSCAmp = rep("NA", nrow(primers)),
                     QCFlag = rep("", nrow(primers)),
                     stringsAsFactors = FALSE)

man <- data.frame(ID = rep("", nrow(primers)),
                     Chr = rep("", nrow(primers)),
		     Pos = rep(0, nrow(primers)),
                     Start = rep(0, nrow(primers)),
                     End = rep(0, nrow(primers)),
                     Lprilen = rep(0, nrow(primers)),
                     Rprilen = rep(0, nrow(primers)),
                     Description = rep("UCSC-isPCR", nrow(primers)),
                     stringsAsFactors = FALSE)


	print("Getting positions from UCSC")

########################### USCS CONNECTION FOR SI SILICO PCR ################################
lastrow=nrow(primers)

for (no in (firstrow:lastrow))

	{
	print("Checking...")
	timestamp <- Sys.time()
	print(timestamp)

	ampid <-paste(primers$Chr[no],primers$Pos[no],sep="_")

	fwd <- primers$Left[no]
	rev <- primers$Right[no]
	
	add1="http://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=345005557&org=Human&db=hg19&wp_target=genome&wp_f="
	add2="&Submit=submit&wp_size=800&wp_perfect=17&wp_good=17&boolshad.wp_flipReverse=0"
	sep="&wp_r="
	
	url=paste(paste(paste(add1,fwd,sep=""),rev,sep=sep),add2,sep="")

	doc2<-readHTMLTable(url)
	parse<-as.character(doc2[[1]][2,1])
	print(ampid)

	PCR$myID[no] <- ampid
	man$ID[no] <- paste(sample,ampid,sep="_")

# They usually change this just to throw the automated scripts out!
#	print(strsplit(parse,split=">")[[1]][4])
	string <- strsplit(parse,split=">")[[1]][4]

	PCR$UCSCchr[no] <- strsplit(string,split=":")[[1]][1]
	man$Chr[no] <- strsplit(string,split=":")[[1]][1]

	# This is the most important SNV position!!
	pos <- primers$Pos[no]
	man$Pos[no] <- as.numeric(pos)

	# This is to split the start and end coords in the format "nnnnn+nnnnn"
	coord <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][1]
	coord1 <- as.character(gsub('([[:punct:]])','\n', coord))
	PCR$UCSCStart[no] <- strsplit(coord1, split="\n")[[1]][1]
	PCR$UCSCEnd[no] <- strsplit(coord1, split="\n")[[1]][2]
	man$Start[no] <- as.numeric(strsplit(coord1, split="\n")[[1]][1])
	man$End[no] <- as.numeric(strsplit(coord1, split="\n")[[1]][2])

	# amplen nnbp (substr removes bp)
	amplen <- strsplit(strsplit(string,split=":")[[1]][2],split=" ")[[1]][2]
	print(amplen)

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


print(PCR)
	######################################## QC ###########################################
	# While waiting perform the QC
	# Enter MiSeq reads
	reads <- Miseq

	system("sleep 15s")

	print(paste("QC for Miseq run of bp =", reads, sep=" "))

	# Match p3 amplicon lenth with UCSC Ampl length

	UCSClen <-  as.numeric(PCR$UCSCAmplen[no])
        p3len <- primers$Length[no]

	if (is.na(UCSClen) | is.na(p3len)) 
		{
		  print("One Length missing")
                  flag5="Unknown"
                  PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag5, sep=" ") 

		} else {
		  
 
        		if ( UCSClen == p3len ) 
				{ 
		  print("Amplicon length MATCH!")
                  flag5="Validated"
                  PCR$QCFlag[no] <- paste(PCR$QCFlag[no], flag5, sep=" ") 
				}

			}

	# Check MiSeq amplification characteristics

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

	print("Writing to file...")
		
	write.csv(PCR, file=PCRfile)
	write.csv(man, file=manfile)

	}

print("DONE!")
timestamp <- Sys.time()
print(timestamp)
endtime <- system("TZ=GMT+7 date | sed 's/GMT/PDT/'")
print(endtime)
print("################################################################")


