# This script reads all the excel files and extracts ALL the Sample IDs from them
library(XLConnect)
library(plyr)
library(dplyr)
library(RCurl)

# This is where all the relevant excel spreadsheets may lie on MONCO
dir="/Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files/"
setwd(dir)

# This is the RNA-seq file
rnaseqfile="RNA-Seq libraries_Aparicio 2012May16 edited.xls"
rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

# Cut and paste from TN inventory Summary 2015June15.xls (first 4 cols and remove subheadings manually)
inventoryfile="TN inventory Summary 2015June15.xls"
collabID="Sample_matches"
codes<-read.table(collabID, sep="\t", header=TRUE, dec=".",stringsAsFactors = FALSE)

sampleIHC="TN_IHC.XLS"
IHC<-readWorksheetFromFile(file=sampleIHC,sheet = c("Sheet1"), header=TRUE)

# This is from TNBC-27 JIRA ticket cut and paste on 31 May 2017
tnbc="TNBC-27_103.txt"
tnbc_seq<-read.table(tnbc, sep="\t", header=FALSE,stringsAsFactors = FALSE)

# TNBC cases from Peter Jul 2017
xenofile="xeno_TNBC_for_seq.xlsx"
naive_xeno<-readWorksheetFromFile(file=xenofile,sheet = c("TNBC_TxNaive"), header=TRUE)

# Savage Dataset from Adrian Jul 2017
savfile="Savage_SA_IDs.xlsx"
savage<-readWorksheetFromFile(file=savfile,sheet = c("DNA_RNA"), header=TRUE)


# This is where the locations are stored (Arusha showed me Jun 1, 2017)
locfile="Triple Negative_TTR_more sample location_2015June20.xls"
#Sheet1 contains "DNA_T_N & RNA"
location<-readWorksheetFromFile(file=locfile,sheet = c("DNA_T_N & RNA"), header=TRUE)
#Sheet2 contains "RNA"
location.RNA<-readWorksheetFromFile(file=locfile,sheet = c("RNA"), header=TRUE)
#Sheet1 contains "WGA-Triple Negative"
location.WGA<-readWorksheetFromFile(file=locfile,sheet = c("WGA-Triple Negative"), header=TRUE)

LL1<-as.character(location$Sample.Name)
LL2<-c(LL1,as.character(location.RNA$Sample.Name))
LL3<-c(LL2,as.character(naive_xeno$Patient_ID))
LL4<-c(LL3,as.character(savage$PDX.ID))
LL<-unique(c(LL4,as.character(location.WGA$Sample.Name)))
univ<-as.data.frame(LL, stringsAsFactors = FALSE)

# This is how many records there might be
rec=nrow(univ)

# defining the data frame with placeholder data
outdf <- data.frame(CollabID = rep("", rec),
		 addCollabID = rep("", rec),	
		        SAID = rep("", rec),	
		        addID = rep("", rec),	
		     add_SAID = rep("", rec),	
#		     loc_rack = rep("NA", rec),	
#		     loc_shelf = rep("NA", rec),	
#		     loc_box = rep("NA", rec),	
		     loc_pos = rep("NA", rec),	
#		     rem_elut = rep("NA", rec),	
		     sam_type = rep("NA", rec),	
		         PROJ = rep("NA", rec),	
		    RNASEQLIB = rep("NA", rec),	
		    WGSSLIB = rep("NA", rec),	
                     stringsAsFactors = FALSE)

# This is the module where I sum up everything from the sheets
for (ri in seq(rec)) {

		# [,1]<-Collborator ID , c(1) <- col in codes to be output (Collab ID - standardized)
		collID<-gsub("[()]/", "_", as.character(codes[grep(univ[[1]][ri], codes$Sample_ID) , c(1)] ))

		if ( length(collID) > 0 )
				{	
					if ( is.na(unique(collID)) ) { next}
					if ( "collID" == univ[[1]][ri] ) { outdf$CollabID[ri] <- collID } 
					else 	{ 
						outdf$addCollabID[ri] <- collID
						outdf$CollabID[ri] <- codes[codes[,1] %in% collID, c(1)] 
						}
				}


		# [,1] <- use standardized Collborator ID , c(1) <- col in codes to be output (Collab ID - standarized)
		said<-codes[codes[,1] %in% outdf$CollabID[ri], c(1)]

		if ( length(said) > 0 )
				{	
					if ( said == outdf$CollabID[ri] ) { outdf$SAID[ri] <- codes[codes[,1] %in% outdf$CollabID[ri], c(2)]  }
					if ( said == outdf$CollabID[ri] ) { outdf$addID[ri] <- codes[codes[,1] %in% outdf$CollabID[ri], c(3)]  }
					if ( said == outdf$CollabID[ri] ) { outdf$add_SAID[ri] <- codes[codes[,1] %in% outdf$CollabID[ri], c(4)]  }
				}


		# SA ID [,2] and tnbc_seq[,1], c(1) <- col in tnbc_seq to be output (SA ID) 
		dnaid<-tnbc_seq[tnbc_seq[,1] %in% outdf$SAID[ri], c(1)]

		if ( length(dnaid) > 0 )
				{	
					outdf$WGSSLIB[ri] <- tnbc_seq[tnbc_seq[,1] %in% outdf$SAID[ri], c(2)] 
				}

		# SA ID [,2] and rnaseq[,2], c(2) <- col in tnbc_seq to be output (SA ID) 
		rnaid<-rnaseq[rnaseq[,2] %in% outdf$SAID[ri], c(2)]

		if ( length(rnaid) > 0 )
				{	
					outdf$RNASEQLIB[ri] <- rnaseq[rnaseq[,2] %in% outdf$SAID[ri], c("Library")] 
					outdf$PROJ[ri] <- rnaseq[rnaseq[,2] %in% outdf$SAID[ri], c(5)] 			
				}


		# Sample ID column match in location[,1], c(6-9) <- col in location to be output (Rack,Shelf,Box, Pos) 
		sam_loc<-unique(location[location[,1] %in%  outdf$CollabID[ri], c(1)])

		if ( length(sam_loc) > 0 )
				{	
					loc_rack.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(6)],sep=","), collapse=", ")
					loc_shelf.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(7)],sep=","), collapse=", ")
					loc_box.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(8)],sep=","), collapse=", ")
					loc_pos.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(9)],sep=","), collapse=", ")

# This module stores the exect location of the samples 
#					if (length(loc_rack.temp) > 0 ) { outdf$loc_rack[ri] <- loc_rack.temp } 
#					if (length(loc_shelf.temp) > 0 ) { outdf$loc_shelf[ri] <- loc_shelf.temp } 
#					if (length(loc_box.temp) > 0 ) { outdf$loc_box[ri] <- loc_box.temp } 
#					if (length(loc_pos.temp) > 0 ) { outdf$loc_pos[ri] <- loc_pos.temp } 

# This modules says that we have samples only (location not given) We use position to indicate that we have physical sample
# AS long as it has a shelf, box or position reference, we assume that we have samples
					if (length(loc_shelf.temp) > 0 ) { outdf$loc_pos[ri] <- "Sample" } 
					if (length(loc_box.temp) > 0 ) { outdf$loc_pos[ri] <- "Sample" } 
					if (length(loc_pos.temp) > 0 ) { outdf$loc_pos[ri] <- "Sample" } 
				}
			}


# Remove incomplete cases which is no SA ID AND no Collaborator ID 
outdf<-outdf[!(outdf$CollabID=="" & outdf$SAID==""),]

# Missing SA IDs from outdf
tnbc_add<-as.data.frame(unique(sort(tnbc_seq[!tnbc_seq$V1 %in% outdf$SAID,c(1)])))
#rna_add1<-unique(sort(rnaseq[!rnaseq$Sample %in% outdf$SAID, c(2)]))
#rna_add2<-unique(sort(rnaseq[!rnaseq$Sample.1 %in% outdf$SAID, c(3) ]))
rna_add<-as.data.frame(c("SA024","SA029","SA030","SA051","SA052","SA053","SA054","SA055","SA056","SA072"))

colnames(tnbc_add)[1]<-"SAID"
colnames(rna_add)[1]<-"SAID"
missingIDs<-rbind(tnbc_add, rna_add)

rec=nrow(missingIDs)
# defining the data frame with placeholder data
# This second data frame contains all the missing SA IDs which supposedly belong to TNBC 
# but were for some reasons missed from the analysis above 
# this list is partially script generated, partially manually curated.
outdf2 <- data.frame(CollabID = rep("", rec),
                 addCollabID = rep("", rec),
                        SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
                     loc_pos = rep("NA", rec),
                     sam_type = rep("NA", rec),
                         PROJ = rep("NA", rec),
                    RNASEQLIB = rep("NA", rec),
                    WGSSLIB = rep("NA", rec),
                     stringsAsFactors = FALSE)

# This is the module where I get specific missing SA IDs from the respective sheets
for (ri in seq(rec)) {

                # Using the Missing IDs as a guide for this data set
		outdf2$SAID[ri]<-as.character(missingIDs$SAID[ri])

                collID<-gsub("[()]/", "_", as.character(codes[grep(missingIDs$SAID[ri], codes$SA_ID) , c(1)] ))

                if ( length(collID) > 0 )
                                {
                                        outdf2$CollabID[ri] <- codes[codes[,2] %in% outdf2$SAID[ri], c(1)]  
                                        outdf2$addID[ri] <- codes[codes[,2] %in% outdf2$SAID[ri], c(3)]  
                                        #outdf2$add_SAID[ri] <- codes[codes[,1] %in% outdf$CollabID[ri], c(4)]  
				}


                # SA ID [,2] and tnbc_seq[,1], c(1) <- col in tnbc_seq to be output (SA ID)
                dnaid<-tnbc_seq[tnbc_seq[,1] %in% outdf2$SAID[ri], c(1)]

                if ( length(dnaid) > 0 )
                                {
                                        outdf2$WGSSLIB[ri] <- tnbc_seq[tnbc_seq[,1] %in% outdf2$SAID[ri], c(2)]
                                }
                # SA ID [,2] and rnaseq[,2], c(2) <- col in tnbc_seq to be output (SA ID)
                rnaid<-rnaseq[rnaseq[,2] %in% outdf2$SAID[ri], c(2)]

                if ( length(rnaid) > 0 )
                                {
                                        outdf2$RNASEQLIB[ri] <- rnaseq[rnaseq[,2] %in% outdf2$SAID[ri], c("Library")]
                                        outdf2$PROJ[ri] <- rnaseq[rnaseq[,2] %in% outdf2$SAID[ri], c(5)]
                                }


                # Sample ID column match in location[,1], c(6-9) <- col in location to be output (Rack,Shelf,Box, Pos)
                sam_loc<-unique(location[location[,1] %in%  outdf2$CollabID[ri], c(1)])

                if ( length(sam_loc) > 0 )
                                {
                                        loc_rack.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(6)],sep=","), collapse=", ")
                                        loc_shelf.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(7)],sep=","), collapse=", ")
                                        loc_box.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(8)],sep=","), collapse=", ")
                                        loc_pos.temp <- paste(as.character(location[location[,1] %in% codes[ri,1],c(9)],sep=","), collapse=", ")

# This modules says that we have samples only (location not given) We use position to indicate that we have physical sample
# AS long as it has a shelf, box or position reference, we assume that we have samples
                                        if (length(loc_shelf.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                        if (length(loc_box.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                        if (length(loc_pos.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                }
                        }


#summary
output<-join(outdf, outdf2, type = "full")

output<-mutate(output, WGSS=ifelse(output$WGSSLIB=="NA", "None", "Seq'd"))
output<-mutate(output, WTSS=ifelse(output$RNASEQLIB=="NA", "None", "Seq'd"))


print("Transcriptome")
count(output$WTSS)

print("Whole Genome")
count(output$WGSS)

print("TNBC cases")
count(output$PROJ)

print("Have samples")
count(output$loc_pos)

print (rnaseqfile)

print(inventoryfile)

print(sampleIHC)

print("TNBC-27_103.txt")

print(locfile)

# File is in this path on MOMAC14
# /Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files
outfile="TNBC_Summary_by_all_SA_onMAC_3.csv"
write.table(output, file=outfile, sep=",")

# Format
# cat (outputfile) | perl -p -e 's/\r\n|\n|\r/\|\n/g' | tr "," "|" | sed 's/^/\|/' | sed 's/\"\"/\ /g' | tr -d '"'

