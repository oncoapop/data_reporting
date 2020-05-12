# This script reads all the excel files and extracts ALL the Sample IDs from them
library(XLConnect)

# This is where all the relevant excel spreadsheets may lie on MONCO
dir="/Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files/"
setwd(dir)

# This is the RNA-seq file
rnaseqfile="RNA-Seq libraries_Aparicio 2012May16.xls"
rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

# Cut and paste from TN inventory Summary 2015June15.xls (first 4 cols and remove subheadings manually)
collabID="Sample_matches"
codes<-read.table(collabID, sep="\t", header=TRUE, dec=".",stringsAsFactors = FALSE)

sampleIHC="TN_IHC.XLS"
IHC<-readWorksheetFromFile(file=samplefile,sheet = c("Sheet1"), header=TRUE)

# This is from TNBC-27 JIRA ticket cut and paste on 31 May 2017
tnbc="TNBC-27_103.txt"
tnbc_seq<-read.table(tnbc, sep="\t", header=FALSE,stringsAsFactors = FALSE)

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
LL<-unique(c(LL2,as.character(location.WGA$Sample.Name)))
univ<-as.data.frame(LL, stringsAsFactors = FALSE)

# This is how many records there might be
rec=nrow(univ)

# Need to get the Sa IDs from all the files
# location[grep("TB02.0318", location$Sample.Name) ,] # an exmaple of grep


# defining the data frame with placeholder data
outdf <- data.frame(CollabID = rep("NA", rec),
		 addCollabID = rep("NA", rec),	
		        SAID = rep("NA", rec),	
		        addID = rep("NA", rec),	
		     add_SAID = rep("NA", rec),	
		     loc_rack = rep("NA", rec),	
		     loc_shelf = rep("NA", rec),	
		     loc_box = rep("NA", rec),	
		     loc_pos = rep("NA", rec),	
		     rem_elut = rep("NA", rec),	
		     sam_type = rep("NA", rec),	
		     sam_type = rep("NA", rec),	
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
				} else {
					outdf$CollabID[ri] <- univ[[1]][ri]
					}
			}


write.table(outdf, file="TNBC_27_matches", sep=",")
