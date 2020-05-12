# This script reads all the excel files and extracts ALL the Sample IDs from them
library(XLConnect)
library(plyr)
library(dplyr)
library(RCurl)

# This is where all the relevant excel spreadsheets may lie on MONCO
dir="/Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files/"
setwd(dir)

# sample pairings
collabID="Sample_matches"
match<-read.table(collabID, sep="\t", header=TRUE, dec=".",stringsAsFactors = FALSE)
	
# This is the RNA-seq file
rnaseqfile="RNA-Seq libraries_Aparicio 2012May16 edited.xls"
rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

# Cut and paste from TN inventory Summary 2015June15.xls (first 4 cols and remove subheadings manually)
inventoryfile="TN inventory Summary 2015June15.xls"
master_IHC<-readWorksheetFromFile(file=inventoryfile,sheet = c("Master_seq_IHC"), header=TRUE)
master_summary<-readWorksheetFromFile(file=inventoryfile,sheet = c("summary"), header=TRUE)

sampleIHC="TN_IHC.XLS"
IHC<-readWorksheetFromFile(file=sampleIHC,sheet = c("Sheet1"), header=TRUE)

# This is from TNBC-27 JIRA ticket cut and paste on 31 May 2017
#tnbc="TNBC-27_103.txt"
tnbc="TNBC_103_LIB.txt"
tnbc_seq<-read.table(tnbc, sep="\t", header=FALSE,stringsAsFactors = FALSE)
colnames(tnbc_seq)[1]<-"SA_ID"
colnames(tnbc_seq)[2]<-"WGSS_LIB_ID"

# This is where the locations are stored (Arusha showed me Jun 1, 2017)
locfile="Triple Negative_TTR_more sample location_2015June20.xls"
#Sheet1 contains "DNA_T_N & RNA"
location<-readWorksheetFromFile(file=locfile,sheet = c("DNA_T_N & RNA"), header=TRUE)
#Sheet2 contains "RNA"
location.RNA<-readWorksheetFromFile(file=locfile,sheet = c("RNA"), header=FALSE)
#Sheet1 contains "WGA-Triple Negative"
location.WGA<-readWorksheetFromFile(file=locfile,sheet = c("WGA-Triple Negative"), header=TRUE)

# SA-ID universe
MM1<-rnaseq[!is.na(rnaseq$Sample), c("Sample")]
MM2<-master_summary[!is.na(master_summary$SA.ID),c("SA.ID")]
MM3<-tnbc_seq[!is.na(tnbc_seq$SA_ID),c("SA_ID")]
saID<-unique(do.call(c, list(MM1,MM2,MM3)))
newsaid<-saID[grep("N",saID, invert = TRUE)]

test1<-as.data.frame(newsaid[!newsaid %in% match$SA_ID])
rec=length(test1[grep("SA",test1[,1]),])

test <-  data.frame(Sample_ID = seq.int(rec),
                        SA_ID = rep("", rec),
                Additional_ID = rep("", rec),
            Additional_SA_IDs = rep("", rec),
             stringsAsFactors = FALSE)

test$SA_ID<-test1[grep("SA",test1[,1]),]
tmp<-match
codes<-rbind(tmp,test)

# CHECK - make sure no SA IDs appear in this
codes$SA_ID[duplicated(codes$SA_ID)]

# Patient ID universe
TN<-tnbc_seq$SA_ID
var1<-unique(master_summary$Sample.ID[1:211])
# removing all the heading with the word "cases"
sum1<-var1[lapply(var1,function(x) length(grep("cases",x,value=FALSE))) == 0]
LL1<-as.character(unique(gsub(" ","",location$Sample.Name)))
LL2<-c(LL1,as.character(unique(gsub(" ","",location.RNA$Sample.Name))))
LL3<-c(LL2,as.character(unique(gsub(" ","",location.WGA$Sample.Name))))
LL4<-c(LL3,codes$Sample_ID)
#LL5<-c(LL4,tnbc_seq[grep("T",tnbc_seq$SA_ID),c("SA_ID")]) # These are the tumours incorrectly labelled with T and some xenos

someID<-unique(do.call(c, list(LL1,LL2,LL3,LL4,sum1)))
univ<-as.data.frame(someID, stringsAsFactors = FALSE)



# Files which contain cellularity and IHC information
calfile="sample_cellularity_ploidy.xls"
ABfile="AlbertaTNSummary.xls"
TBfile="TN_MB_clinicalFile.xls"
TNBCfile="TNBCa_in seq queue_IHC_May012014.xls"
VBAfile="Frozen H&E review for BrCa samples with matched normals.xlsx"

#Sheet1 contains the calculated cellularity
cell.cal<-readWorksheetFromFile(file=calfile,sheet = c("sample_cellularity_ploidy.txt"), header=TRUE)
cell.AB<-readWorksheetFromFile(file=ABfile,sheet = c("GT file"), header=TRUE, startRow = 2)
cell.TB<-readWorksheetFromFile(file=TBfile,sheet = c("Clinical Data"), header=TRUE, startRow = 23)
cell.TNBC<-readWorksheetFromFile(file=TNBCfile,sheet = c("Sheet1"), header=TRUE, startRow = 1)
cell.VBA<-readWorksheetFromFile(file=VBAfile,sheet = c("Histo_GT"), header=TRUE, startRow = 1)


# This is how many records there might be
rec=nrow(univ)

# defining the data frame with placeholder data
outdf <-  data.frame(CollabID = rep("", rec),
                  addCollabID = rep("", rec),
                         SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
#                    loc_rack = rep("NA", rec),
#                   loc_shelf = rep("NA", rec),
#                     loc_box = rep("NA", rec),
                      loc_pos = rep("NA", rec),
#                    rem_elut = rep("NA", rec),
#                    sam_type = rep("NA", rec),
                Tumour.Normal = rep("NA", rec),
                       ER_TMA = rep("NA", rec),
                       PR_TMA = rep("NA", rec),
                     HER2_TMA = rep("NA", rec),
                     EGFR_TMA = rep("NA", rec),
                    CK5.6_TMA = rep("NA", rec),
                    RNASeq_Mb = rep("NA", rec),
                    GSC_EXCAP = rep("NA", rec),
                    SOLID_WGS = rep("NA", rec),
                 Add_SOLIDLIB = rep("NA", rec),
                 SNP_ARRAY_ID = rep("NA", rec),
                  METABRIC_ID = rep("NA", rec),
                         PROJ = rep("NA", rec),
                       SOURCE = rep("NA", rec),
                    RNASEQLIB = rep("NA", rec),
                      WGSSLIB = rep("NA", rec),
                 NORM_WGSSLIB = rep("NA", rec),

                       Normal = rep("NA", rec),
                       Tumour = rep("NA", rec),
                RNASEQLIB.seq = rep("NA", rec),
                  WGSSLIB.seq = rep("NA", rec),
             NORM_WGSSLIB.seq = rep("NA", rec),

                     Cell.cal = rep("", rec),
                     Cell.obs = rep("", rec),
             stringsAsFactors = FALSE)

# This is the module where I sum up everything from the sheets
for (ri in seq(rec)) {

                # Using univ$someID to grep names from codes$Sample_ID (which is the Collaborator ID wanted) from codes
                collID<-gsub("[()]/", "_", as.character(codes[grep(univ$someID[ri], codes$Sample_ID) , c("Sample_ID")] ))

                if ( sum(nchar(collID)) > 0 )
                                {
                                        if ( is.na(unique(collID)) ) 	{ next } 
                                        if ( "collID" == univ$someID[ri] ) { outdf$CollabID[ri] <- collID }
                                        else    {
                                                outdf$addCollabID[ri] <- collID
                                                outdf$CollabID[ri] <- codes[codes$Sample_ID %in% collID, c("Sample_ID")]
                                                outdf$SAID[ri] <- codes[codes$Sample_ID %in% collID, c("SA_ID")]
                                                }
                                }


                # Used the standardized Collborator ID (Sample_ID) , c(1) <- col in codes to be output (Sample_ID)
                said<-codes[codes$Sample_ID %in% outdf$CollabID[ri], c("Sample_ID")]

                if ( sum(nchar(said)) > 0 )
                                {
                                        if ( said == outdf$CollabID[ri] ) { outdf$SAID[ri] <- codes[codes$Sample_ID %in% outdf$CollabID[ri], c("SA_ID")]  }
                                        if ( said == outdf$CollabID[ri] ) { outdf$addID[ri] <- codes[codes$Sample_ID %in% outdf$CollabID[ri], c("Additional_ID")]  }
                                        if ( said == outdf$CollabID[ri] ) { outdf$add_SAID[ri] <- codes[codes$Sample_ID %in% outdf$CollabID[ri], c("Additional_SA_IDs")]  }
                                }

                # For UPDATE: use WGSS library from TNBC-27 JIRA (103 genomes) and confirm that they are TNBC
                dnaid<-tnbc_seq[tnbc_seq$SA_ID %in% outdf$SAID[ri], c("SA_ID")]

                if ( sum(nchar(dnaid)) > 0 )
                                {
                                        outdf$WGSSLIB[ri] <- tnbc_seq[tnbc_seq$SA_ID %in% outdf$SAID[ri], c("WGSS_LIB_ID")]
                                        outdf$NORM_WGSSLIB[ri] <- tnbc_seq[tnbc_seq$SA_ID %in% paste(outdf$SAID[ri],"N",sep=""), c("WGSS_LIB_ID")]
                                        outdf$PROJ[ri] <- "TNBC_103"
                                }

		# SA ID [,2] and rnaseq$Sample, c("Sample") <- col in tnbc_seq to be output (SA ID) 
		rnaid<-rnaseq[rnaseq$Sample %in% outdf$SAID[ri], c("Sample")]

		if ( sum(nchar(rnaid)) == 0 )	
				{ 
					rnaid<-rnaseq[rnaseq$Sample.1 %in% outdf$SAID[ri], c("Sample.1")] 
					if ( sum(nchar(rnaid)) > 0 )	
						{ 
							outdf$RNASEQLIB[ri] <- rnaseq[rnaseq$Sample.1 %in% outdf$SAID[ri], c("Library")] 
							outdf$PROJ[ri] <- rnaseq[rnaseq$Sample.1 %in% outdf$SAID[ri], c("Description")] 			
							rnaid<-NULL 			
						}
				}
		if ( sum(nchar(rnaid)) > 0 )	
				{ 
					outdf$RNASEQLIB[ri] <- rnaseq[rnaseq$Sample %in% outdf$SAID[ri], c("Library")] 
					outdf$PROJ[ri] <- rnaseq[rnaseq$Sample %in% outdf$SAID[ri], c("Description")] 			
				}
	
		# Collab ID [,1] and location$Sample.Name, c("Sample.Name") <- col in location to be output (Collab ID) 
		locid<-unique(location[location$Sample.Name %in% outdf$CollabID[ri], c("Sample.Name")])

		if ( sum(nchar(locid)) > 0 )
				{	
					source <- unique(location[location$Sample.Name %in% outdf$CollabID[ri], c("Source")])
					if ( sum(nchar(source)) > 1) 
						{
						outdf$SOURCE[ri] <- paste(unique(location[location$Sample.Name %in% outdf$CollabID[ri], c("Source")]), collapse=",")
						}
				}

## SPLIT HERE
		# Sample ID column match in location[,1], c(6-9) <- col in location to be output (Rack,Shelf,Box, Pos) 
		sam_loc<-unique(location[location$Sample.Name %in%  outdf$CollabID[ri], c(1)])

		if ( sum(nchar(sam_loc)) > 0 )
				{	
					loc_rack.temp <- paste(as.character(location[location$Sample.Name %in% outdf$CollabID[ri], c("Rack..")],sep=","), collapse=", ")
					loc_shelf.temp <- paste(as.character(location[location[,1] %in% outdf$CollabID[ri], c("Shelf..")], sep=","), collapse=", ")
					loc_box.temp <- paste(as.character(location[location[,1] %in% outdf$CollabID[ri], c("Box..")], sep=","), collapse=", ")
					loc_pos.temp <- paste(as.character(location[location[,1] %in% outdf$CollabID[ri], c("Position..")], sep=","), collapse=", ")

		# This is the new module of extracting cellularity information
		said<-cell.cal[cell.cal$sample_id %in% outdf$SAID[ri], c("sample_id")]

		if ( sum(nchar(said)) > 0 )
				{	
					if ( said == outdf$SAID[ri] ) { outdf$Cell.cal[ri]<-cell.cal[cell.cal$sample_id %in% outdf$SAID[ri], c("cellularity")] }

				}

		# from AB summary file
		cID<-cell.AB[cell.AB$Matching.Tumour.sample.name %in% outdf$CollabID[ri], c("Matching.Tumour.sample.name")]

		if ( sum(nchar(cID)) > 0 )
					{	
					if ( cID == outdf$CollabID[ri] ) 	{ 
										outdf$Cell.obs[ri]<-cell.AB[cell.AB$Matching.Tumour.sample.name %in% outdf$CollabID[ri], c("Cellularity")] 
										}
					}

		# from TB summary file
		patid<-gsub("_",".",cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf$CollabID[ri], c("TBID")])
		if ( sum(nchar(patid)) > 0 )
				{	
					if ( patid == outdf$CollabID[ri] ) { 
										outdf$Cell.obs[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf$CollabID[ri], c("cellularity")] 
										outdf$ER_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf$CollabID[ri], c("GenefuSubtype")] 
										outdf$HER2_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf$CollabID[ri], c("GenefuSubtype")] 
										outdf$CK5.6_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf$CollabID[ri], c("Pam50Subtype_overall")] 
										}
				}

		# from VBA summary file
		patid<-cell.VBA[cell.VBA$Sample.Name %in% outdf$CollabID[ri],c("Sample.Name")]
		if ( sum(nchar(patid)) > 0 )
				{	
					if ( patid == outdf$CollabID[ri] ) { 
										outdf$Cell.obs[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf$CollabID[ri], c("Cellularity")] 
										outdf$ER_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf$CollabID[ri], c("ER_clin")] 
										outdf$PR_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf$CollabID[ri], c("PR_clin")] 
										outdf$HER2_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf$CollabID[ri], c("HER2_clin")] 
										}
				}

		# from master_IHC file
		patid<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("Sample.ID")]
		if ( sum(nchar(patid)) > 0 )
				{	
					if ( patid == outdf$CollabID[ri] ) { 
										outdf$ER_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("ER_TMA")]
										outdf$PR_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("PR_TMA")]
										outdf$HER2_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("HER2_TMA")]
										outdf$EGFR_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("EGFR_TMA")]
										outdf$CK5.6_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf$CollabID[ri],c("CK5.6_TMA")]
										}
				}
### SPLIT HERE
		# from TNBC summary file
		patid<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("Sample.ID")]
		if ( sum(nchar(patid)) > 0 )
				{	
					if ( patid == outdf$CollabID[ri] ) 	{ 
										outdf$Cell.obs[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("cellularity")] 
										outdf$ER_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("ER_TMA")] 
										outdf$PR_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("PR_TMA")] 
										outdf$HER2_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("HER2_TMA")] 
										outdf$EGFR_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("EGFR_TMA")] 
										outdf$CK5.6_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf$CollabID[ri], c("CK5.6_TMA")] 
										}
				}


# This module stores the exact location(s) of the samples 
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

## END SEQ

# Remove incomplete cases which is no SA ID AND no Collaborator ID 
outdf<-unique(outdf[!(outdf$CollabID=="" & outdf$SAID==""),])

## CHECK
outdf[c("CollabID", "SAID","PROJ","Cell.cal","Cell.obs")]
####################################

# File is in this path on MOMAC14
# /Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files
tempfile="TNBC_Summary_int.csv"
write.table(outdf, file=tempfile, sep=",")

########################################
# Missing SA IDs from outdf
tnbc_add<-as.data.frame(unique(sort(tnbc_seq[!tnbc_seq$SA_ID %in% outdf$SAID,c("SA_ID")])), stringsAsFactors = FALSE)
rna_add1<-as.data.frame(unique(sort(rnaseq[!rnaseq$Sample %in% outdf$SAID, c(2)])), stringsAsFactors = FALSE)
rna_add2<-as.data.frame(unique(sort(rnaseq[!rnaseq$Sample.1 %in% outdf$SAID, c(3) ])), stringsAsFactors = FALSE)
rna_add<-as.data.frame(c("SA024","SA029","SA030","SA051","SA052","SA053","SA054","SA055","SA056","SA072"), stringsAsFactors = FALSE)

colnames(tnbc_add)[1]<-"SAID"
colnames(rna_add1)[1]<-"SAID"
colnames(rna_add2)[1]<-"SAID"
colnames(rna_add)[1]<-"SAID"
missingIDs<-rbind(tnbc_add, rna_add, rna_add1,rna_add2)

# make them into a list so that I can remove those that we already have data for
MID<-unlist(missingIDs, use.name=FALSE)
# These are the normals
NOR<-grep('N', MID, value=TRUE)
# Those which are not present in the outdf dataset
MIE<-MID[which(!MID %in% NOR)]
missingIDs<-MIE[MIE != " "]

rec=length(missingIDs)

# defining the data frame with placeholder data
# This second data frame contains all the missing SA IDs which supposedly belong to TNBC 
# but were for some reasons missed from the analysis above 
# this list is partially script generated

# defining the data frame with placeholder data
outdf2 <-  data.frame(CollabID = rep("", rec),
                  addCollabID = rep("", rec),
                         SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
#                    loc_rack = rep("NA", rec),
#                   loc_shelf = rep("NA", rec),
#                     loc_box = rep("NA", rec),
                      loc_pos = rep("NA", rec),
#                    rem_elut = rep("NA", rec),
#                    sam_type = rep("NA", rec),
                Tumour.Normal = rep("NA", rec),
                       ER_TMA = rep("NA", rec),
                       PR_TMA = rep("NA", rec),
                     HER2_TMA = rep("NA", rec),
                     EGFR_TMA = rep("NA", rec),
                    CK5.6_TMA = rep("NA", rec),
                    RNASeq_Mb = rep("NA", rec),
                    GSC_EXCAP = rep("NA", rec),
                    SOLID_WGS = rep("NA", rec),
                 Add_SOLIDLIB = rep("NA", rec),
                 SNP_ARRAY_ID = rep("NA", rec),
                  METABRIC_ID = rep("NA", rec),
                         PROJ = rep("NA", rec),
                       SOURCE = rep("NA", rec),
                    RNASEQLIB = rep("NA", rec),
                      WGSSLIB = rep("NA", rec),
                 NORM_WGSSLIB = rep("NA", rec),

                       Normal = rep("NA", rec),
                       Tumour = rep("NA", rec),
                RNASEQLIB.seq = rep("NA", rec),
                  WGSSLIB.seq = rep("NA", rec),
             NORM_WGSSLIB.seq = rep("NA", rec),

                     Cell.cal = rep("", rec),
                     Cell.obs = rep("", rec),
             stringsAsFactors = FALSE)

# This is the module where I sum up everything from the sheets
for (ri in seq(rec)) {

                # Using univ$someID to grep names from codes$Sample_ID (which is the Collaborator ID wanted) from codes
                collID<-gsub("[()]/", "_", as.character(codes[grep(univ$someID[ri], codes$Sample_ID) , c("Sample_ID")] ))

                if ( sum(nchar(collID)) > 0 )
                                {
                                        if ( is.na(unique(collID)) )    { next }
                                        if ( "collID" == univ$someID[ri] ) { outdf2$CollabID[ri] <- collID }
                                        else    {
                                                outdf2$addCollabID[ri] <- collID
                                                outdf2$CollabID[ri] <- codes[codes$Sample_ID %in% collID, c("Sample_ID")]
                                                outdf2$SAID[ri] <- codes[codes$Sample_ID %in% collID, c("SA_ID")]
                                                }
                                }


                # Used the standardized Collborator ID (Sample_ID) , c(1) <- col in codes to be output (Sample_ID)
                said<-codes[codes$Sample_ID %in% outdf2$CollabID[ri], c("Sample_ID")]

                if ( sum(nchar(said)) > 0 )
                                {
                                        if ( said == outdf2$CollabID[ri] ) { outdf2$SAID[ri] <- codes[codes$Sample_ID %in% outdf2$CollabID[ri], c("SA_ID")]  }
                                        if ( said == outdf2$CollabID[ri] ) { outdf2$addID[ri] <- codes[codes$Sample_ID %in% outdf2$CollabID[ri], c("Additional_ID")]  }
                                        if ( said == outdf2$CollabID[ri] ) { outdf2$add_SAID[ri] <- codes[codes$Sample_ID %in% outdf2$CollabID[ri], c("Additional_SA_IDs")]  }
                                }

                # For UPDATE: use WGSS library from TNBC-27 JIRA (103 genomes) and confirm that they are TNBC
                dnaid<-tnbc_seq[tnbc_seq$SA_ID %in% outdf2$SAID[ri], c("SA_ID")]

                if ( sum(nchar(dnaid)) > 0 )
                                {
                                        outdf2$WGSSLIB[ri] <- tnbc_seq[tnbc_seq$SA_ID %in% outdf2$SAID[ri], c("WGSS_LIB_ID")]
                                        outdf2$NORM_WGSSLIB[ri] <- tnbc_seq[tnbc_seq$SA_ID %in% paste(outdf2$SAID[ri],"N",sep=""), c("WGSS_LIB_ID")]
                                        outdf2$PROJ[ri] <- "TNBC_103"
                                }

                # SA ID [,2] and rnaseq$Sample, c("Sample") <- col in tnbc_seq to be output (SA ID)
                rnaid<-rnaseq[rnaseq$Sample %in% outdf2$SAID[ri], c("Sample")]

                if ( sum(nchar(rnaid)) == 0 )
                                {
                                        rnaid<-rnaseq[rnaseq$Sample.1 %in% outdf2$SAID[ri], c("Sample.1")]
                                        if ( sum(nchar(rnaid)) > 0 )
                                                {
                                                        outdf2$RNASEQLIB[ri] <- rnaseq[rnaseq$Sample.1 %in% outdf2$SAID[ri], c("Library")]
                                                        outdf2$PROJ[ri] <- rnaseq[rnaseq$Sample.1 %in% outdf2$SAID[ri], c("Description")]
                                                        rnaid<-NULL
                                                }
                                }
                if ( sum(nchar(rnaid)) > 0 )
                                {
                                        outdf2$RNASEQLIB[ri] <- rnaseq[rnaseq$Sample %in% outdf2$SAID[ri], c("Library")]
                                        outdf2$PROJ[ri] <- rnaseq[rnaseq$Sample %in% outdf2$SAID[ri], c("Description")]
                                }

                # Collab ID [,1] and location$Sample.Name, c("Sample.Name") <- col in location to be output (Collab ID)
                locid<-unique(location[location$Sample.Name %in% outdf2$CollabID[ri], c("Sample.Name")])

                if ( sum(nchar(locid)) > 0 )
                                {
                                        source <- unique(location[location$Sample.Name %in% outdf2$CollabID[ri], c("Source")])
                                        if ( sum(nchar(source)) > 1)
                                                {
                                                outdf2$SOURCE[ri] <- paste(unique(location[location$Sample.Name %in% outdf2$CollabID[ri], c("Source")]), collapse=",")
                                                }
                                }

## SPLIT HERE
                # Sample ID column match in location[,1], c(6-9) <- col in location to be output (Rack,Shelf,Box, Pos)
                sam_loc<-unique(location[location$Sample.Name %in%  outdf2$CollabID[ri], c(1)])

                if ( sum(nchar(sam_loc)) > 0 )
                                {
                                        loc_rack.temp <- paste(as.character(location[location$Sample.Name %in% outdf2$CollabID[ri], c("Rack..")],sep=","), collapse=", ")
                                        loc_shelf.temp <- paste(as.character(location[location[,1] %in% outdf2$CollabID[ri], c("Shelf..")], sep=","), collapse=", ")
                                        loc_box.temp <- paste(as.character(location[location[,1] %in% outdf2$CollabID[ri], c("Box..")], sep=","), collapse=", ")
                                        loc_pos.temp <- paste(as.character(location[location[,1] %in% outdf2$CollabID[ri], c("Position..")], sep=","), collapse=", ")

                # This is the new module of extracting cellularity information
                said<-cell.cal[cell.cal$sample_id %in% outdf2$SAID[ri], c("sample_id")]

                if ( sum(nchar(said)) > 0 )
                                {
                                        if ( said == outdf2$SAID[ri] ) { outdf2$Cell.cal[ri]<-cell.cal[cell.cal$sample_id %in% outdf2$SAID[ri], c("cellularity")] }

                                }

                # from AB summary file
                cID<-cell.AB[cell.AB$Matching.Tumour.sample.name %in% outdf2$CollabID[ri], c("Matching.Tumour.sample.name")]

                if ( sum(nchar(cID)) > 0 )
                                        {
                                        if ( cID == outdf2$CollabID[ri] )        {
                                                                                outdf2$Cell.obs[ri]<-cell.AB[cell.AB$Matching.Tumour.sample.name %in% outdf2$CollabID[ri], c("Cellularity")]
                                                                                }
                                        }

                # from TB summary file
                patid<-gsub("_",".",cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf2$CollabID[ri], c("TBID")])
                if ( sum(nchar(patid)) > 0 )
                                {
                                        if ( patid == outdf2$CollabID[ri] ) {
                                                                                outdf2$Cell.obs[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf2$CollabID[ri], c("cellularity")]
                                                                                outdf2$ER_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf2$CollabID[ri], c("GenefuSubtype")]
                                                                                outdf2$HER2_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf2$CollabID[ri], c("GenefuSubtype")]
                                                                                outdf2$CK5.6_TMA[ri]<-cell.TB[gsub("_",".", cell.TB$TBID) %in% outdf2$CollabID[ri], c("Pam50Subtype_overall")]
                                                                                }
                                }

                # from VBA summary file
                patid<-cell.VBA[cell.VBA$Sample.Name %in% outdf2$CollabID[ri],c("Sample.Name")]
                if ( sum(nchar(patid)) > 0 )
                                {
                                        if ( patid == outdf2$CollabID[ri] ) {
                                                                                outdf2$Cell.obs[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf2$CollabID[ri], c("Cellularity")]
                                                                                outdf2$ER_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf2$CollabID[ri], c("ER_clin")]
                                                                                outdf2$PR_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf2$CollabID[ri], c("PR_clin")]
                                                                                outdf2$HER2_TMA[ri]<-cell.VBA[cell.VBA$Sample.Name %in% outdf2$CollabID[ri], c("HER2_clin")]
                                                                                }
                                }

                # from master_IHC file
                patid<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("Sample.ID")]
                if ( sum(nchar(patid)) > 0 )
                                {
                                        if ( patid == outdf2$CollabID[ri] ) {
                                                                                outdf2$ER_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("ER_TMA")]
                                                                                outdf2$PR_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("PR_TMA")]
                                                                                outdf2$HER2_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("HER2_TMA")]
                                                                                outdf2$EGFR_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("EGFR_TMA")]
                                                                                outdf2$CK5.6_TMA[ri]<-master_IHC[master_IHC$Sample.ID %in% outdf2$CollabID[ri],c("CK5.6_TMA")]
                                                                                }
                                }
### SPLIT HERE
                # from TNBC summary file
                patid<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("Sample.ID")]
                if ( sum(nchar(patid)) > 0 )
                                {
                                        if ( patid == outdf2$CollabID[ri] )      {
                                                                                outdf2$Cell.obs[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("cellularity")]
                                                                                outdf2$ER_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("ER_TMA")]
                                                                                outdf2$PR_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("PR_TMA")]
                                                                                outdf2$HER2_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("HER2_TMA")]
                                                                                outdf2$EGFR_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("EGFR_TMA")]
                                                                                outdf2$CK5.6_TMA[ri]<-cell.TNBC[cell.TNBC$Sample.ID %in% outdf2$CollabID[ri], c("CK5.6_TMA")]
                                                                                }
                                }


# This module stores the exact location(s) of the samples
#                                       if (length(loc_rack.temp) > 0 ) { outdf2$loc_rack[ri] <- loc_rack.temp }
#                                       if (length(loc_shelf.temp) > 0 ) { outdf2$loc_shelf[ri] <- loc_shelf.temp }
#                                       if (length(loc_box.temp) > 0 ) { outdf2$loc_box[ri] <- loc_box.temp }
#                                       if (length(loc_pos.temp) > 0 ) { outdf2$loc_pos[ri] <- loc_pos.temp }

# This modules says that we have samples only (location not given) We use position to indicate that we have physical sample
# AS long as it has a shelf, box or position reference, we assume that we have samples
                                        if (length(loc_shelf.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                        if (length(loc_box.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                        if (length(loc_pos.temp) > 0 ) { outdf2$loc_pos[ri] <- "Sample" }
                                }
                        }

## END SEQ

# Remove incomplete cases which is no SA ID AND no Collaborator ID
outdf2<-unique(outdf2[!(outdf2$CollabID=="" & outdf2$SAID==""),])



#Join data sets
outdf3_1<-join(outdf, outdf2, type = "full")
outdf3<-outdf3_1[!duplicated(outdf3_1),]

#outdf3<-mutate(output, WGSS=ifelse(output$WGSSLIB=="NA", "None", "Seq'd"))
#outdf3<-mutate(output, WTSS=ifelse(output$RNASEQLIB=="NA", "None", "Seq'd"))

############# Use the completed data to cross-reference Sequencing Master file 2017 ############
# This is the master sequencing file
seqfile="Aparicio_Library_master_list__January132017.xlsx"
seq<-readWorksheetFromFile(file=seqfile,sheet = c("master page"), header=TRUE)

#### NOTE CHNAGED OUTDF STRUCTURE ############
#### ADD AT END AS MERGING BY ROWS ##########

rec=nrow(outdf3)
outdf4 <- data.frame(CollabID = rep("", rec),
                  addCollabID = rep("", rec),
                         SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
                      loc_pos = rep("", rec),
                Tumour.Normal = rep("", rec),
                       ER_TMA = rep("", rec),
                       PR_TMA = rep("", rec),
                     HER2_TMA = rep("", rec),
                     EGFR_TMA = rep("", rec),
                    CK5.6_TMA = rep("", rec),
                    RNASeq_Mb = rep("", rec),
                    GSC_EXCAP = rep("", rec),
                    SOLID_WGS = rep("", rec),
                 Add_SOLIDLIB = rep("", rec),
                 SNP_ARRAY_ID = rep("", rec),
                  METABRIC_ID = rep("", rec),
                         PROJ = rep("", rec),
                       SOURCE = rep("", rec),
                    RNASEQLIB = rep("", rec),
                      WGSSLIB = rep("", rec),
                 NORM_WGSSLIB = rep("", rec),

                       Normal = rep("", rec),
                       Tumour = rep("", rec),
                RNASEQLIB.seq = rep("", rec),
                  WGSSLIB.seq = rep("", rec),
             NORM_WGSSLIB.seq = rep("", rec),

                     Cell.cal = rep("", rec),
                     Cell.obs = rep("", rec),
             stringsAsFactors = FALSE)

# I assmume that everything that gets SA ID for seq 
# Therefore it makes sense that I only use the SA ID to probe the spreadsheets

for (ri in seq(rec)) {

		# First get all the SA ID and if SA ID is not there label not suitable	 or missing
		test_said<-outdf3$SAID[ri]

                if ( sum(nchar(test_said)) > 1 )
				{
				if ( test_said == "not suitable for sequencing")  	{ test_said<-"not suitable" 
											outdf4$SAID[ri]<-test_said
											next
											} else { test_said<-"not assigned" 
												outdf4$SAID[ri]<-test_said
												}

				# Then transfer entire line from outdf3 to outdf4
        	        	outdf4[ri,]<-outdf3[ri,]

				collID<-gsub(" ","", master_summary[master_summary$SA.ID %in% outdf3$SAID[ri], c("Sample.ID")])

                if ( sum(nchar(collID)) > 0 & length(collID) == 1 )
                                {
				 # Fix Collaborator ID to remove spaces 
                                 outdf4$CollabID[ri] <- collID
				
                                }

		# Getting the IHC results from master_IHC speadsheet
                ihcid<-master_IHC[master_IHC$SA_id %in% outdf4$SAID[ri], c("SA_id")]

                if ( sum(nchar(length(ihcid)) > 0 & length(ihcid) == 1))
                                {
				# Get the normals and put them in the respective columns
				outdf4[ri,c("Normal")]<-unique(master_IHC[master_IHC$SA_id %in% paste(outdf4$SAID[ri],"N",sep=""), c("Tumour.Normal")])
				outdf4[ri,c("NORM_WGSSLIB.seq","SOLID_WGS","Add_SOLIDLIB")]<-unique(master_IHC[master_IHC$SA_id %in% paste(outdf4$SAID[ri],"N",sep=""), 
					c("GSC_WGSS","SOLID_WGS","Addditional.SOLiD.libraries")])

				# get the tumours and put them in the respective columns
				outdf4[ri,c("Tumour","ER_TMA","PR_TMA","HER2_TMA","EGFR_TMA","CK5.6_TMA")]<-master_IHC[master_IHC$SA_id %in% outdf4$SAID[ri], 
					c("Tumour.Normal","ER_TMA","PR_TMA","HER2_TMA","EGFR_TMA","CK5.6_TMA")]
				outdf4[ri,c("RNASEQLIB.seq","RNASeq_Mb","WGSSLIB.seq","GSC_EXCAP","SOLID_WGS","Add_SOLIDLIB","PROJ","SOURCE")]<-master_IHC[master_IHC$SA_id %in% outdf4$SAID[ri], 
					c("GSC_RNASeq","RNASeq_Mb","GSC_WGSS","GSC_EXCAP","SOLID_WGS","Addditional.SOLiD.libraries","PROJECT","Source")]
                                }


				 }

			
			}

### END HERE FOR TEST
## CHECK
colSums(is.na(outdf4))
## Incorporate Paul Savage / Morag Park TNBC cases with SA IDs
parkmetafile="GCRC PDX metadata Overview 17-02-28 (Aparicio).xlsx"
parkmeta<-readWorksheetFromFile(file=parkmetafile,sheet = c("Sheet1"), header=TRUE)

parkmeta[parkmeta=="<NA>"]=NA

parkseqfile="Paul Savage SA IDs.xlsx"
parkseq<-readWorksheetFromFile(file=parkseqfile,sheet = c("Sheet1"), header=TRUE)

parkseq[parkseq=="<NA>"]=NA

rec=nrow(parkseq)
outdf5 <- data.frame(CollabID = rep("", rec),
                  addCollabID = rep("", rec),
                         SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
                      loc_pos = rep("", rec),
                Tumour.Normal = rep("", rec),
                       ER_TMA = rep("", rec),
                       PR_TMA = rep("", rec),
                     HER2_TMA = rep("", rec),
                     EGFR_TMA = rep("", rec),
                    CK5.6_TMA = rep("", rec),
                    RNASeq_Mb = rep("", rec),
                    GSC_EXCAP = rep("", rec),
                    SOLID_WGS = rep("", rec),
                 Add_SOLIDLIB = rep("", rec),
                 SNP_ARRAY_ID = rep("", rec),
                  METABRIC_ID = rep("", rec),
                         PROJ = rep("Savage", rec),
                       SOURCE = rep("Xenos", rec),
                    RNASEQLIB = rep("", rec),
                      WGSSLIB = rep("", rec),
                 NORM_WGSSLIB = rep("", rec),

                       Normal = rep("", rec),
                       Tumour = rep("", rec),
                RNASEQLIB.seq = rep("", rec),
                  WGSSLIB.seq = rep("", rec),
             NORM_WGSSLIB.seq = rep("", rec),

             stringsAsFactors = FALSE)

for (rj in seq(rec)) {

#			outdf5$CollabID[rj]<-substr(parkseq$PDX.ID[rj],1,8)
			outdf5$CollabID[rj]<-parkseq$PDX.ID[rj]

                        # get the tumours and put them in the respective columns
                        outdf5[rj,c("ER_TMA","PR_TMA","HER2_TMA","Tumour.Normal")]<-parkmeta[substr(parkmeta$PDX.ID,1,8) %in% outdf5$CollabID[rj],
                                       c("ER.","PR.","HER2.","Histology")]

                        outdf5[rj,c("SAID","WGSSLIB.seq","RNASEQLIB.seq","NORM_WGSSLIB.seq","RNASEQLIB","WGSSLIB")]<-parkseq[parkseq$PDX.ID %in% outdf5$CollabID[rj],
                                       c("SA.ID","WGS.PDX","RNAseq.PDX","WGS.Germline","RNAseq.Patient","WGS.Patient")]



                        }



outdf6<-join(outdf4, outdf5, type = "full")

# Incorporate Peter's Treatment Naive xenos
xeno_seqfile="xeno_TNBC_for_seq.xlsx"
xenoseq<-readWorksheetFromFile(file=xeno_seqfile,sheet = c("TNBC_TxNaive"), header=TRUE)

norm_xenofile="Prep Matched normal samples for DNA prep_BC_MS.xlsx"
normxeno<-readWorksheetFromFile(file=norm_xenofile,sheet = c("Sheet1"), header=TRUE, startRow = 2)

# Get all the CollabIDs
LL1<-as.character(unique(xenoseq$Patient_ID))
LL2<-c(LL1,as.character(unique(normxeno$Number)))
univ<-unique(LL2)

rec=length(univ)
outdf7 <- data.frame(CollabID = rep("", rec),
                  addCollabID = rep("", rec),
                         SAID = rep("", rec),
                        addID = rep("", rec),
                     add_SAID = rep("", rec),
                      loc_pos = rep("", rec),
                Tumour.Normal = rep("", rec),
                       ER_TMA = rep("", rec),
                       PR_TMA = rep("", rec),
                     HER2_TMA = rep("", rec),
                     EGFR_TMA = rep("", rec),
                    CK5.6_TMA = rep("", rec),
                    RNASeq_Mb = rep("", rec),
                    GSC_EXCAP = rep("", rec),
                    SOLID_WGS = rep("", rec),
                 Add_SOLIDLIB = rep("", rec),
                 SNP_ARRAY_ID = rep("", rec),
                  METABRIC_ID = rep("", rec),
                         PROJ = rep("TxNaiveXenos", rec),
                       SOURCE = rep("Peter", rec),
                    RNASEQLIB = rep("", rec),
                      WGSSLIB = rep("", rec),
                 NORM_WGSSLIB = rep("", rec),

                       Normal = rep("", rec),
                       Tumour = rep("", rec),
                RNASEQLIB.seq = rep("", rec),
                  WGSSLIB.seq = rep("", rec),
             NORM_WGSSLIB.seq = rep("", rec),

             stringsAsFactors = FALSE)

for (rj in seq(rec)) {

			outdf7$CollabID[rj]<-univ[rj]

                        # get the tumours and put them in the respective columns

			coid<-normxeno[normxeno$Number %in% outdf7$CollabID[rj],c("Number")]

                	if ( sum(nchar(coid)) > 1 & length(coid) == 1 )
                                {

                        outdf7[rj,c("Normal","addID")]<-normxeno[normxeno$Number %in% outdf7$CollabID[rj],
                                       c("Type..salive..buffy.coat..etc.","Vial.Label")]
				}

			said<-xenoseq[xenoseq$Patient_ID %in% outdf7$CollabID[rj],c("SA_ID")]

                	if ( sum(nchar(said)) > 1 & length(said) == 1 )
                                {
                        outdf7[rj,c("SAID","Tumour.Normal")]<-xenoseq[xenoseq$Patient_ID %in% outdf7$CollabID[rj],
                                      c("SA_ID","Pathology")]
				}



                        }



outdf8<-join(outdf6, outdf7, type = "full")
# File is in this path on MOMAC14
# /Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files
tempfile="TNBC_Summary_outdf8.csv"
write.table(outdf8, file=tempfile, sep=",")
################################


############# SUMMARY ################

# Need RNA-seq if have WGSS but no RNA-seq


# SPlit into columns for RNA

RNASeq_Need <-outdf8[,c("CollabID", "SAID", "loc_pos", "RNASeq_Mb","RNASEQLIB","RNASEQLIB.seq")]
DNASeq_Need <-outdf8[,c("CollabID", "SAID", "loc_pos", "GSC_EXCAP","SOLID_WGS", "WGSSLIB","WGSSLIB.seq")]
Norm_DNASeq_Need <-outdf8[,c("CollabID", "SAID", "loc_pos", "NORM_WGSSLIB","NORM_WGSSLIB.seq")]
IHC_Need <- outdf8[,c("CollabID", "SAID", "loc_pos", "ER_TMA","PR_TMA", "HER2_TMA")]
Histo_Need <- outdf8[,c("CollabID", "SAID", "loc_pos", "Cell.cal","Cell.obs")]
Complete <- outdf8[,c("CollabID", "SAID", "RNASEQLIB", "RNASEQLIB.seq","WGSSLIB", "WGSSLIB.seq", "NORM_WGSSLIB","NORM_WGSSLIB.seq")]

# File is in this path on MOMAC14
# /Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files
outfile="TNBC_Summary_by_RNASeq_Need.csv"
write.table(RNASeq_Need, file=outfile, sep=",")

outfile="TNBC_Summary_by_DNASeq_Need.csv"
write.table(DNASeq_Need, file=outfile, sep=",")

outfile="TNBC_Summary_by_Norm_DNASeq_Need.csv"
write.table(Norm_DNASeq_Need, file=outfile, sep=",")

outfile="TNBC_Summary_by_IHC_Need.csv"
write.table(IHC_Need, file=outfile, sep=",")

outfile="TNBC_Summary_by_Histo_Need.csv"
write.table(Histo_Need, file=outfile, sep=",")

outfile="TNBC_Summary_by_CompleteCases.csv"
write.table(Complete, file=outfile, sep=",")





############

#outdf4$addID<-paste(outdf4$addID,outdf4$add_SAID,sep="")

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
outfile="TNBC_Summary_by_all_SA_onMAC_4.csv"
write.table(output, file=outfile, sep=",")

# Format
# cat (outputfile) | perl -p -e 's/\r\n|\n|\r/\|\n/g' | tr "," "|" | sed 's/^/\|/' | sed 's/\"\"/\ /g' | tr -d '"'

