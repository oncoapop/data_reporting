library(XLConnect)
library(plyr)

# This is where all the relevant excel spreadsheets may lie on MONCO
dir="/Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files/"
setwd(dir)

# This is the file where all that have been sequenced before and have data on Lustre (Archive)
archive="file_Luster_Jun17"
all_saids<-read.table(archive,header=TRUE,stringsAsFactors = FALSE,sep="\t")
saids<-data.frame(as.character(unique(all_saids$sid)), stringsAsFactors = FALSE)
colnames(saids)[1]<-"SA_IDs"
# count(all_saids$sid)

# This is the RNA-seq file
rnaseqfile="RNA-Seq libraries_Aparicio 2012May16 edited.xls"
rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

# This is the 103 genome data from Tyler
#genomefile="TNBC_WGSS_103.csv"
genomefile="TNBC_27_matches.csv"

# This is how many records there might be
rec=nrow(saids)

# wgss files
wgss<-read.csv(genomefile,header=TRUE,stringsAsFactors = FALSE)

# defining the data frame with placeholder data
outdf <- data.frame(       SAID = rep("NA", rec),
		       add_SAID = rep("", rec),
		       CollabID = rep("", rec),
		   add_CollabID = rep("", rec),
			   WGSS = rep("", rec),
	   		   WTSS = rep("", rec),
			Project = rep("", rec),
                       WGSSLIB = rep("", rec),
                      RNASEQLIB = rep("", rec),
#	              RNASEQLIB2 = rep("", rec),
                     stringsAsFactors = FALSE)



for (ri in seq(rec)) {
                outdf$SAID[ri] <- saids[ri,1]

		dnaid<-wgss[wgss[,3] %in% outdf$SAID[ri],c(3)]

		if ( length(dnaid) > 0 )
				{	
					if (dnaid == outdf$SAID[ri] ) { outdf$WGSSLIB[ri] <- wgss[wgss[,3] %in% saids[ri,1],c(7)]
									outdf$add_SAID[ri] <- wgss[wgss[,3] %in% saids[ri,1],c(5)]
									outdf$CollabID[ri] <- wgss[wgss[,3] %in% saids[ri,1],c(2)]
									outdf$add_CollabID[ri] <- wgss[wgss[,3] %in% saids[ri,1],c(4)]
#                      							outdf$RNASEQLIB2[ri] <- wgss[wgss[,3] %in% saids[ri,1],c(6)]

									 } 
				}

		rnaid<-rnaseq[rnaseq[,2] %in% outdf$SAID[ri],c(2)]

		if ( length(rnaid) > 0 )
				{	
					if (rnaid == outdf$SAID[ri] ) { outdf$RNASEQLIB[ri] <- rnaseq[rnaseq[,2] %in% saids[ri,1],c(1)] }
					if (rnaid == outdf$SAID[ri] ) { outdf$Project[ri] <- rnaseq[rnaseq[,2] %in% saids[ri,1],c(5)] }  
					
				}

 
			}
# use to total all_saids file
# check if there are RNA-seq if yes, annotate all_saids
# do the same from WGSS.


#MS<-mutate(MS, Interactors=ifelse(MS$Gene.name %in% Int$V1, "Known Interactors (34/97)", "Unknown"))
outdf<-mutate(outdf, WGSS=ifelse(nchar(outdf$WGSSLIB) > 1 & !is.na(outdf$WGSSLIB), "Yes", "None"))
outdf<-mutate(outdf, WTSS=ifelse(nchar(outdf$RNASEQLIB) > 1 & !is.na(outdf$RNASEQLIB), "Yes", "None"))


