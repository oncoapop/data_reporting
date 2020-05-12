library(XLConnect)
library(plyr)

# This is where all the relevant excel spreadsheets may lie on MONCO
dir="/Volumes/Monco/Aparicio Lab - Projects/BreastCancerOutcomesUnitRelated/TNBC files/"
setwd(dir)

# This is the file where all that have been sequenced before and have data on Lustre (Archive)
archive="file_Luster_Jun17"
all_saids<-read.table(archive,header=TRUE,stringsAsFactors = FALSE,sep="\t")
saids<-unique(all_saids$sid)
count(all_saids$sid)

# This is the RNA-seq file
rnaseqfile="RNA-Seq libraries_Aparicio 2012May16.xls"
rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

# This is the 103 genome data from Tyler
genomefile="TNBC_WGSS_103.csv"

# This is how many records there might be
rec=length(saids)

# wgss files
wgss<-read.csv(genomefile,header=FALSE,stringsAsFactors = FALSE)

# defining the data frame with placeholder data
outdf <- data.frame(SAID = rep("NA", rec),
                     WGSSPATH = rep("NA", rec),
                     RNASEQLIB = rep("NA", rec),
                     stringsAsFactors = FALSE)

for (ri in seq(rec)) {
                outdf$SAID[ri] <- saids[ri]
		outdf$WGSSPATH[ri] <- wgss[ri,2]

		dnaid<-wgss[wgss[,1] %in% outdf$SAID[ri],c(1)]

		if ( length(dnaid) > 0 )
				{	
					if (dnaid == outdf$SAID[ri] ) { outdf$WGSSPATH[ri] <- wgss[wgss[,1] %in% outdf$SAID[ri],c(2)] } 
					
				}

		rnaid<-rnaseq[rnaseq[,2] %in% outdf$SAID[ri],c(2)]

		if ( length(rnaid) > 0 )
				{	
					if (rnaid == outdf$SAID[ri] ) { outdf$RNASEQLIB[ri] <- rnaseq[rnaseq[,2] %in% outdf$SAID[ri],c(1)] } 
					
				}

 
			}
