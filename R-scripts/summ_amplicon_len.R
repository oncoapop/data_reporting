# Script to plot amplicon lengths
# Input is Amplen.ProjID (col 5) of combine order file
# Output is Amplicon length plot

primers<-read.table(file="/home/dyap/Projects/TNBC/primer3/TNBC-combined_order.txt", sep=",", 
header=TRUE)

amplicons <- data.frame ( Sam= rep("", nrow(primers)),
                         Amplen = rep(0, nrow(primers)),
                     stringsAsFactors = FALSE)

for (ri in seq(nrow(primers)))
	{
	amplicons$Sam[ri] <- as.character(primers[ri,1])
	amplicons$Amplen[ri] <- as.numeric(gsub("bp", "", strsplit(as.character(primers[ri,5]), 
split=":")[[1]][1]))
	}


pdf("TNBC_Amplicons.pdf", width=6, height=10)

boxplot(Amplen ~ Sam, data=amplicons, main="Size distribution of TNBC Amplicons for MiSeq 250x2bp 
run", ylab="Amplicon size (bp)", xlab="Sample ID (10 total)", cex.axis=0.8, las=2)

dev.off()

