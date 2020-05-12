# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library('BSgenome.Hsapiens.UCSC.hg19')
snvdf <- read.table(file = "/home/dyap/Projects/Single_Cell/positions/VBA0038_28plus7_SNVs-hg19.txt", stringsAsFactors = FALSE)

# defining the data frame with placeholder data
outdf <- data.frame(Chr = rep("", nrow(snvdf)),
                     Pos = rep(0, nrow(snvdf)),
                     WT = rep("N", nrow(snvdf)),
                     SeqUp = rep("", nrow(snvdf)),
                     SeqDown = rep("", nrow(snvdf)),
                     stringsAsFactors = FALSE)
                     
# This is the offset on each side of the SNV that is to be returned
offset <- 5
      
# Note that the output up and downstream sequences will not contain the SNV nucleotide     

		for (ri in seq(nrow(snvdf))) {
  		chr <- snvdf[ri,1]
  		position <- snvdf[ri,2]
  		
  		WT <- paste(getSeq(Hsapiens,chr,position,position),
  				sep='')
  		
  		seqUp <- paste(getSeq(Hsapiens,chr,position-offset,position-1),
              sep='')

  		seqDown <- paste(getSeq(Hsapiens,chr,position+1,position+offset),
              sep='')
              
  		outdf$Chr[ri] <- chr
  		outdf$Pos[ri] <- position
  		outdf$WT[ri] <- WT
  		outdf$SeqUp[ri] <- seqUp
  		outdf$SeqDown[ri] <- seqDown
}	



write.csv(outdf, file = "/home/dyap/Projects/Single_Cell/positions/VBA0038_35_SNV_5.txt")

# file 2
snvdf <- read.table(file = "/home/dyap/Projects/Single_Cell/positions/13_CoNAn-SNV-hg19.txt", stringsAsFactors = FALSE)

# defining the data frame with placeholder data
outdf <- data.frame(Chr = rep("", nrow(snvdf)),
                     Pos = rep(0, nrow(snvdf)),
                     WT = rep("N", nrow(snvdf)),
                     SeqUp = rep("", nrow(snvdf)),
                     SeqDown = rep("", nrow(snvdf)),
                     stringsAsFactors = FALSE)
                     
# This is the offset on each side of the SNV that is to be returned
offset <- 5
      
# Note that the output up and downstream sequences will not contain the SNV nucleotide     

		for (ri in seq(nrow(snvdf))) {
  		chr <- snvdf[ri,1]
  		position <- snvdf[ri,2]
  		
  		WT <- paste(getSeq(Hsapiens,chr,position,position),
  				sep='')
  		
  		seqUp <- paste(getSeq(Hsapiens,chr,position-offset,position-1),
              sep='')

  		seqDown <- paste(getSeq(Hsapiens,chr,position+1,position+offset),
              sep='')
              
  		outdf$Chr[ri] <- chr
  		outdf$Pos[ri] <- position
  		outdf$WT[ri] <- WT
  		outdf$SeqUp[ri] <- seqUp
  		outdf$SeqDown[ri] <- seqDown
  		  									}	



write.csv(outdf, file = "/home/dyap/Projects/Single_Cell/positions/VBA0038_13_SNV_5.txt")


