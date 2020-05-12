# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library('BSgenome.Hsapiens.UCSC.hg19')
snvdf <- read.table(file = "/share/lustre/backup/dyap/Projects/Single_Cell/positions/SNV/VBA0038_pos.txt", 
stringsAsFactors = FALSE)

outdf <- data.frame(Chr = rep("", nrow(snvdf)),
                     Pos = rep(0, nrow(snvdf)),
                     Seq = rep("", nrow(snvdf)),
                     stringsAsFactors = FALSE)
                     
offset <- 150
                     
for (ri in seq(nrow(snvdf))) {
  chr <- strsplit(snvdf[ri,1], split=":")[[1]][1]
  position <- as.numeric(strsplit(strsplit(snvdf[ri,1], split=":")[[1]][2], split = "-")[[1]][1])

  seq <- paste(getSeq(Hsapiens,chr,position-offset,position),
              getSeq(Hsapiens,chr,position+1,position+offset),
              sep='')
  outdf$Chr[ri] <- chr
  outdf$Pos[ri] <- position
  outdf$Seq[ri] <- seq
  }
		  									
write.csv(outdf, file = 
"/share/lustre/backup/dyap/Projects/Single_Cell/positions/VBA0038_positions.csv")


