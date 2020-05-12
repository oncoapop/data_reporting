# These commands must be specifed in order for this script to work
#source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library('BSgenome.Hsapiens.UCSC.hg19')
biocLite("SNPlocs.Hsapiens.dbSNP.20120608")

getSeqHardMasked <-
  function(BSg,GR,maskList,letter) {
### PURPOSE: return list of DNAString sequences extracted from the
### BSgenome <BSg> corresponding to each location in GenomicRange
### <GR>, and masked with <letter> according to the masks named in
### <maskList> (which are encoded following BSParams convention).
###
### USE CASE - write fasta file of hard masked regions, using
###            RepeatMasker (RM) and Tandem Repeat Finder (TRF):
###
### GR <- GRanges('chr2L',IRanges(c(1,100),c(15,125)))
### writeFASTA(getSeqHardMasked(BSgenome, GR, c(RM=TRUE,TRF=TRUE), "N")
###            ,"myExtractForGR.fa"
###            ,paste(seqnames(GR),start(GR),end(GR),strand(GR),sep=':')
###            )
###
### NB: The implementation was coded 'pay(ing) attention to memory
### management' following suggestions from Herve in:
### https://stat.ethz.ch/pipermail/bioconductor/2011-March/038143.html.
### In particular, the inidividual chromosomes and their
### subseq(uences) should be garbage collectable after the function
### exits and they go out of scope, (although the chromosomes _are_
### all simultaneously loaded which I think is unavoidable if the
### results are to preserve the arbitrary order of GR).
###
### NB: My initial implementation FAILed as it used bsapply & BSParams
### whose FUN can not 'know' the name of the sequence (which was
### needed to know which subseqs to extract).
    ']]' <-
      ## utility to subscript b by a.
      function(a,b) b[[a]]
    Vsubseq <-
      ## vectorized version of subseq.
      Vectorize(subseq)
    VinjectHardMask <-
      ## vectorized version of injectHardMask.
      Vectorize(injectHardMask)
    activeMask <-
      ## A logical vector indicating whether each mask should be ON or
      ## OFF to be applied to each chromosome in BSg.
      masknames(BSg) %in% names(maskList[which(maskList)])
    BSg_seqList <-
      ## BSg as a list of named MaskedDNAString (one per chromosome)...
      sapply(seqnames(BSg),']]',BSg)
    BSg_seqList <-
      ## ... with the masks for each chromosome activated.
      sapply(BSg_seqList,function(x) {active(masks(x)) <- activeMask;x})
    GR_seq <-
      ## the MaskedDNAString corresponding to GR
      sapply(as.character(seqnames(GR)),']]',BSg_seqList)
    VinjectHardMask(Vsubseq(GR_seq,start(GR),end(GR)),letter=letter)
}



snvdf <- read.table(file = "/home/dyap/Projects/Single_Cell/positions/VBA0038_28plus7_SNVs-hg19.txt", stringsAsFactors = FALSE)

# defining the data frame with placeholder data
outdf <- data.frame(Chr = rep("", nrow(snvdf)),
                     Pos = rep(0, nrow(snvdf)),
                     WT = rep("N", nrow(snvdf)),
                     SeqUp = rep("", nrow(snvdf)),
                     SeqDown = rep("", nrow(snvdf)),
                     stringsAsFactors = FALSE)

#inject snps into the sequence
SnpHsapiens <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20120608")
                     
# This is the offset on each side of the SNV that is to be returned
offset <- 100
      
# Note that the output up and downstream sequences will not contain the SNV nucleotide     

for (ri in seq(nrow(snvdf))) {
  		chr <- snvdf[ri,1]
  		position <- snvdf[ri,2]
  		
# 		WT <- paste(getSeq(Hsapiens,chr,position,position), sep='')
  		GR <- GRanges(chr,IRanges(position,position))
		tmp <- getSeqHardMasked(Hsapiens, GR, c(RM=TRUE,TRF=TRUE), "N")
		WT <- paste(tmp$chr,sep='')
		  		
#   		seqUp <- paste(getSeq(Hsapiens,chr,position-offset,position-1), sep='')
  		GR <- GRanges(chr,IRanges(position-offset,position-1))
		tmp <- getSeqHardMasked(Hsapiens, GR, c(RM=TRUE,TRF=TRUE), "N")
		seqUp <- paste(tmp$chr,sep='')

#              seqDown <- paste(getSeq(Hsapiens,chr,position+1,position+offset), sep='')
                GR <- GRanges(chr,IRanges(position+1,position+offset))
		tmp <- getSeqHardMasked(Hsapiens, GR, c(RM=TRUE,TRF=TRUE), "N")
		seqDown <- paste(tmp$chr,sep='')

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


