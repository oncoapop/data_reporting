# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library('BSgenome.Hsapiens.UCSC.hg19')

# Please change a total of 4 (four) positions where the name of the sample appears

rearr <- read.csv(file="/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/SVs_selected.csv",  stringsAsFactors = FALSE)

# For same chromosome deletions, get the upstream and downstream sequences of the breakpoints

# cluster_id,sample_gene_name1_gene_location_1_gene_name2_gene_location2



outdf <- data.frame(Id = rep("id", nrow(rearr)),
                    Des = rep("des", nrow(rearr)),
                    Ups = rep(0, nrow(rearr)),
                    Dps = rep(0, nrow(rearr)),
                    Upseq = rep("", nrow(rearr)),                    
                    Downseq = rep("", nrow(rearr)),                    
                    Myseq = rep("", nrow(rearr)), 
                    WTseq = rep("", nrow(rearr)),                   
                    Chimera = rep("", nrow(rearr)),                    
                    stringsAsFactors = FALSE)
                    
# Include design space of n bp on either side of breakpoint                     
offset <- 400

# Add primer3 markers to ensure amplification of breakpoints, n bp either side 
protect <- 50
                     
for (ri in seq(nrow(rearr))) {
  
  id <- rearr$cluster_id[ri]
  des <- paste(rearr$Sample[ri],rearr$gene_name1[ri],rearr$gene_location1[ri],rearr$gene_name2[ri],rearr$gene_location2[ri], sep="_")
  chr <- paste("chr",rearr$chromosome_1[ri], sep="")
  
  cr1 <- rearr$chromosome_1[ri]
  cr2 <- rearr$chromosome_2[ri]
  bk1 <- rearr$break_1[ri]
  bk2 <- rearr$break_2[ri]
  
  if (cr1==cr2) { if (bk1 < bk2) ups <- bk1
  	dps <- bk2 }
  	
  if (cr1==cr2) { if (bk1 > bk2) ups <- bk2
    dps <- bk1 }
                
  
  chimera <- rearr$sequence[ri]

  upseq <- paste(getSeq(Hsapiens,chr,ups-offset,ups-protect),"<",
              getSeq(Hsapiens,chr,ups-protect+1,ups+protect),">",
              getSeq(Hsapiens,chr,ups+protect+1,ups+offset),
              sep='')
              
  downseq <- paste(getSeq(Hsapiens,chr,dps-offset,dps-protect),"<",
              getSeq(Hsapiens,chr,dps-protect+1,dps+protect),">",
              getSeq(Hsapiens,chr,dps+protect+1,dps+offset),
              sep='')
              
  myseq <- paste(getSeq(Hsapiens,chr,ups-offset,ups-protect),"<",
              getSeq(Hsapiens,chr,ups-protect+1,ups),
              getSeq(Hsapiens,chr,dps,dps+protect),">"
              getSeq(Hsapient,chr,dps+protect+1,dps+offset),
              sep='')
              
  wtseq <- getSeq(Hsapiens,chr,ups-offset,dps+offset)
              
            
# chromosome positions needs to be in format "chrn", where n=1-23 or X,Y 
             
             
  outdf$Id[ri] <- id
  outdf$Des[ri] <- des
  outdf$Ups[ri] <- ups
  outdf$Dps[ri] <- dps
  outdf$Upseq[ri] <- upseq
  outdf$Downseq[ri] <- downseq
  outdf$Myseq[ri] <- myseq
  outdf$WTseq[ri] <- wtseq
  outdf$Chimera[ri] <- chimera
  }

write.csv(outdf, file = "/share/lustre/backup/dyap/Projects/Tumour_Evol/Rearrangements/SVs_selected_designspace.csv")



