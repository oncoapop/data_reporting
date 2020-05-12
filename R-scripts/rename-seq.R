### Rename gene names
### Note:  RBP_CG_Assay G-CUSTOM-176687.csv saved as .csv from Excel, then edited to get rid of two top header rows and bottom cruft
###        Can do this in excel (but double check in a Real Editor such as Emacs) CAVEAT
indf <- read.table(file = "RBP_CG_Assay G-CUSTOM-176687.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
indfnoc <- indf[!grepl("ON-TARGET", indf$Gene.Symbol), ]
indfnoc$Gene.Symbol
gsrle <- rle(indfnoc$Gene.Symbol)
                   
table(gsrle$values, useNA = "always")
table(table(gsrle$values))
seqnos <- unlist(lapply(gsrle$lengths, function(x) seq(x)))
gnms <- unlist(lapply(seq(along = gsrle$values), function(x) {rep(gsrle$values[x], gsrle$lengths[x]) } ) )
gnmsnos <- paste(gnms, seqnos, sep = "-")

indf$Gene.Symbols.N <- rep(NA_character_, nrow(indf))
indf[!grepl("ON-TARGET", indf$Gene.Symbol), ]$Gene.Symbols.N <- gnmsnos

write.csv(indf, file = "RBP_CG_Assay G-CUSTOM-176687_UniqueGN.csv")
