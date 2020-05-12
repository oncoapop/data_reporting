# Reference : http://genomicsclass.github.io/book/pages/GEOquery.html

#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")

library(GEOquery)

### This will download a 20 Mb
gse <- getGEO("GSE21653", GSEMatrix = TRUE)
show(gse)

filePaths = getGEOSuppFiles("GSE21653")
filePaths


dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])

df1 <- getGSEDataTables("GSE3494")
lapply(df1, head)


