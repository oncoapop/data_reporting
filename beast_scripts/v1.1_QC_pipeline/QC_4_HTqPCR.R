# This R-script uses the HTqPCR 
# package to analyse the Ct value and compute RQ values

# myR 
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# source("http://bioconductor.org/biocLite.R")
# biocLite("HTqPCR")
library("HTqPCR")

path<-"/home/dyap/Projects/Takeda_T3/CG_factors/qPCR_summary"

# Find the file names with the same name and pattern in the working dir
setwd(path)

#ls -lr | awk -F" " '{print $9}' > files.txt
# For incomplete plates (<384 wells) uncomment the following command

# Read all the files name of the experiment in to a file (files.txt)
files<-read.delim(file.path(path, "files.txt"))

# Reads the file name into Q (in a list)
QC_list = lapply(as.character(files$File), read.csv, header = TRUE)

# Read each individual file into qPCRset manually
# Read each file separately into data frame
# CG Panel has 55 primer sets 
primset=55

# No of Targets in run#1
# ABI format (old machine HT7900 format)
# reformatted, sample names always in col=2
sam1=6
primset1=51
raw1 <- readCtData(as.character(files$File[1]), path = path, n.features = primset1, n.data = sam1, column.info = list(flag = 6, feature = 3, type = 4, Ct = 5, position = 1), sep = ",", header=TRUE)
sampleNames(raw1)
featureNames(raw1)

# No of Targets in run#2
# New LifeTech format (new machine)
sam2=18
raw2 <- readCtData(as.character(files$File[2]), path = path, n.features = primset, n.data = sam2, column.info = list(feature = 3,Ct = 4, position = 1), sep = ",", header=TRUE)
sampleNames(raw2)
featureNames(raw2)

# No of Targets in run#3
sam3=18
raw3 <- readCtData(as.character(files$File[3]), path = path, n.features = primset, n.data = sam3, column.info = list(feature = 3, Ct = 4, position = 1), sep = ",", header=TRUE)
sampleNames(raw3)
featureNames(raw3)

# No of Targets in run#4
sam4=6
raw4 <- readCtData(as.character(files$File[4]), path = path, n.features = primset, n.data = sam4, column.info = list(feature = 3, Ct = 4, position = 1), sep = ",", header=TRUE)
sampleNames(raw4)
featureNames(raw4)

############################

# Check the samples!
# phenoData(raw1)
# pData(raw1)
# featureData(raw1)
# head(fData(raw1))

# Identify replicates (manually done!)
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(pData(raw1))),split='-',fixed=TRUE)))
pData(raw1)[2]<-as.character(foo$X1)
pData(raw1)[3]<-rep(1:2, 3)
colnames(pData(raw1))[2]<-"Genotype"
colnames(pData(raw1))[3]<-"Replicate"

#sampleNames(raw1)<-pData(raw1)$Genotype
# for raw 1 only
g <- featureNames(raw1)
plotCtOverview(raw1, genes = g, groups = sampleNames(raw1), conf.int = TRUE, ylim = c(0, 55))
plotCtOverview(raw1, genes = g, groups = sampleNames(raw1), calibrator = "siNT")

## XX
#raw1.mix <- raw1
#exprs(raw1.mix)[,5] <- sample(exprs(raw1[,5]))
#plotCtVariation(raw1.mix, variation = "sd", log = TRUE, main = "SD of replicated features", ylim=c(0,10000), col = "lightgrey")


# Feature categories and filtering
raw1.cat <- setCategory(raw1, groups = sampleNames(raw1), quantile = 0.8)
plotCtCategory(raw1.cat)
plotCtCategory(raw1.cat, by.feature = TRUE, cexRow = 0.6, cexCol=0.8)

# Normalisation
 q.norm <- normalizeCtData(raw1.cat, norm = "quantile")
# too few replicates for these normalization
# sr.norm <- normalizeCtData(raw1.cat, norm = "scale.rank")
# nr.norm <- normalizeCtData(raw1.cat, norm = "norm.rank")

g.norm <- normalizeCtData(raw1.cat, norm = "geometric.mean")
plot(exprs(raw1), exprs(q.norm), pch = 20, main = "Quantile normalisation", col = rep(brewer.pal(6, "Spectral"), each = 384))

# Filtering and subsetting the data
qFilt <- filterCtData(q.norm, remove.category = "Undetermined", remove.type = "Endogenous Control")

# Quality assessment
plotCtCor(raw1, main = "Ct correlation")
plotCtCor(qFilt, main = "Ct correlation")
summary(raw1)
plotCtHeatmap(raw1, gene.names = featureNames(raw1), dist = "euclidean")
plotCtHeatmap(qFilt, dist = "euclidean")
plotCVBoxes(raw1, stratify = "type")
plotCVBoxes(qFilt, stratify = "type")


# Clustering
clusterCt(q.norm, type = "samples")
cluster.list <- clusterCt(q.norm, type = "genes", n.cluster = 6, cex = 0.5)
c6 <- cluster.list[[6]]
print(c6)

plotCtPCA(raw1)
plotCtPCA(raw1, features = FALSE)
plotCtPCA(qFilt, features = FALSE)



#############
# 10 Differential expression

getCtHistory(q.norm)
getCtHistory(qFilt)

#siCPEB3
qDE.mwtest.siCPEB3 <- mannwhitneyCtData(q.norm[,c(1,2,6)], groups = as.factor(pData(raw1)$Genotype[c(1,2,6)]), calibrator = "siNT")
# siELAVL1
qDE.mwtest.siELAVL1 <- mannwhitneyCtData(q.norm[,c(3,4,6)], groups = as.factor(pData(raw1)$Genotype[c(3,4,6)]), calibrator = "siNT")
# siELAVL3
qDE.mwtest.siELAVL3 <- mannwhitneyCtData(q.norm[,c(5,6)], groups = as.factor(pData(raw1)$Genotype[c(5,6)]), calibrator = "siNT")

plotCtRQ(qDE.mwtest.siCPEB3, genes = 1:43, main="siCPEB3",  p.val = 0.05)
plotCtRQ(qDE.mwtest.siELAVL1, genes = 1:43, cexCol=0.8)
plotCtRQ(qDE.mwtest.siELAVL3, genes = 1:43, cexCol=0.8)



# multiple samples


design <- model.matrix(~0 + pData(raw1)$Genotype)
colnames(design) <- unique(pData(raw1)$Genotype)
contrasts <- makeContrasts(siCPEB3 - siNT, siELAVL1 - siNT, siELAVL3 - siNT, (siELAVL3 + siELAVL1 + siCPEB3)/3 - siNT, levels = design)

print(contrasts)
q.norm2 <- q.norm[order(featureNames(q.norm)), ]

qDE.limma <- limmaCtData(q.norm2, design = design, contrasts = contrasts, ndups = 2, spacing = 1)



