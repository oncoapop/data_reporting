# This R-script uses the HTqPCR
# package to analyse the Ct value and compute RQ values

# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# source("http://bioconductor.org/biocLite.R")
# biocLite("HTqPCR")
library("HTqPCR")

path<-"/home/dyap/Projects/Takeda_T3/siRNA"

# Find the file names with the same name and pattern in the working dir
setwd(path)

# command line command
# ls -lr | grep 151104 | awk -F" " '{print $9" "$10}' > files.txt

# Read all the files name of the experiment in to a file (files.txt)
files<-read.table(file.path(path, "files.txt"), sep=",")

# Reads the file name into Q (in a list)
# Deleted the lines before and after the formatted text manually

# Read each individual file into qPCRset manually
# Read each file separately into data frame

filename=as.character(files$V1[2])

raw<-read.table(file=filename, header=FALSE, sep="\t")

# two time points (each =2)
Genotype = rep(c("siSRRM1_SRRM1","siSRRM2_SRRM1","siNT_SRRM1",
                 "siSRRM1_SRRM2","siSRRM2_SRRM2","siNT_SRRM2",
                 "siSRRM1_GAPDH","siSRRM2_GAPDH","siNT_GAPDH"), each=2)

# features (3 detectors)
#features=length(levels(raw$V6))
features=1
# Set of data in each time point
data=54

# flag = flagged column
# feature = no of targets
# type = control or target


raw1 <- readCtData(as.character(files$V1[2]), path = path, n.features = features, n.data = data, 
	column.info = list(flag = "Flag", feature = "Detector", type = "Task", Ct = "Ct", position = "Pos"), sep = "\t", 
	header=TRUE)
#sampleNames(raw1) <- gsub("Expt","",gsub("siRNA_", "", gsub(" ","_",gsub(".*]","",levels(raw$V5)))))
sampleNames(raw1) <- Genotype

# Structure the data into replicates, genotypes
pData(raw1) <- data.frame(Genotype = Genotype)

names<-sampleNames(raw1)
sampleNames(raw1)
featureNames(raw1)
head(fData(raw1))

# Overall Ct distirbution but treatment (QC)
plotCtOverview(raw1, genes = featureNames(raw1), xlim = c(0, 50), groups=sampleNames(raw1),
	conf.int = TRUE, ylim = c(0, 55))

#plotCtOverview(raw1, genes = featureNames(raw1), xlim = c(0, 50), groups = levels(QC_list[[1]]$Sample), 
#	calibrator = "GAPDH SYBR")


# position controls (by sample - assumes it is one sample one plate)
plotCtCard(raw1, col.range = c(1, 35), well.size = 2.6)


# duplicates - no duplicates in this case
plotCtReps(raw1, card = 2, percent = 20)

# Variation across samples
par(mar=c(10, 4, 2, 2) + 0.1)
raw.mix <- raw1
exprs(raw.mix)[, 6] <- sample(exprs(raw1[, 6]))

plotCtVariation(raw.mix, variation = "sd", log = TRUE,
	main = "SD of replicated features", col = "lightgrey", las=2, names=names,
	cex.axis=0.8)

# detail of variation
raw.variation <- plotCtVariation(raw.mix, type = "detail",
	add.featurenames = TRUE, pch = 20, cex = 0.8)

names(raw.variation)

raw.variation[["Var"]]
raw.variation[["Mean"]]
apply(raw.variation[["Var"]][, 2:15], 2, summary)
colSums(raw.variation[["Var"]][, 2:15] > 20)

# Feature categories and filtering of dataset

raw.cat <- setCategory(raw1, groups = names,
	 quantile = 0.8)

plotCtCategory(raw.cat)
#plotCtCategory(raw.cat, stratify = "class")
plotCtCategory(raw.cat, by.feature = TRUE, cexRow = 0.1)

# NORMALIZATION

# quantile Will make the distribution of Ct values more or less identical across samples
 q.norm <- normalizeCtData(raw.cat, norm = "quantile")

 sr.norm <- normalizeCtData(raw.cat, norm = "scale.rank") 
  nr.norm <- normalizeCtData(raw.cat, norm = "norm.rank")

# Normalize
 d.norm <- normalizeCtData(raw.cat, norm = "deltaCt",
  deltaCt.genes = c("GAPDH SYBR"))
g.norm <- normalizeCtData(raw.cat, norm = "geometric.mean")

plot(exprs(raw1), exprs(q.norm), pch = 20, main = "Quantile normalization",
	col = rep(brewer.pal(6, "Spectral"), each = 14))
plot(exprs(raw1), exprs(d.norm), pch = 20, main = "dCt normalization",
	col = rep(brewer.pal(6, "Spectral"), each = 14))


# Filtering and subsetting the data
# cannot perform them as rank test fails

# quality assessment

plotCtCor(raw1, main = "Ct correlation")

# Distribution of Ct values

summary(raw1)
plotCtDensity(q.norm)
plotCtHistogram(q.norm)

plotCtDensity(d.norm)
plotCtHistogram(d.norm)

# Scatter across all samples
plotCtPairs(q.norm, col = "type", diag = TRUE)
plotCtPairs(d.norm, col = "type", diag = TRUE)


# Ct heatmaps
plotCtHeatmap(raw1, gene.names = featureNames(raw1), dist = "euclidean")
plotCtHeatmap(q.norm, gene.names = featureNames(d.norm), dist = "euclidean")
plotCtHeatmap(d.norm, gene.names = featureNames(d.norm), dist = "euclidean")

# Coefficients of variation
plotCVBoxes(raw1, stratify = "type")

# Clustering

clusterCt(q.norm, type = "samples")

cluster.list <- clusterCt(q.norm, type = "genes",
	n.cluster = 6, cex = 0.5)

# Principal components analysis
plotCtPCA(raw1)
plotCtPCA(raw1, features = FALSE)

# Differential expression
getCtHistory(q.norm)
getCtHistory(d.norm)

# Two sample t-tests
 qDE.ttest <- ttestCtData(q.norm, groups = levels(QC_list[[1]]$Sample),
 calibrator = "Control")

# Multiple sample types - limma

design <- model.matrix(~0 + levels(QC_list[[1]]$Sample)[3:14])
colnames(design) <- gsub("Expt","",gsub("siRNA_", "", gsub(" ","_",gsub(".*]","",colnames(design) ))))
print(design)
colnames(design) 

contrasts <- makeContrasts(	SRRM1_24hr_1 - siNT_24hr_1, SRRM2_24hr_1 - siNT_24hr_1,
				SRRM1_24hr_2 - siNT_24hr_2, SRRM2_24hr_2 - siNT_24hr_2,
				SRRM1_48hr_1 - siNT_48hr_1, SRRM2_48hr_1 - siNT_48hr_1,
				SRRM1_48hr_2 - siNT_48hr_2, SRRM2_48hr_2 - siNT_48hr_2,
				levels = design)
#colnames(contrasts) <- c("SRRM1 RQ", "SRRM2 RQ", "Control")

# Removal of blank and RT inactivated samples ie Samples 1,2
q.norm2 <- q.norm[order(featureNames(q.norm)),3:14]
d.norm2 <- d.norm[order(featureNames(q.norm)),3:14]

sampleNames(q.norm2)
featureNames(q.norm2)
head(fData(q.norm2))

qDE.limma <- limmaCtData(q.norm2, design = design,
	contrasts = contrasts, ndups = 1, spacing = 1)

qDE.limma <- limmaCtData(d.norm2, design = design,
	contrasts = contrasts, ndups = 3, spacing = 1)

# ok can't normalize so try on raw data
#######################################
#  install.packages("statmod")
library("statmod")
design <- model.matrix(~0 + levels(QC_list[[1]]$Sample)[3:14])
colnames(design) <- gsub("Expt","",gsub("siRNA_", "", gsub(" ","_",gsub(".*]","",colnames(design) ))))
print(design)
colnames(design)

contrasts <- makeContrasts(     SRRM1_24hr_1 - siNT_24hr_1, SRRM2_24hr_1 - siNT_24hr_1,
                                SRRM1_24hr_2 - siNT_24hr_2, SRRM2_24hr_2 - siNT_24hr_2,
                                SRRM1_48hr_1 - siNT_48hr_1, SRRM2_48hr_1 - siNT_48hr_1,
                                SRRM1_48hr_2 - siNT_48hr_2, SRRM2_48hr_2 - siNT_48hr_2,
                                levels = design)
#colnames(contrasts) <- c("SRRM1 RQ", "SRRM2 RQ", "Control")

# Removal of blank and RT inactivated samples ie Samples 1,2
#raw2 <- raw1[order(featureNames(raw1)),3:14]
raw2 <- raw1[order(featureNames(raw1))]

sampleNames(raw2)<-names
featureNames(raw2)
head(fData(raw2))

qDE.limma <- limmaCtData(raw2, design = design,
        contrasts = contrasts, ndups = 1, spacing=0 )

