library(Biobase)
library(lattice)
library(RColorBrewer)
library(ddCt)

# Worked example
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/ddCt/inst/doc/rtPCR.pdf

datadir <- function(x) system.file("extdata", x, package="ddCt")
savedir <- function(x) file.path(tempdir(), x)
file.names <- c(datadir("Experiment1.txt"),datadir("Experiment2.txt"))
info <- datadir("sampleData.txt")
warningFile <- savedir("warnings.txt")

# Reference sample and housekeeping gene
name.reference.sample <- c("Sample1", "Sample2")
name.reference.gene <- c("Gene2")

# Read in data
CtData <- SDMFrame(file.names)
sampleInformation <- read.AnnotatedDataFrame(info,header=TRUE, row.names=NULL)

# Apply the ddCt method
result <- ddCtExpression(CtData,
	calibrationSample=name.reference.sample,
	housekeepingGene=name.reference.gene,
	sampleInformation=sampleInformation,
	warningStream=warningFile)

# Check values
Ct(result)
CtErr(result)
dCt(result)
dCtErr(result)
ddCt(result)
ddCtErr(result)

# RQ values
exprs(result)
levelErr(result)

# QC
numberNA(result)


# Visualization
br <- errBarchart(result)
print(br)

# Write result to text file
elistWrite(result,file=savedir("allValues.txt"))

##################################

# My samples input 

datadir <- "/home/dyap/Projects/Takeda_T3/siRNA"
savedir <- function(x) file.path(tempdir(), x)
setwd(datadir)
file.names <- c("151104_SRRM1_2_siKD_qPCR.sdm")
info <- datadir("sampleData.txt")
warningFile <- savedir("warnings.txt")


# Reference sample(s) and housekeeping gene
name.reference.sample <- c("siNT 24hr Expt1","siNT 48hr Expt1","SRRM2 siRNA 48hr Expt1")
name.reference.gene <- c("GAPDH SYBR")

# Read in data
CtData <- SDMFrame(file.names)
#sampleInformation <- read.AnnotatedDataFrame(info,header=TRUE, row.names=NULL)

# Apply the ddCt method
data <- ddCtExpression(CtData,
        calibrationSample=name.reference.sample,
        housekeepingGene=name.reference.gene,
#        sampleInformation=sampleInformation,
        warningStream=warningFile)

# Check values
Ct(data)
CtErr(data)
dCt(data)
dCtErr(data)
ddCt(data)
ddCtErr(data)

# RQ values
exprs(data)
levelErr(data)

# QC


# Visualization
br <- errBarchart(data)
print(br)

# Write result to text file
elistWrite(result,file=savedir("allValues.txt"))
