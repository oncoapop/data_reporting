# This R-script uses the ddCt package to analyse the Ct value and 
# compute RQ values

# myR 
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# source("http://bioconductor.org/biocLite.R")
# biocLite("ddCt")

library(Biobase)
library(lattice)
library(RColorBrewer)
library(ddCt)
datadir <- "/home/dyap/Projects/Takeda_T3/CG_factors/qPCR_summary"
savedir <- datadir

# Find the file names with the same name and pattern in the working dir
setwd(datadir)
pat="siPanel_Plate*.*summary.csv"
file.names = list.files(pattern = pat)

info <- datadir("sampleData.txt")
warningFile <- savedir("warnings.txt")

CtData <- SDMFrame(file.names)

