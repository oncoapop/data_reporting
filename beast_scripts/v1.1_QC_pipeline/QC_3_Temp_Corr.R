# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

#install.packages("directlabels")
library(lattice)
library(directlabels)
library(calibrate)
require(MASS)
library(ddCt)

#######################################################################
SA="siPanel_Plate2"
# if run directly uncomment the sample name
# Command line `Rscript QC_3_Temp_Corr.R --no-save --no-restore --args siPanel_Plate2`

# This takes the 4th argument (see str above) which is sample name
#args <- commandArgs(trailingOnly = TRUE)
#SA <- args[1]

# This is where the input files from QC_qPCRmeltTemp3.sh script are deposited
# This is where the output of this R-script will be
path = "/home/dyap/Projects/Takeda_T3/CG"
setwd(path)

#########################################################################


############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################


# Find the file names with the same name and pattern in the working dir
pat=paste(SA,"*.*_Temp.csv",sep="-")
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.csv, header = TRUE)

# Format
# Sample,Primer,RQ,actemp,actemp2,actemp3,caltemp,ampleng

# Read and plot the Temp Correlations for each file and out each as a PDF
for (i in 2:length(file_names)) {  
	print(i)
        # Filter out incomplete cases for flagging
        filt <- QC_list[[i]][!complete.cases(QC_list[[i]]$Actual_Temp),]
        # Include only complete cases with no NA values
        cases <- QC_list[[i]][complete.cases(QC_list[[i]]$Actual_Temp), ] 
#        cases <- (QC_list[[i]])
        # Filter out duplicated primer set for flagging
        dups <- as.character(QC_list[[i]][,2][duplicated(QC_list[[i]][,2])])

        # Gets the sample type T,N,Xn etc from the name
        sample <- unique(QC_list[[i]][1])
        name <- paste(paste("Temp_Corr",SA,sep="_"), sample[[1]], sep="-")
        fname <- paste(path,name,sep="/")

        # Output each as a separate PDF
        pdf (file=paste(fname, ".pdf", sep=""))

        # Formats the plots
        title <- paste(paste("Corr. of Amplicon Melting Temps for", SA, sep=" "), sample[[1]], sep="-")

	# This modules removes the outliers
	# reference : http://stackoverflow.com/questions/4666590/remove-outliers-from-correlation-coefficient-calculation
	res<-cases$actemp
	# Manual inspection to remove outliers
	want <- which(res >= 70)
	selected<-cases[want,]
	plot(selected$actemp, selected$caltemp, xlab="Actual Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main=title, pch=19, ylim=c(40,100), xlim=c(40,100), type="n")

	points(selected$actemp, selected$caltemp, col = "black", pch = 21, bg = "black", cex = 0.8)
        abline(lm(selected$caltemp ~ selected$actemp))

	donotwant <- which(res < 70)
	outliers<-cases[donotwant,]
	points(outliers$actemp, outliers$caltemp, col = "red", pch = 21, bg = "red", cex = 0.8)

        text(70,95, "pearson's, r =")
        text(105,95, signif(cor(selected$caltemp[want], selected$actemp[want], method = "pearson",use="pairwise.complete.obs")),3)

        # This section prints out the wells that failed (missing data or no data)
        createCounter <- function(value) { function(j) { value <<- value+j} }

	if (length(filt$Primer)>0) {
	        counter <- createCounter(50)
               text(80,60, "Primers with missing data / no RQ data:")
              for (j in 1:length(filt$Primer)) {x <- counter(round(50/length(filt$Primer))); text(x, 56, as.character(filt$Primer[j]),srt=45, pos=2,cex=0.6); print(x)}
					}

	if (length(outliers$Primer)>0) {
	        counter <- createCounter(70)
               text(50,100, "Primers filtered by Melt curve:", cex=0.8)
              for (j in 1:length(outliers$Primer)) {y <- counter(2); text(60, y, as.character(outliers$Primer[j]),srt=0, pos=2,cex=0.6); print(y)}
	        counter <- createCounter(70)
              for (j in 1:length(outliers$RQ)) {y <- counter(2); text(73, y, as.character(outliers$RQ[j]),srt=0, pos=2,cex=0.6); print(y)}
					}
#	textxy(outliers$actemp, outliers$caltemp, outliers$Primer, cx=0.8, m=c(mean(jitter(outliers$actemp)),mean(jitter(outliers$caltemp))))
	dev.off() 
# Output the netx step
#source("https://bioconductor.org/biocLite.R")
#biocLite("ddCt")

}


