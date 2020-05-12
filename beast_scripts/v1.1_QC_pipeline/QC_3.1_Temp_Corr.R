# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

#install.packages("directlabels")
library(lattice)
library(directlabels)
library(calibrate)
require(MASS)
library(ddCt)

#######################################################################
# SA="siPanel_Plate2"
# SA="siPanel_Plate3"
SA="siPanel_Plate4"
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
QC_list = lapply(file_names, read.csv, header = TRUE, na.strings="Undetermined")

# Format
#   Sample.Name                       Target.Name      Ct.Mean    Tm1	Tm2  Tm3 caltemp ampleng

outpath="/home/dyap/Projects/Takeda_T3/CG_factors/Temp_Corr"
setwd(outpath)

# Read and plot the Temp Correlations for each file and out each as a PDF
for (i in 1:length(file_names)) {  
######################################################################  Manual
	print(i)
        # Filter out incomplete cases for flagging
        filt <- QC_list[[i]][!complete.cases(QC_list[[i]]$Ct.Mean),]
        # Include only complete cases with no NA values
        cases <- QC_list[[i]][complete.cases(QC_list[[i]]$Ct.Mean), ] 
#        cases <- (QC_list[[i]])
        # Filter out duplicated primer set for flagging
        dups <- as.character(QC_list[[i]][,2][duplicated(QC_list[[i]][,2])])

        # Gets the sample type T,N,Xn etc from the name
        sample <- unique(QC_list[[i]][1])
	if(is.na(sample)) {
            	next
		}

        name <- paste(paste("Temp_Corr",SA,sep="_"), sample[[1]], sep="-")
        fname <- paste(path,name,sep="/")

######################################################################  Manual
        # Output each as a separate PDF
        pdf (file=paste(fname, ".pdf", sep=""))
######################################################################  Manual

        # Formats the plots
        title <- paste(paste("Corr. of Amplicon Melting Temps for", SA, sep=" "), sample[[1]], sep="-")

	# This modules removes the outliers
	# reference : http://stackoverflow.com/questions/4666590/remove-outliers-from-correlation-coefficient-calculation
	res<-cases$Tm1
	# Manual inspection to remove outliers
	want <- which(res >= 70)
	selected<-cases[want,]
	plot(selected$Tm1, selected$caltemp, xlab="Actual Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main=title, pch=19, ylim=c(40,100), xlim=c(40,100), type="n")

	points(selected$Tm1, selected$caltemp, col = "black", pch = 21, bg = "black", cex = 0.8)
        abline(lm(selected$caltemp ~ selected$Tm1))

	donotwant <- which(res < 70)
	outliers<-cases[donotwant,]
	points(outliers$Tm1, outliers$caltemp, col = "red", pch = 21, bg = "red", cex = 0.8)

        text(70,95, "pearson's, r =")
        text(105,95, signif(cor(selected$caltemp[want], selected$Tm1[want], method = "pearson",use="pairwise.complete.obs")),3)

        # This section prints out the wells that failed (missing data or no data)
        createCounter <- function(value) { function(j) { value <<- value+j} }

	if (length(filt$Target.Name)>0) {
	        counter <- createCounter(50)
               text(80,60, "Primers with missing data / no CT data:")
              for (j in 1:length(filt$Target.Name)) {x <- counter(round(50/length(filt$Target.Name))); text(x, 56, as.character(filt$Target.Name[j]),srt=45, pos=2,cex=0.6); print(x)}
					}

	if (length(outliers$Target.Name)>0) {
	        counter <- createCounter(70)
               text(50,100, "Primers filtered by Melt curve:", cex=0.8)
              for (j in 1:length(outliers$Target.Name)) {y <- counter(2); text(60, y, as.character(outliers$Target.Name[j]),srt=0, pos=2,cex=0.6); print(y)}
	        counter <- createCounter(70)
              for (j in 1:length(outliers$Ct.Mean)) {y <- counter(2); text(73, y, as.character(outliers$Ct.Mean[j]),srt=0, pos=2,cex=0.6); print(y)}
					}
######################################################################  Manual
	dev.off() 
}

# Output the next step

#source("https://bioconductor.org/biocLite.R")
#biocLite("ddCt")




