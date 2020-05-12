# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

library(calibrate)

#######################################################################
SA="SA535"
# if run directly uncomment the sample name
# Command line `Rscript QC_4_Temp_Corr.R --no-save --no-restore --args DAH54`

# This takes the 4th argument (see str above) which is sample name
#args <- commandArgs(trailingOnly = TRUE)
#SA <- args[1]

# This is where the input files from QC_qPCRmeltTemp3.sh script are deposited
setwd("/home/dyap/dyap_temp/")
# This is where the output of this R-script will be
path = paste("/home/dyap/Projects/Tumour_Evol/Xenodrug",SA,sep="/")

#########################################################################


############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################

# Find the file names with the same name and pattern in the working dir
pat=paste(paste("Temp",SA,sep="_"),"*.*.txt", sep="-")
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.csv, header = TRUE)


# Read and plot the Temp Correlations for each file and out each as a PDF
for (i in 1:length(file_names)) {  

	# Removes null CT column nd removes null data from plate1 (not there)
		QC_list[[i]]$CT<-NULL
		QC_list[[i]]<-QC_list[[i]][!is.na(QC_list[[i]]$Actual_Temp),]
        # Filter out incomplete cases for flagging
        filt <- QC_list[[i]][!complete.cases(QC_list[[i]]),]
        # Include only complete cases with no NA values
        cases <- na.omit(QC_list[[i]])
        # Filter out duplicated primer set for flagging
        dups <- as.character(QC_list[[i]][,2][duplicated(QC_list[[i]][,2])])

        # Gets the sample type T,N,Xn etc from the name
        sample <- strsplit(colnames(QC_list[[i]]), split="_")[[2]][2]
        name <- paste(paste("Temp_Corr",SA,sep="_"), sample, sep="-")
        fname <- paste(path,name,sep="/")

        # Output each as a separate PDF
        pdf (file=paste(fname, ".pdf", sep=""))

        # Formats the plots
        title <- paste(paste("Correlation of Amplicon Melting Temperatures for", SA, sep=" "), sample, sep="-")
        plot(cases$Actual_Temp, cases$Cal_Temp, xlab="Actual Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main=title, pch=19, ylim=c(60,100), xlim=c(60,100))

        abline(lm(cases$Cal_Temp ~ cases$Actual_Temp))

        text(66,95, "r =")
        text(70,95, round(cor(cases$Cal_Temp, cases$Actual_Temp, method = "pearson"),3))

        text(66,90, "No. Amplicons =")
        text(72,90, nrow(cases))

        # This section prints out the wells that failed (missing data or no data)
        createCounter <- function(value) { function(j) { value <<- value+j} }

	if (length(filt$PCR_well) > 0) {
	        counter <- createCounter(60)
                text(70,68, "Wells with missing data / no data:")
                for (j in 1:length(filt$PCR_well)) {x <- counter(round(40/length(filt$PCR_well))); text(x, 65, filt$PCR_well[j], cex=0.6); print(x)}
					}

	if (length(dups) > 0)	{
	        counter <- createCounter(55)
                text(67,62, "Duplicate primer sets:")
                for (j in 1:length(dups)) {x <- counter(round(18/length(dups))); text(x, 60, dups[j], cex=0.6); print(x)}
				}

	# This labels every other point on the plot for QC purposes
	textxy(cases$Actual_Temp, cases$Cal_Temp, cases$PCR_well, cex=0.8, m=c(mean(cases$Actual_Temp),mean(cases$Cal_Temp)))

	dev.off() 

	# Output each as a separate PDF
	#	pdf (file=paste(fname, "_pairs.pdf", sep=""))
	#	pairs(cases)
	#	dev.off()
	print("end")
}

q()


