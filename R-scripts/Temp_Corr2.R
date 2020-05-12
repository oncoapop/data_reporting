# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

library(calibrate)

#######################################################################
# SA="SA499"
# if run directly uncomment the sample name
# Command line `Rscript Temp_Corr2.R --no-save --no-restore --args SA499`

# This takes the 4th argument (see str above) which is sample name
args <- commandArgs(trailingOnly = TRUE)
SA <- args[4]

# This is where the input files from QC_qPCRmeltTemp2.sh script are deposited
setwd("/home/dyap/dyap_temp/")
# This is where the output of this R-script will be
path = "/share/lustre/backup/dyap/Projects/Tumour_Evol/Primer3_outputs"

#########################################################################


############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################

# Find the file names with the same name and pattern in the working dir
pat=paste(paste("Temp",SA,sep="_"),"*.*.txt", sep="-")
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.csv, header = TRUE)

# Read and plot the Temp Correlations for each file and out each as a PDF
for (i in 1:length(file_names)) {  
	
	# Filter out incomplete cases for flagging                          
	filt <- QC_list[[i]][!complete.cases(QC_list[[i]]),]
	# Include only complete cases with no NA values
	cases <- na.omit(QC_list[[i]])

	# Gets the sample type T,N,Xn etc from the name
	sample <- strsplit(colnames(QC_list[[i]]), split="_")[[2]][2]
	name <- paste(paste("Temp_Corr",SA,sep="_"), sample, sep="-")
	fname <- paste(path,name,sep="/")

	# Output each as a separate PDF
	pdf (file=paste(fname, ".pdf", sep=""))

	# Formats the plots								
	title <- paste(paste("Correlation of Amplicon Melting Temperatures for", SA, sep=" "), sample, sep="-")
	plot(cases$Actual_Temp, cases$Cal_Temp, xlab="Actual Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main=title, pch=19, ylim=c(60,100), xlim=c(60,100))

	# This section prints out the wells that failed (missing data or no data)
	createCounter <- function(value) { function(j) { value <<- value+j} }
	counter <- createCounter(60)
	
		text(70,68, "Wells with missing data / no data:")
		for (j in 1:length(filt$PCR_well)) {x <- counter(round(40/length(filt$PCR_well))); text(x, 65, filt$PCR_well[j], cex=0.6); print(x)}

	# This labels every other point on the plot for QC purposes
	textxy(cases$Actual_Temp, cases$Cal_Temp, cases$PCR_well, cx=0.8, m=c(mean(cases$Actual_Temp),mean(cases$Cal_Temp)))

	dev.off() 

	print("end")

}

q()


