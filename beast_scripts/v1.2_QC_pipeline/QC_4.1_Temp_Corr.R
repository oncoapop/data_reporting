# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name between #### in this script to make it run correctly

library(calibrate)
library(Hmisc)
library(ggrepel)

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

 a<-rcorr(cases$Actual_Temp, cases$Cal_Temp, type="pearson")
 r<-signif(a$r[1,2], 3)
 P<-signif(a$P[1,2],3)
 attributes(a)

        file=paste(fname, ".pdf", sep="")

        title <- paste(paste("Correlation of Amplicon Melting Temperatures for", SA, sep=" "), sample, sep="-")
	ylab <- expression("Calculated Amplicon Tm " ( degree~C))
	xlab <- expression("Actual Amplicon Tm " ( degree~C))

		xleft<-round(range(cases$Actual_Temp)[1])-1
		xright<-round(range(cases$Actual_Temp)[2])+1
		ybottom<-round(range(cases$Cal_Temp)[1])-1
		ytop<-round(range(cases$Cal_Temp)[2])+1

p = ggplot(cases, aes(cases$Actual_Temp, cases$Cal_Temp), pch=20) +
            geom_point(data=cases, aes(colour = "black"))+
            geom_smooth(colour = "red", fill = "lightgreen", method = 'lm') +
            geom_text_repel(data=cases, aes(label=PCR_well)) +

            annotate("label", x =  xleft+5, y = ytop-1, label = paste("Pearsons~rho== ",r),parse = TRUE, fontface = "bold") +
            annotate("label", x =  xleft+5, y= ytop-3 , label = paste("No~of~Amplicon == ", nrow(cases)),parse = TRUE, fontface = "bold") +
	    annotate("label", x = xright-4, y = ybottom+round(1/8*(ytop-ybottom)), label = "Wells with missing data / no data:") + 
	    annotate("label", x = (round(xright)-length(filt$PCR_well)):(round(xright)-1), y = ybottom+round(1/8*(ytop-ybottom))-1, label = filt$PCR_well) +
	    annotate("label", x = xright-4, y = ybottom+round(1/4*(ytop-ybottom)), label = "Duplicated primer sets:") + 
	    annotate("label", x = (round(xright)-length(dups)):(round(xright)-1), y =  ybottom+round(1/4*(ytop-ybottom))-1, label = dups) +

	  theme(legend.position = "none") + ggtitle(title) +
  	theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  	theme(axis.text=element_text(size=12),
        	axis.title=element_text(size=14,face="bold")) +
	labs(x=xlab) +  
	labs(y=ylab)  

	p
	ggsave(file, width = 10, height = 10)
	print(fname)


}



