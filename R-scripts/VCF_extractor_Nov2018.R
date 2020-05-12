source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
library("VariantAnnotation")
library("IRanges")
library("GenomicRanges")
library(foreign)
#####

# To run this script change the setwd()
setwd("/home/dyap/Projects/EIF4A3_paper/followup")
outfile="ouput.csv"
pdffile="pdffile.pdf"

# Change the number of samples

# check with "list(vcf_list)"
samples <- 3

# read from the pattern
#file_names = list.files(pattern = 'SA49*.*vcf')
file_names = system("ls /share/lustre/projects/takeda_EIF4A3/pipeline_inbox/TAK-82/TAK-162/SING_MUTA/SA79*/results/TASK_9_CBIO_FILTER_SA*_MutationSeq.annotSnpEff.annotMA.flagDBsnp.flag1000gen.flagCosmic.filtered.vcf",intern=TRUE)

# Extract all the VCFs into a concatenated VCF list
#vcf_list = lapply(file_names, readVcf, "hg19", sep = "\t")
vcf_list = lapply(file_names, readVcf, "GRCh37", sep = "\t")

getwd()

createCounter <- function(value) { function(i) { value <<- value+i} }
count <- createCounter(1)

#########################################################
# THIS IS THE CURRENT WORKING MODULE TO EXPERIMENT WITH #
#########################################################

sumdf <- data.frame(	Sample_ID = rep("", samples),
			Variants = rep(0, samples),
			Rows = rep("", samples),
			stringsAsFactors = FALSE)

for (rj in seq(samples)) 	{
			sid <- rownames(vcf_list[[rj]])
			varn <- nrow(vcf_list[[rj]])

sumdf$Samples_ID[rj] <- sid
sumdf$Variants[rj] <- varn

# For each of the list of samples
len <- nrow(vcf_list[[rj]])

	if ( len > 0 ){

# What do we want to get out of it
d.frame <- data.frame(	     ID = rep("", len),
			     Varalfreq = rep (0, len),
#			     Seqdepth = rep(0, len),
			     stringsAsFactors = FALSE)

for (ri in seq(len) ) 	{

		d.frame$ID[ri] <- rownames(vcf_list[[rj]][ri])
#		d.frame$Calls[ri] <- info(vcf_list[[rj]][ri])$GT
		d.frame$Varalfreq[ri] <- info(vcf_list[[rj]][ri])$VF
#		d.frame$Seqdepth[ri] <- geno(vcf_list[[rj]][ri])$DP

			}

# unique freq val or seq depth for summary plot
# must be the same length as names(d.frame)
names(d.frame)[2] <- paste(sid)

			} else { d.frame <- "NULL" };

assign(paste("Nuclei", rj, sep=""), d.frame)

# for first value 
	if ( rj == 1 && d.frame != "NULL"  ) {
		first <- d.frame }
 
# combining successive data.frames
	if ( rj == 2 && d.frame != "NULL"  ) { 
		sum1 <- merge(first, d.frame, by="ID", all=TRUE) }

# combining successive data.frames
	if ( rj > 2 && d.frame != "NULL"  ) { 
		a <- count(1)
		sum1 <- merge(sum1, d.frame, by="ID", all=TRUE) 
		} else { print("skip") }

}

write.table(sum1,file=outfile,sep=",",row.names=FALSE,col.names=TRUE)
# write.table(obj_name,file="C:\\Users\\dyap_000\\Documents\\R\\SA494\\SA494.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)

###################################

# Drawing heatmap

# Must be convert into a data.matrix (non-numeric converted to N/A)
ef <- data.matrix(sum1[2:ncol(sum1)])

# col headers - unique nuclei
names(sum1)

# Label rownames with ID
rownames(ef) <- sum1$ID

# heatmap(ef, Rowv=NA, Colv=NA, col = heat.colors(1024), scale="column", margins=c(5,10))

pdf(pdffile, width=6, height=6)
heatmap(ef, Rowv=NA, Colv=NA, na.rm=TRUE, main="Variant Freq-SA494 main run (MiSeq Reporter)", xlab="SA494 Samples", ylab="Position", cexCol=0.8, col=rev(heat.colors(1000)))
dev.off()
########################################################




# Next step - combing them 
# Ref: http://ryouready.wordpress.com/2009/01/23/r-combining-vectors-or-data-frames-of-unequal-length-into-one-data-frame/

#myList <- list(paste("Nuclei",seq(1:samples),sep=""))
myList <- list(Nuclei1, Nuclei2, Nuclei3)

l <- myList
do.call(rbind, lapply(lapply(l, unlist), "[",
        unique(unlist(c(sapply(l,names))))))

# merge two data frames by ID (includes all rows, introduced N/A is absent)
sum <- merge(Nuclei1,Nuclei2,by="ID", all=TRUE)

#######################################################



#######################################################
# Not useful if one file has 0 length
dataFiles <- lapply(Sys.glob("SA494*.vcf"), read.table)


# read txt files into a list (assuming separator is a tab)
data_list = lapply(file_names, read.csv, sep = "\t")
# not very useful -> all in one column!

# fl <- system.file("extdata", dataFiles[1], package="VariantAnnotation")
# cannot get this to work


#####################

for (rj in seq(samples)) 
 	{
obj_name <- data.frame	(ID 	= rep("", samples),
			length 	= rep("", len),
		     Varalfreq  = rep(0 , len),
			stringsAsFactors = FALSE)
			
len <- nrow(vcf_list[[rj]])


#####################
works for indiv files\

vcf <- readVcf("SA494-gDNA-Control-001_S52.vcf", "hg19")
# This works!

# Check out
info(vcf)
geno(vcf)

rowData(vcf_list[[1]])$REF[1]

# some examples
geno <- geno(vcf)$GT
varalfreq <- geno(vcf)$VF
seqdepth <- geno(info)$DP
aldepth <- geno(vcf)$DP
geno(vcf)$AD[1]
# There are two allelic freq
