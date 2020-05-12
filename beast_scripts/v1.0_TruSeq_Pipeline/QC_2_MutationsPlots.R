# R Script to display mutants in samples
# based on the script written by Derek Chu, Dec 2016
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rlibstree")
#require(Rlibstree)

library(ggplot2)
library(dplyr)

# Set working directory
dir="/home/dyap/Projects/ctDNA/AT094/variants"
setwd(dir)

# Use one category only
category="TUMOURonly"
# category="germline"
# category="WTonly"
pdffile=paste(category,"pdf",sep=".")

files.plot <- list.files(pattern = category)
plottitle=paste("QC from MSR VSFs: Distribution of Counts in", category, sep=" ")

# Combined data from vcf 
comb.dat <- sapply(files.plot, function(x) {
  tryCatch({
    dat <- read.csv(x, header = FALSE, sep=",", stringsAsFactors = FALSE)
    names(dat) <- c("CHAR", "POS", "IND", "REF", "ALT")
    return(dat)
  }, error = function(e)
    data.frame(matrix(NA_character_, ncol = 5,
                      dimnames = list(NULL, c("CHAR", "POS", "IND", "REF", "ALT")))))
}, simplify = FALSE) %>% 
  Map(cbind, ., Sample = names(.)) %>%
  plyr::rbind.fill()

if (category == "TUMOURonly") {
				comb.dat$ALT[1]<-"T"
				comb.dat$ALT[23]<-"T"
				}
comb.dat$ID<-paste(comb.dat$CHAR, comb.dat$POS, sep="_")
comb.dat$MUT<-paste(paste(comb.dat$REF, comb.dat$ALT, sep=">"),"\n",sep="")

nearfinal<-comb.dat[complete.cases(comb.dat),]

# module to get the chr # only from the chrn 
cc       <- strsplit(as.vector(nearfinal$CHAR),'r')
ch    <- unlist(cc)[2*(1:length(nearfinal$CHAR))-1]
nearfinal$CHRN    <- unlist(cc)[2*(1:length(nearfinal$CHAR))  ]

# Sorting using chromosome and pos 
final <- nearfinal[order(as.numeric(nearfinal$CHRN),as.numeric(nearfinal$POS)),]

# Get count data for every combination of sample * CHAR
count.dat <- with(final, expand.grid(Sample = unique(Sample), 
                                          ID = unique(ID) )) %>% 
  merge(., final, all.x = TRUE) %>% 
  plyr::ddply(., ~Sample + ID, function(x) {
    if (any(is.na(x$POS)))
      return(0)
    else
      return(nrow(x))
  }) %>% 
  dplyr::rename(Count = V1)

# Checking from original vcfs
len<-nrow(count.dat)
chk <- data.frame(ID = rep("", len),
                     Sample = rep("sample", len),
                     Chr = rep("chr", len),
                     Pos = rep(0, len),
                     Ref = rep("ATCG", len),
                     Alt = rep("ATCG", len),
                     Vcf = rep("vcfs", len),
                     Count = rep(0, len),
                     stringsAsFactors = FALSE)


        print("Checking positions against original vcfs")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

for (no in (1:len))

        {
		chk$ID[no]		<-	as.character(count.dat$ID)[no]
		chk$Sample[no]		<-	as.character(count.dat$Sample)[no]
		chk$Count[no]		<-	as.character(count.dat$Count)[no]
		chk$Ref[no]		<-	as.character(comb.dat[which(comb.dat$ID == chk$ID[no]),4 ])
		chk$Alt[no]		<-	as.character(comb.dat[which(comb.dat$ID == chk$ID[no]),5 ])
		chk$Chr[no]		<-	strsplit(as.character(count.dat$ID[no]),'_')[[1]][1]
		chk$Pos[no]		<-	strsplit(as.character(count.dat$ID[no]),'_')[[1]][2]

			# Checking module
			x	<-	chk$Sample[no]
			file	<-	substrRight(x, 2)
			grep  	<-	paste(paste("grep",chk$Pos[no],sep=" "),"*.vcf", sep=" ")
        	        grep2	<-	paste("grep",paste(file,"_",sep=""),sep=" ")
			awk 	<- 	paste(paste("awk '$4 == ", chk$Ref[no], sep='"'),"'",sep='"')
			awk2 	<- 	paste(paste("awk '$5 == ", chk$Alt[no], sep='"'),"'",sep='"')
			command	<-	paste(paste(paste(grep,grep2,sep=" | "),awk,sep=" | "),awk2,sep=" | ")

			print(command)

			match	<-	system(command, intern = TRUE)
			#match <- system("grep 24145675 *.vcf | grep 47_", intern = TRUE)

			dd    <- strsplit(match,'.vcf')
			dd2    <- unlist(dd)[2*(1:length(match))-1]
		
			print(dd2)

		chk$Vcf[no]<-paste(unique(dd2), collapse="\n")
		count.dat$Vcf[no]<-paste(unique(dd2), collapse="\n")

			match<-NULL
			dd<-NULL
			dd2<-NULL
	}


# Heatmap
pdf(file=pdffile, width=15, height=10)
ggplot(count.dat, aes(x = ID, y = Sample, fill = Count)) +
  geom_tile(color = "red") +
  scale_fill_continuous(breaks = sort(unique(count.dat$Count))) +
  geom_text(aes(label = count.dat$Vcf), cex=2.5) +
  labs(title = plottitle) +
  scale_fill_gradient(low="white", high="red")+
  theme(axis.text.x = element_text(angle = 90))
dev.off()


