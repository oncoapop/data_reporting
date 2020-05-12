# Script to GO of siCG factors

# ON beast
library(plyr)
library(ggrepel)
library(ggplot2)

wd="/share/lustre/backup/dyap/Projects/Takeda_T3/CG_factors"
setwd(wd)

# Submit this file to http://www.pantherdb.org/geneListAnalysis.do
# remove the "'" from 3' as it causes import failure
#filein="siCFfactor_GO_BP.txt"
filein="CG_RNA_motifs_full_GO_BP.txt"
fileout=paste(filein,"processed", sep="_")

# cat CG_RNA_motifs_full_GO_BP.txt  | sed 's/'\''/\-prime/' > CG_RNA_motifs_full_GO_BP.txt_processed
GO<-read.table(file=fileout, header=TRUE, skip=10, sep="\t")

####################

GOin<-GO
colnames(GOin)[1]<-"PantherGO"
colnames(GOin)[6]<-"Fold_Enrichment"
colnames(GOin)[7]<-"Adjusted.P.value"
GO <- GOin[order(-GOin$Adjusted.P.value),]


# Set value >5 to 5
GO$Enrichment <- as.numeric(gsub(" > ","", GO$Fold_Enrichment))
GO$PantherGO<-gsub("\\s*\\([^\\)]+\\)","",as.character(GO$PantherGO))

#GO$Label <- do.call(paste, c(GO[c("si.factors.with.motif.den.change.in.CG..32.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))
GO$Label <- do.call(paste, c(GO[c("CG.enriched_factors..46.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))

# annotation of no / total in category

# Core wrapping function

wrap.it <- function(x, len)
{
  sapply(x, function(y) paste(strwrap(y, len),
                              collapse = "\n"),
         USE.NAMES = FALSE)
}

# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

GO$value<-round(as.numeric(-log10(GO$Adjusted.P.value)),1)

#subsetting the data
gosub <-subset(GO, Enrichment > 20 & Enrichment != "Inf")
#gosub <-subset(GO, value > 10 & Enrichment != "Inf")
gosub$lab<-wrap.labels(gosub$PantherGO,20)

gosub$Fold_Enrichment<-round(gosub$Enrichment,1)

q<-ggplot(data=gosub, aes(x=reorder(lab,value), y=value)) +
  geom_bar(width=0.8, stat="identity")+
  geom_text(data=gosub, aes(x=lab, y=value, hjust=-0.15, label=Label)) +
  geom_text(data=gosub, aes(x=lab, y=value, hjust=1.2, label=Fold_Enrichment), color="white") +
#  scale_y_continuous(limits=c(0,28)) +
#  ggtitle("CLK2 interactors by Gene Ontology") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  labs(y = "-log10(Adjusted p-value)", x="Top GO Biological Processes by Fold Enrichment") +
  coord_flip()



#pdf(file="SupplFig_siCGfactors_byGO_BP.pdf", useDingbats=FALSE)
pdf(file="SupplFig_CGfactors_full_byGO_BP.pdf", useDingbats=FALSE)
q
dev.off()


#library(gridExtra)
#grid.arrange(q, r, nrow=2)

############################

# Getting the genes involved in 
get=c("3-prime-UTR binding","3-prime-end processing")

#filein="siCGfactor_GOList.txt"
filein="CG_RNA_motifs_full_GOList.txt"

GO<-read.table(file=filein, header=FALSE, skip=0, sep="\t")

# grep "end processing" siCGfactor_GOList.txt | awk -F"\t" '{print $2}' | awk -F";" '{print $2}' > CG_3endproc_V7.txt
#system("grep \"3-prime-end processing\"siCGfactor_GOList.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CG_3endproc_V7.txt")
#system("grep \"3-prime-UTR binding\"siCGfactor_GOList.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CG_3UTRbind_V7.txt")
#system("grep \"poly(A) RNA binding\"siCGfactor_GOList.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CG_polyAbind_V7.txt")

####################

#GO<-mutate(GO, UTR=ifelse(grepl("3-prime-UTR binding", GO$V13), "3-prime-UTR binding (n=6/32)",""))
#GO<-mutate(GO, END=ifelse(grepl("3-prime-end processing", GO$V14), "3-prime-end Processing (n=7/32)",""))
#GO<-mutate(GO, polyA=ifelse(grepl("poly\\(A\\) RNA binding", GO$V13), "poly(A) RNA binding (n=29/32)",""))
GO<-mutate(GO, UTR=ifelse(grepl("3-prime-UTR binding", GO$V8), "3-prime-UTR binding (n=7/46)",""))
GO<-mutate(GO, END=ifelse(grepl("3-prime-end processing", GO$V7), "3-prime-end Processing (n=7/46)",""))
GO<-mutate(GO, polyA=ifelse(grepl("poly\\(A\\) RNA binding", GO$V8), "poly(A) RNA binding (n=39/46)",""))

table(GO$END)
table(GO$UTR)
table(GO$polyA)

GOcom<-GO
GOcom$END[GOcom$END == "3-prime-end Processing (n=7/46)"]  <- "YES"
GOcom$UTR[GOcom$UTR == "3-prime-UTR binding (n=7/46)"]  <- "YES"
GOcom$polyA[GOcom$polyA == "poly(A) RNA binding (n=39/46)"]  <- "YES"

names(GOcom)
colnames(GOcom)[9]
colnames(GOcom)[9]<-"Involved_in_3'UTR_binding"
colnames(GOcom)[10]
colnames(GOcom)[10]<-"Involved_in_3'-end_Processing"
colnames(GOcom)[11]
colnames(GOcom)[11]<-"Involved_in_poly(A)_RNA_binding"

#foo <- data.frame(do.call('rbind', strsplit(as.character(GOcom$V1),'|',fixed=TRUE)))
GO1<-within(GOcom, UniProtID<-gsub("UniProtKB=","",as.character(do.call('rbind', strsplit(as.character(GOcom$V1), '|', fixed=TRUE))[,3])))
GO2<-within(GO1, GeneID<-as.character(do.call('rbind', strsplit(as.character(GO1$V3), ';', fixed=TRUE))[,2]))
GO3<-within(GO2, Description<-as.character(do.call('rbind', strsplit(as.character(GO2$V3), ';', fixed=TRUE))[,1]))

colnames(GO3)
suptab <- GO3[c(12,13,14,9,10,11)]
#suptab <- GO3[c(20,19,21,17,16,18)]


sapply(suptab,class)

write.table(suptab, file="CG_RNA_motifs_full_3endfun.tsv", quote=FALSE, sep="\t")
#write.table(suptab, file="CGfactors_3endfun.tsv", quote=FALSE, sep="\t")

####################
#Old method

U<-GO[grepl("3-prime-UTR binding", GO$V13),]
E<-GO[grepl("3-prime-end processing", GO$V14),]
A<-GO[grepl("poly\\(A\\) RNA binding", GO$V13),]


UTR<-U[c(1,2)]
END<-E[c(1,2)]
polyA<-A[c(1,2)]

write.table(UTR, file = "3'UTR", append = FALSE, quote = FALSE, sep = ",")
write.table(END, file = "3'END", append = FALSE, quote = FALSE, sep = ",")
write.table(polyA, file = "polyA", append = FALSE, quote = FALSE, sep = ",")
