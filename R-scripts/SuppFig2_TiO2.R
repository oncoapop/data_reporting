# Script to analyse MS data
# R
# sudo R CMD INSTALL (pacakges downloaded from CRAN)
# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree" 
# Pleione CENTOS 6.7
# library(dplyr) 
# install.packages("seqRFLP")
# install.packages("png")

#

# ON beast
library(png)
library(seqRFLP)
library(plyr)
library(ggrepel)
library(ggplot2)
library(gridExtra)

wd="/home/dyap/Projects/Takeda_T3/TiO2"
pat="PeptideSeq_TiO2Expt.csv"
setwd(wd)
file_names = list.files(pattern = pat)

phos<-as.data.frame(read.csv(file_names[1],header=TRUE))

colnames(phos)
# reorder and subset data using [ ]
phospep <- phos[c("GeneName","Sequence","BHadj_OmnibusPval","LinTrendSlope","Effect_t3_500nM","Effect_t3_5000nM","AGDSMName")]


# Highlight known targets of CLK2
clk2<-as.data.frame(read.table(file="Predicted_CLK2_substrates.tsv", header=FALSE, skip=1, sep="\t"))

# Nice labelling (non-overlapping)
#phospep = mutate(phospep, Legend=ifelse(phospep$BHadj_OmnibusPval<0.005, "p<0.005 (n=505)", "Not Significant (n=3152)"))
#phospep = mutate(phospep, Legend=ifelse(phospep$BHadj_OmnibusPval<0.0005, "p<0.0005 (n=293)", "Not Significant (n=3364)"))
#phospep = mutate(phospep, Legend=ifelse(phospep$BHadj_OmnibusPval<0.05, "p<0.05 (n=906)", "Not Significant (n=2751)"))

phospep = mutate(phospep, Legend=ifelse(phospep$BHadj_OmnibusPval<0.05 & LinTrendSlope>1 , "p<0.05 and Trend>1 (n=204)", 
		ifelse(phospep$BHadj_OmnibusPval<0.05 & LinTrendSlope<(-1) , "p<0.05 and Trend<-1 (n=67)", "Not selected (n=3386)" ) ))

sig<-phospep[grep("p<0.05",phospep$Legend),]


phospep = mutate(phospep, Targets=ifelse(phospep$GeneName %in% clk2$V1, "Predicted targets by gene name", "Unknown"))

count(phospep$Legend)
count(phospep$Targets)
#labels
#########
p = ggplot(phospep, aes(LinTrendSlope, -log10(BHadj_OmnibusPval)), pch=20) +
  geom_point(aes(col=Legend)) +
#  theme_light()+geom_text_repel(data = phospep[grep("p<0.05",phospep$Legend),], aes(label=GeneName), fontface = "bold") +
  theme_light()+geom_text_repel(data = sig[grep("SRSF",sig$GeneName),], aes(label=GeneName), fontface = "bold") +
#  theme_light()+geom_text_repel(data = sig[grep("RM",sig$GeneName),], aes(label=GeneName), fontface = "bold", force = 7, 
#  nudge_x = -1.7, nudge_y=0.23) +
#  geom_point(aes(col="Targets")) + 
  scale_color_manual(values=c("grey90", "grey50", "grey50","transparent")) +
  theme(legend.position = c(0.15, 0.85)) + ggtitle("CLK2 targets by TiO2 enrichment MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  labs(x = "linear trend magnitude: DMSO control to T3 500nM to T3 5000nM", y="-log10(Adjusted p-value)")+
  geom_vline(xintercept=c(-1,1), na.rm = FALSE, show.legend = TRUE, linetype = "dotted") 

q = ggplot(phospep, aes(LinTrendSlope, -log10(BHadj_OmnibusPval)), pch=20) +
  geom_point(aes(col=Legend)) +
#  theme_light()+geom_text_repel(data = phospep[grep("p<0.05",phospep$Legend),], aes(label=GeneName), fontface = "bold") +
  theme_light()+geom_text_repel(data = sig[grep("SRRM1",sig$GeneName),], aes(label=GeneName), fontface = "bold") +
#  theme_light()+geom_text_repel(data = sig[grep("RM",sig$GeneName),], aes(label=GeneName), fontface = "bold", force = 7,
#  nudge_x = -1.7, nudge_y=0.23) +
#  geom_point(aes(col="Targets")) +
  scale_color_manual(values=c("grey90", "grey50", "grey50","transparent")) +
  theme(legend.position = c(0.15, 0.85)) + ggtitle("CLK2 targets by TiO2 enrichment MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  labs(x = "linear trend magnitude: DMSO control to T3 500nM to T3 5000nM", y="-log10(Adjusted p-value)")+
  geom_vline(xintercept=c(-1,1), na.rm = FALSE, show.legend = TRUE, linetype = "dotted")

r = ggplot(phospep, aes(LinTrendSlope, -log10(BHadj_OmnibusPval)), pch=20) +
  geom_point(aes(col=Legend)) +
#  theme_light()+geom_text_repel(data = phospep[grep("p<0.05",phospep$Legend),], aes(label=GeneName), fontface = "bold") +
  theme_light()+geom_text_repel(data = sig[grep("SRRM2",sig$GeneName),], aes(label=GeneName), fontface = "bold") +
#  theme_light()+geom_text_repel(data = sig[grep("RM",sig$GeneName),], aes(label=GeneName), fontface = "bold", force = 7,
#  nudge_x = -1.7, nudge_y=0.23) +
#  geom_point(aes(col="Targets")) +
  scale_color_manual(values=c("grey90", "grey50", "grey50","transparent")) +
  theme(legend.position = c(0.15, 0.85)) + ggtitle("CLK2 targets by TiO2 enrichment MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
#  labs(x = "linear trend magnitude: DMSO control to T3 500nM to T3 5000nM", y="-log10(Adjusted p-value)")+
  geom_vline(xintercept=c(-1,1), na.rm = FALSE, show.legend = TRUE, linetype = "dotted")


pdf(file="X1_T3_phosphopeptide.pdf")
#grid.arrange(p, q, r, ncol=3)
p
dev.off()



#######
phospep = mutate(phospep, Legend=ifelse(phospep$BHadj_OmnibusPval<0.05 & Effect_t3_500nM>1 , "p<0.05 and Trend>1 (n=204)", 
		ifelse(phospep$BHadj_OmnibusPval<0.05 & Effect_t3_500nM<(-1) , "p<0.05 and Trend<-1 (n=67)", "Not selected (n=3386)" ) ))
sig<-phospep[grep("p<0.05",phospep$Legend),]

q = ggplot(phospep, aes(Effect_t3_500nM, -log10(BHadj_OmnibusPval)), pch=20) +
  geom_point(aes(col=Legend)) +
#  theme_light()+geom_text_repel(data = phospep[phospep$GeneName %in% clk2$V1,], aes(label=GeneName), fontface = "bold") +
#  theme_light()+geom_text_repel(data = sig[grep("SRSF",sig$GeneName),], aes(label=GeneName), fontface = "bold") +
  theme_light()+geom_text_repel(data = sig[grep("SRRM1",sig$GeneName),], aes(label=GeneName), fontface = "bold") +
#  geom_point(aes(col=Targets)) + 
#  force = 7, nudge_x = -1.7, nudge_y=0.23) +
  scale_color_manual(values=c("grey90", "grey50", "grey50","red","transparent")) +
  theme(legend.position = c(0.15, 0.85)) + ggtitle("CLK2 targets by TiO2 enrichment MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  labs(x = "log2(fold change): T3 500nM relative to DMSO control", y="-log10(Adjusted p-value)")+
  geom_vline(xintercept=c(-1,1), na.rm = FALSE, show.legend = TRUE, linetype = "dotted") 
q









##################################################
targetpeptides<-phospep[grep("p<0.05",phospep$Legend),]

write.csv(targetpeptides, file="CLK2_targets.csv", quote=FALSE)

targetpeptides$fastaID<-gsub(" ","",make.unique(as.character(targetpeptides$GeneName)))
depl<-targetpeptides[targetpeptides$LinTrendSlope<(-1), c("fastaID","Sequence")]
enrich<-targetpeptides[targetpeptides$LinTrendSlope>1, c("fastaID","Sequence")]

mememinus<-depl[nchar(as.character(depl$Sequence)) > 8,]
memeplus<-enrich[nchar(as.character(enrich$Sequence)) > 8,]

dataframe2fas(mememinus, file = "CLK2_target_Seq_depleted.fa")
dataframe2fas(memeplus, file = "CLK2_target_Seq_enriched.fa")

#clk2<-as.list(as.character(targetpeptides$Sequence))
#lapply(clk2, write, "CLK2_target_sequence.txt", append=TRUE, ncolumns=1000)


#################
# Submit this file to http://www.pantherdb.org/geneListAnalysis.do
# remove the "'" from 3' as it causes import failure
GO<-read.table(file="CLK2-bind_MS-GO.txt", header=TRUE, skip=10, sep="\t")

colnames(GO)[1]<-"PantherGO"
colnames(GO)[6]<-"Fold_Enrichment"
colnames(GO)[7]<-"Adjusted.P.value"

# Set value >5 to 5
GO$Enrichment <- as.numeric(gsub(" > ","", GO$Fold_Enrichment))

q = ggplot(GO, aes(Enrichment, -log10(Adjusted.P.value)), pch=20) 

# Highlight known Uniprot interactors of CLK2
Int<-read.table(file="Uniprot_CLK2_interactors.txt", header=FALSE, skip=1, sep="\t")

# Nice labelling (non-overlapping)
MS = mutate(MS, Legend=ifelse(MS$All_p.adj<0.05 & MS$all.mean>1, "p<0.05 (n=180)",
        "Not Significant (n=1084)"))
MS = mutate(MS, Interactors=ifelse(MS$Gene %in% Int$V1, "Uniprot Interactors (11/26)", "Unknown"))

hl<-c("CLK1","CLK2","CLK3","CLK4")
MS = mutate(MS, CLKs=ifelse(MS$Gene %in% hl, "CLKs", "Unknown"))

table(MS$Interactors)
table(MS$CLKs)
known<-MS[MS$Gene %in% Int$V1,]

# Label Uniprot binary interactions

p = ggplot(MS, aes(all.mean, -log10(All_p.adj)), pch=20) +
  geom_point(aes(col=Legend)) +
  scale_color_manual(values=c("light grey", "dark grey")) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% Int$V1,], aes(label=Gene), fontface = "bold", 
  force = 7, nudge_x = -1.7, nudge_y=0.23) +
  geom_point(aes(col=Interactors)) + 
  scale_color_manual(values=c("grey90", "grey50", "red","transparent")) +
  theme(legend.position = c(0.15, 0.85)) + ggtitle("CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  labs(x = "log2(Mean enrichment ratios)", y="-log10(Adjusted p-value)") 

pdf(file="SupplFig_CLK2_Interactors.pdf")
p
dev.off()


#####################################
# label all those with GO splicing
Splice<-read.table(file="CLK2_int_Splicing.txt", header=FALSE)
# Highlight known Uniprot interactors of CLK2
Int<-read.table(file="Uniprot_CLK2_interactors.txt", header=FALSE, skip=1, sep="\t")
RBin<-read.table(file="CLK2_int_mRNAbinding.txt", header=FALSE)
allbind<-as.character(allbinding$Gene)

MS = mutate(MS, Legend=ifelse(MS$All_p.adj<0.05 & MS$all.mean>1, "p<0.05 (n=180)",
        "Not Significant (n=1084)"))
MS = mutate(MS, Function=ifelse(MS$Gene %in% Splice$V1, "RNA splicing (n=76)", "Unknown"))
MS = mutate(MS, Interactors=ifelse(MS$Gene %in% Int$V1, "Uniprot Interactors (11/26)", "Unknown"))

print("Check n=, for the following values correspond to the legend")
table(MS$Interactors)
table(MS$Function)

##############################################################################
#For Supplemental Table
# List the selected proteins by cut of and p value - raw data to 2 sig figs (in excel)
# subset from MS

MS = mutate(MS, Binders=ifelse(MS$Gene %in% RBin$V1, "RNA binding (n=34)", "Unknown"))

table(MS$Binders)

sig<-MS[MS$Legend == "p<0.05 (n=180)",]
sigcom<-sig[complete.cases(sig),]
sigcom$Binders[sigcom$Binders == "RNA binding (n=34)"]  <- "YES"
sigcom$Function[sigcom$Function == "RNA splicing (n=76)"] <- "YES"
sigcom$Interactors[sigcom$Interactors == "Uniprot Interactors (11/26)"] <- "YES"

################
print("check")

colnames(sigcom)[17]
colnames(sigcom)[17]<-"Interaction_with_CLK2_Uniprot"
colnames(sigcom)[18]
colnames(sigcom)[18]<-"Involved_in_RNA_Splicing"
colnames(sigcom)[19]
colnames(sigcom)[19]<-"Involved_in_RNA_Binding"

colnames(sigcom)[4]
colnames(sigcom)[4]<-"log2(NtermIP#1/ControlIP)"
colnames(sigcom)[5]
colnames(sigcom)[5]<-"log2(NtermIP#2/ControlIP)"
colnames(sigcom)[6]
colnames(sigcom)[6]<-"log2(CtermIP#1/ControlIP)"
colnames(sigcom)[7]
colnames(sigcom)[7]<-"log2(CtermIP#2/ControlIP)"
colnames(sigcom)[8]
colnames(sigcom)[8]<-"log2(CLK2IP/ControlIP)"
colnames(sigcom)[11]
colnames(sigcom)[11]<-"Adjusted_p_value(All_Expts)"
colnames(sigcom)[12]
colnames(sigcom)[12]<-"Unique_Peptides"
colnames(sigcom)[15]
colnames(sigcom)[15]<-"Mean(Log2(Ratios_of_all_Expts))"

# reorder and subset data using [ ]
suptab <- sigcom[c(1,2,3:8,15,11,12,17:19)]

write.table(suptab, file="SupplTable_Sig_Rawdata.tsv", col.names = TRUE, quote=FALSE, sep="\t")



#####
left=c("SRSF11","EIF4A3","CPSF1","PABPC1","CLK1","SRSF1","SRSF6","CLK3","SRSF10","SMN1",
       "LUC7L3","LUC7L2","CPSF1","MAGOHB")
right=c("DDX46","SRRM2","HNRNPC","PRPF8","SRSF2","CPSF7","CLK4","ELAVL1","SRSF3")
top=c("U2AF2","CLK2","NUDT21","CPSF6","U2AF1","SRRM1")

p = ggplot(MS, aes(all.mean, -log10(All_p.adj)), pch=20) +
  geom_point(aes(col=Legend)) +
#  scale_color_manual(values=c("light grey", "dark grey")) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% left,], aes(label=Gene),
  fontface = "bold", force = 7, nudge_x = -1.3, nudge_y=0.23) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% right,], aes(label=Gene),
  fontface = "bold", force = 7, nudge_x = 0.56, nudge_y=-0.23) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% top,], aes(label=Gene),
  fontface = "bold", force = 7, nudge_y=0.23) +
  geom_point(aes(col=Function)) +
#  geom_point(aes(col=Binders)) +
  geom_point(aes(col=Interactors)) +
#  scale_color_manual(values=c("grey90", "grey50","blue","dark green","red","transparent")) +
  scale_color_manual(values=c("grey90", "grey50","blue","red","transparent")) +
  annotate("text", x = 3.25, y = 0.25, label = "Involved in RNA splicing (selected)", fontface = "bold") +
  theme(legend.position = c(0.20, 0.85)) + ggtitle("CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "Mean (Log2 enrichment ratios of five experiments)", y="-log10(Adjusted p-value)") +
  geom_vline(xintercept=1, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")

##########################################
# plots all interactors to identify binders
o= ggplot(MS, aes(all.mean, -log10(All_p.adj)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% allbind,], aes(label=Gene))+  
  theme(legend.position = c(0.20, 0.85)) + ggtitle("CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "Mean (log2 enrichment ratios of five experiments)", y="-log10(Adjusted p-value)") +
  geom_vline(xintercept=1, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")


pdf(file="SupplFig_CLK2_Interactors_Splicing_V2.pdf")
p
dev.off()

p = ggplot(MS, aes(all.mean, -log10(All_p.adj)), pch=20) +
  geom_point(aes(col=Legend)) +
  scale_color_manual(values=c("light grey", "dark grey")) +
  theme_light()+geom_text_repel(data = MS[MS$Gene %in% left,], aes(label=Gene),
  fontface = "bold", force = 7, nudge_x = -1.3, nudge_y=0.23) +

##############################
# Submit this file to http://www.pantherdb.org/geneListAnalysis.do
# remove the "'" from 3' as it causes import failure
GOinq<-read.table(file="CLK2-bind_MS-GO.txt", header=TRUE, skip=10, sep="\t")
GOinr<-read.table(file="CLK2-bind_MS-GOMolFun.txt", header=TRUE, skip=10, sep="\t", quote='"')

colnames(GOinq)[1]<-"PantherGO"
colnames(GOinq)[6]<-"Fold_Enrichment"
colnames(GOinq)[7]<-"Adjusted.P.value"
GOq <- GOinq[order(-GOinq$Adjusted.P.value),]

colnames(GOinr)[1]<-"PantherGO"
colnames(GOinr)[6]<-"Fold_Enrichment"
colnames(GOinr)[7]<-"Adjusted.P.value"
GOr <- GOinr[order(-GOinr$Adjusted.P.value),]

# Set value >5 to 5
GOq$Enrichment <- as.numeric(gsub(" > ","", GOq$Fold_Enrichment))
GOq$PantherGO<-gsub("\\s*\\([^\\)]+\\)","",as.character(GOq$PantherGO))
GOr$Enrichment <- as.numeric(gsub(" > ","", GOr$Fold_Enrichment))
GOr$PantherGO<-gsub("\\s*\\([^\\)]+\\)","",as.character(GOr$PantherGO))

GOq$Label <- do.call(paste, c(GOq[c("Client.Text.Box.Input..178.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))
GOr$Label <- do.call(paste, c(GOr[c("CLK2.Interactors..180.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))

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

GOq$value<-as.numeric(-log10(GOq$Adjusted.P.value))
GOr$value<-as.numeric(-log10(GOr$Adjusted.P.value))

gt3q <-subset(GOq, Enrichment > 3)
gt3r <-subset(GOr, Enrichment > 3 & value > 10)
gt3q$lab<-wrap.labels(gt3q$PantherGO,20)
gt3r$lab<-wrap.labels(gt3r$PantherGO,20)

q<-ggplot(data=gt3q, aes(x=reorder(lab,value), y=value)) +
  geom_bar(width=0.8, stat="identity")+
  geom_text(data=gt3q, aes(x=lab, y=value, hjust=-0.15, label=Label)) +
  scale_y_continuous(limits=c(0,30)) +
#  ggtitle("CLK2 interactors by Gene Ontology") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  labs(y = "-log10(Adjusted p-value)", x="Top GO Biological Processes") +
  coord_flip()

r<-ggplot(data=gt3r, aes(x=reorder(lab,value), y=value)) +
  geom_bar(width=0.8, stat="identity")+
  geom_text(data=gt3r, aes(x=lab, y=value, hjust=-0.15, label=Label)) +
  scale_y_continuous(limits=c(0,67)) +
#  ggtitle("CLK2 interactors by Gene Ontology") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  labs(y = "-log10(Adjusted p-value)", x="Top GO molecular functions") +
  coord_flip()


pdf(file="SupplFig_CLK2_Interactors_byGO2.pdf")
library(gridExtra)
grid.arrange(q, r, nrow=2)
dev.off()




#################
# List genes by GO - supplemental table

GOG<-read.delim(file="CLK2_interactors_byGO.txt", header=FALSE)

# Only get the GO codes within brackets
match1<-gsub("\\s*\\^*.+\\(","",as.character(GOin$PantherGO[1]))
match<-gsub("\\)","",match1)

GOG[GOG$V13 %in% match,]




##################
# working on it have to change the layout if required
### plotting both side by side
library(gridExtra)
grid.arrange(p, q, ncol=2)

