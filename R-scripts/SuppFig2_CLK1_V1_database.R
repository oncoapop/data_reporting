# Script to analyse MS data

# ON beast
library(plyr)
library(ggrepel)
library(ggplot2)

# On beast
wd="/home/dyap/Projects/Takeda_T3/MS data/Expt3_CLK1"
pat="CLK_resMat_PSM_GeneLevel_vsn_Vars_v02.csv"
setwd(wd)
file_names = list.files(pattern = pat)

MS<-as.data.frame(read.csv(file_names[1],header=TRUE))
MS<-read.csv(file_names[1],header=TRUE)
results<-read.csv(file_names[1],header=TRUE)

# Highlight known Uniprot/BioGrid/IntAct UBI database interactors of CLK2
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, sep="\t")

##############################################
###   Version for Enrichment lg2 = 0 - V7  ###
### Highlight 3' processing factors in red ###
##############################################

## PRELIM

allbinding=subset(MS, Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect < 0.05 & 
	log2.average.treatment.effect.difference.from.control.> 0, select=c(Accession,Gene.name,Descriptions))
#,pepNum))


# Nice labelling (non-overlapping)
MS<-mutate(MS, Legend=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0, "p<0.05 (n=270)", "Not Significant (n=957)"))
MS<-mutate(MS, Interactors=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Int$V1, "Known Interactors (28/97)", "Unknown"))
#MS<-mutate(MS, Interactors=ifelse(MS$Gene.name %in% Int$V1, "Known Interactors (34/97)", "Unknown"))


count(MS$Legend)
count(MS$Interactors)
nrow(allbinding)

write.csv(allbinding, file="CLK2_interactors_MS_V5.csv", quote=FALSE)
clk2<-as.list(as.character(allbinding$Gene))
system("rm -f CLK2_int_genes_V5.txt")
lapply(clk2, write, "CLK2_int_genes_V5.txt", append=TRUE, ncolumns=nrow(allbinding))
# Submit this file to http://www.pantherdb.org/geneListAnalysis.do
# remove the "'" from 3' as it causes import failure
GO<-read.table(file="CLK2-bind_MS-GO_new_V5.txt", header=TRUE, skip=10, sep="\t")

####################

# Highlight known Uniprot/BioGrid/IntAct UBI database interactors of CLK2
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, sep="\t")
# label all those with GO splicing

Splice<-read.table(file="CLK2_inv_Splicing_V5.txt", header=FALSE)
CG.enrich<-read.table(file="/home/dyap/Projects/Takeda_T3/CG_factors/CG_RNA_motifs_full", header=FALSE)
EndPro <-read.table(file="CLK2_3endproc_V7.txt", header=FALSE)

# Nice labelling (non-overlapping)
MS<-read.csv(file_names[1],header=TRUE)

MS<-mutate(MS, Legend=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0, "Above threshold (n=270)", "Below threshold (n=957)"))
MS<-mutate(MS, Interactors=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Int$V1, "Known Interactors (28/97)", ""))
MS<-mutate(MS, Function=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1, "RNA splicing (n=98)", ""))
MS<-mutate(MS, Both=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1 & MS$Gene.name %in% Int$V1, "Both (n=24)", ""))
MS<-mutate(MS, EndProc=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% EndPro$V1, "3'-end Processing (n=31)", ""))
MS<-mutate(MS, Enriched=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% CG.enrich$V1, "CG Enriched (n=10)", ""))

table(MS$Legend)
table(MS$Interactors)
table(MS$Function)
table(MS$Both)
table(MS$Enriched)
table(MS$EndProc)

# Plot all important genes for later filtering and arrangement

p = ggplot(MS, aes(log2.average.treatment.effect.difference.from.control., 
           -log10(Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+
  geom_label_repel(data = MS[(MS$Enriched == "CG Enriched (n=10)" & MS$EndProc != "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "bold", force = 7, nudge_y=c(20,-5,-5,5,5,23), nudge_x =c(-1.5,-0.15,-1.2,-1.2,-0.1,-1.23), color="black") +
  geom_label_repel(data = MS[(MS$Enriched == "CG Enriched (n=10)" & MS$EndProc == "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "bold", force = 7, nudge_y=c(15,35,33,10), nudge_x =c(-1.75,-1.3,-1.6,0), color="red") +
  geom_text_repel(data = MS[(MS$Enriched != "CG Enriched (n=10)" & MS$EndProc == "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "plain", force = 7, nudge_y=70, nudge_x =-0.5, color="red") +
  geom_text_repel(data = MS[(MS$Both ==  "Both (n=24)"  & MS$EndProc != "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "plain", force = 7, nudge_y=10, nudge_x =1, color="blue") +

  geom_point(aes(col=Function)) + 
  geom_point(aes(col=Interactors)) + 
  geom_point(aes(col=EndProc)) + 
  geom_point(aes(col=Enriched)) + 
  scale_color_manual(values=c("transparent","red","grey50", "grey90","orange","blue","darkgreen")) +
  theme(legend.position = c(0.18, 0.75)) + ggtitle("CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  annotate("text", x = 2.53, y = 147, label = "Factors with motifs enriched in CG", fontface = "bold") +
  geom_vline(xintercept=0, na.rm = FALSE, show.legend = TRUE, linetype = "dotted") +
  labs(x = "log2(Mean Enrichment)", y="-log10((Benjamini-Hochberg adjusted Wald p-value)") 

p

#pdf(file="SupplFig_CLK2_Interactors_RNAknownBinder.pdf")
pdf(file="SupplFig_CLK2_Interactors_invol_splice_3end.pdf")
p
dev.off()


###########################
## FOR PAPER IP-MS FIGURE ##
###########################

## LOAD libraries in header

# Reload fresh MS data
MS<-as.data.frame(read.csv(file_names[1],header=TRUE))

# Download CLK2_interactors_GOList_V5 from PantherGO
system("grep \"splic\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_inv_Splicing_V5.txt")

Splice<-read.table(file="CLK2_inv_Splicing_V5.txt", header=FALSE)
CG.enrich<-read.table(file="/home/dyap/Projects/Takeda_T3/CG_factors/CG_RNA_motifs_full", header=FALSE)

# Highlight known interactors of CLK2 (UBI databases)
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, skip=0, sep="\t")
Uniprot<-read.table(file="Uniprot_Interactors.txt", header=FALSE, skip=0, sep="\t")
Biogrid<-read.table(file="BioGrid_Interactors.txt", header=FALSE, skip=0, sep="\t")
IntAct<-read.table(file="IntAct_Interactors.txt", header=FALSE, skip=1, sep="\t")

# grep "RNA binding" CLK2_interactors_GOList.txt | awk -F"\t" '{print $2}' | awk -F";" '{print $2}' > new_CLK2_int_mRNAbinding.txt
system("grep \"RNA binding\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_inv_mRNAbinding_V5.txt")

RBin<-read.table(file="CLK2_inv_mRNAbinding_V5.txt", header=FALSE)

# grep "end-processing" CLK2_interactors_GOList_V5.txt | awk -F"\t" '{print $2}' | awk -F";" '{print $2}' > CLK2_3endproc_V7.txt
system("grep \"end-processing\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_3endproc_V7.txt")
system("grep \"end processing\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' >> CLK2_3endproc_V7.txt")
system("grep \"poly(A) binding\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' >> CLK2_3endproc_V7.txt")
system("grep \"termination of RNA polymerase II transcription\" CLK2_interactors_GOList_V5.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_Pol2term_V7.txt")

EndPro <-read.table(file="CLK2_3endproc_V7.txt", header=FALSE)
Pol2Term<- read.table(file="CLK2_Pol2term_V7.txt", header=FALSE)

##########

MS<-mutate(MS, Legend=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0, "Above threshold (n=270)", "Below threshold (n=957)"))
MS<-mutate(MS, Interactors=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Int$V1, "Known Interactors (28/97)", ""))
MS<-mutate(MS, Function=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1, "RNA splicing (n=98)", ""))
MS<-mutate(MS, Unpi=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Uniprot$V1, 
	"Uniprot Interactors (n=10/26)", ""))
MS<-mutate(MS, Bigi=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Biogrid$V1, 
	"Biogrid Interactors (n=27/69)", ""))
MS<-mutate(MS, Inai=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% IntAct$V1, 
	"IntAct Interactors (n=13/52)", ""))
MS<-mutate(MS, Binders=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% RBin$V1, "RNA binding (n=149)", ""))
MS<-mutate(MS, EndProc=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% EndPro$V1, "3'-end Processing (n=31)", ""))
MS<-mutate(MS, Enriched=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% CG.enrich$V1, "Motif density change in CG (n=10)", ""))
#MS<-mutate(MS, Transterm=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
#	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Pol2Term$V1, "termination of RNA polymerase II transcription (n=xx)", ""))
MS<-mutate(MS, Z.both=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 &
        MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% CG.enrich$V1 & MS$Gene.name %in% EndPro$V1, "Motif change & 3'-end Processing (n=4)", ""))

print("Check n=, for the following values correspond to the legend")
table(MS$Legend)
table(MS$Interactors)
table(MS$Function)
table(MS$Unpi)
table(MS$Bigi)
table(MS$Inai)
table(MS$Binders)
table(MS$EndProc)
#table(MS$Transterm)
table(MS$Enriched)
table(MS$Z.both)

##############################################################################
#For Supplemental Table
# List the selected proteins by cut of and p value - raw data to 2 sig figs (in excel)

sig<-MS[MS$Legend == "Above threshold (n=270)",]
sigcom<-sig[complete.cases(sig),]
sigcom$EndProc[sigcom$EndProc == "3'-end Processing (n=31)"]  <- "YES"
sigcom$Binders[sigcom$Binders == "RNA binding (n=149)"]  <- "YES"
sigcom$Function[sigcom$Function == "RNA splicing (n=98)"] <- "YES"
sigcom$Interactors[sigcom$Interactors == "Known Interactors (28/97)"] <- "YES"
sigcom$Unpi[sigcom$Unpi == "Uniprot Interactors (n=10/26)"] <- "YES"
sigcom$Bigi[sigcom$Bigi == "Biogrid Interactors (n=27/69)"] <- "YES"
sigcom$Inai[sigcom$Inai == "IntAct Interactors (n=13/52)"] <- "YES"
sigcom$Enriched[sigcom$Enriched == "Motif density change in CG (n=10)"] <- "YES"

names(sigcom)
print("check")

####################################

colnames(sigcom)[28]
colnames(sigcom)[28]<-"Known_Interactors"
colnames(sigcom)[29]
colnames(sigcom)[29]<-"Involved_in_RNA_Splicing"
colnames(sigcom)[30]
colnames(sigcom)[30]<-"Uniprot_Interactors"
colnames(sigcom)[31]
colnames(sigcom)[31]<-"BioGrid_Interactors"
colnames(sigcom)[32]
colnames(sigcom)[32]<-"IntAct_Interactors"
colnames(sigcom)[33]
colnames(sigcom)[33]<-"Involved_in_RNA_binding"
colnames(sigcom)[34]
colnames(sigcom)[34]<-"Involved_in_3'-end_Processing"
colnames(sigcom)[35]
colnames(sigcom)[35]<-"Binding_Motif_Density_change_in_CG"
colnames(sigcom)[36]
colnames(sigcom)[36]<-"Motif_change_CG_&_Poly(A)_3'_Processing"

# reorder and subset data using [ ]
suptab <- sigcom[c("Gene.name","Accession","Descriptions","Sequence","Annotated.Sequence","Modifications",
	"Absolute.fold.change.direction","Fold.change.for.average.treatment.effect",
	"log2.average.treatment.effect.difference.from.control.","Average.treatment.effect.t.statistic.Wald.based.p.value",
	"Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect",
	"Binding_Motif_Density_change_in_CG","Known_Interactors","Uniprot_Interactors","BioGrid_Interactors","IntAct_Interactors",
	"Involved_in_RNA_Splicing","Involved_in_RNA_binding","Involved_in_3'-end_Processing")]

write.table(suptab, file="SupplTable_V7.tsv", col.names = TRUE, quote=FALSE, sep="\t")


##############
# FOR Figure #
##############
#MS<-as.data.frame(read.csv(file_names[1],header=TRUE))
#MS<-read.csv(file_names[1],header=TRUE)

# Highlight known Uniprot/BioGrid/IntAct UBI database interactors of CLK2
#Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, sep="\t")

#############
# Nice labelling (non-overlapping)
# inherit from Supp Table

label3end=c("CPSF4","DDX3X","EIF4A3","RBM8A","SRSF3","SRSF5","SRSF7","SRSF11","NUDT21","SRRM2")
highlight=c("CLK2")
highlight2=c("CLK1")

p = ggplot(MS, aes(log2.average.treatment.effect.difference.from.control., 
           -log10(Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+

# CG enriched factors (label in not involed in end processing  ELAVL1  HNRNPC  HNRNPH1 MATR3   SFPQ    ZNF638)
  geom_label_repel(data = MS[(MS$Enriched == "Motif density change in CG (n=10)" & MS$EndProc != "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "bold", force = 3, nudge_y=c(40,1,10,8,5,30), nudge_x =c(-2.2,-2,-1.8,-1.7,-0.1,-2.1), color="black") +

# CG enriched factors (label in red if invol in 3' end processing -  KHDRBS1 SRSF1   SRSF6   U2AF2 )
  geom_label_repel(data = MS[(MS$Enriched == "Motif density change in CG (n=10)" & MS$EndProc == "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "bold", force = 3, nudge_y=c(20,35,33,30), nudge_x =c(-1.8,-1.3,-1.6,-0.02), color="red") +

# manual label those involved in 3' end processing (selected) - DDX3X  EIF4A3 NUDT21 RBM8A  SRRM2  SRSF11 SRSF3  SRSF5 SRSF7
  geom_text_repel(data = MS[(MS$Gene.name %in% label3end),], aes(label=Gene.name), fontface = "bold", 
  force = 3, nudge_y=c(-10,10,25,5,-10,5,5,-5,38,7), nudge_x =c(1.8,-0.4,0.01,0.5,2,-0.02,0.35,0.8,-2.1,-0.05) , color="red") +

#  geom_text_repel(data = MS[(MS$Enriched != "Motif density change in CG (n=10)" & MS$EndProc == "3'-end Processing (n=31)" & MS$Interactors != "Known Interactors (28/97)"),], 
#  aes(label=Gene.name), fontface = "plain", force = 7, nudge_y=70, nudge_x =-0.5, color="red") +

# Involved in splicing and known interactors (dark green) (ACIN1   BCLAF1  CLK3    DDX5    HNRNPM  LUC7L   LUC7L2  PNN    
#							 PRPF38A   RBM39   SNRNP70 SON     THRAP3  TNPO3   TRA2A   YTHDC1)
  geom_text_repel(data = MS[(MS$Enriched != "Motif density change in CG (n=10)" & MS$EndProc != "3'-end Processing (n=31)" & !(MS$Gene.name %in% highlight) 
    & MS$Interactors == "Known Interactors (28/97)" & MS$Function == "RNA splicing (n=98)"),], aes(label=Gene.name), fontface = "bold", 
    force = 3, nudge_y=c(0.01,10,8,5,-10,10,10,10,-5,15,15,-5,-5,15,-4,-10), nudge_x =c(-0.2,0.2,1.2,0.05,-1,1.2,-0.2,0.4,-0.7,0.01,-0.4,-0.1,-0.2,0.6,0.8,1.2),   
    color="black") +

# Known interactors involved in 3'-end processing (CPSF1  CPSF2  CPSF7  FIP1L1 RNPS1  SRRM1  U2AF1  WDR33)
  geom_text_repel(data = MS[(MS$Interactors == "Known Interactors (28/97)"  & MS$EndProc == "3'-end Processing (n=31)"),], aes(label=Gene.name),  
  fontface = "bold", force = 3, nudge_y=c(-5,-5,0,20,-2,5,-2,-20), nudge_x =c(0.2,0.1,0.2,0.1,1.2,0.2,0.4,0.7), color="red") +

# highlights
  geom_text_repel(data = MS[(MS$Gene.name %in% highlight2),], aes(label=Gene.name),  
  fontface = "bold", force = 3, nudge_y=c(-2), nudge_x =c(-0.15), color="black") +
  geom_label_repel(data = MS[(MS$Gene.name %in% highlight),], aes(label=Gene.name),  
  fontface = "bold", force = 3, nudge_y=c(10), nudge_x =c(-0.02), color="black", box.padding = unit(0.25, "lines"), fill="orange") +

  geom_point(aes(col=Function)) + 
  geom_point(aes(col=EndProc)) + 
  geom_point(aes(col=Interactors)) + 
  geom_point(aes(col=Enriched)) + 
  scale_color_manual(values=c("transparent","red","grey50", "grey90","orange","blue","darkgreen")) +
  theme(legend.position = c(0.30, 0.85)) + ggtitle("CLK2 associated factors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "log2(Mean enrichment)", y="-log10(Benjamini-Hochberg adjusted Wald p-value)") +
  geom_vline(xintercept=0, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")


pdf(file="Fig_CLK2_Interactors_Splicing_V7.pdf", useDingbats=FALSE)
p
dev.off()


##############################


#########################################
library(ggplot2)
# On beast
wd="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
setwd(wd)

# Submit this file to http://www.pantherdb.org/geneListAnalysis.do
# remove the "'" from 3' as it causes import failure
GOinq<-read.table(file="CLK2-bind_MS-GO_V5.txt", header=TRUE, skip=10, sep="\t")
GOinr<-read.table(file="CLK2-bind_MS-GOMolFun_V5.txt", header=TRUE, skip=10, sep="\t", quote='"')

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

GOq$Label <- do.call(paste, c(GOq[c("CLK2_int_genes_V5..271.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))
GOr$Label <- do.call(paste, c(GOr[c("CLK2_int_genes_V5..271.", "Homo.sapiens...REFLIST..20814.")], sep = "/"))

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

GOq$value<-round(as.numeric(-log10(GOq$Adjusted.P.value)),1)
GOr$value<-round(as.numeric(-log10(GOr$Adjusted.P.value)),1)

gt3q <-subset(GOq,  Adjusted.P.value < 0.05 & Enrichment > 23.5 & Enrichment != "Inf")
gt3r <-subset(GOr,  Adjusted.P.value < 0.05 & Enrichment != "Inf")
gt3q$lab<-wrap.labels(gt3q$PantherGO,20)
gt3r$lab<-wrap.labels(gt3r$PantherGO,20)

gt3q$Fold_Enrichment<-round(gt3q$Enrichment,1)

q<-ggplot(data=gt3q, aes(x=reorder(lab,value), y=value)) +
  geom_bar(width=0.8, stat="identity")+
  geom_text(data=gt3q, aes(x=lab, y=value, hjust=-0.15, label=Label)) +
  geom_text(data=gt3q, aes(x=lab, y=value, hjust=1.2, label=Fold_Enrichment), color="white") +
  scale_y_continuous(limits=c(0,28)) +
#  ggtitle("CLK2 interactors by Gene Ontology") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  labs(y = "-log10(Adjusted p-value)", x="Top GO Biological Processes by Fold Enrichment") +
  coord_flip()

r<-ggplot(data=gt3r, aes(x=reorder(lab,value), y=value)) +
  geom_bar(width=0.8, stat="identity")+
  geom_text(data=gt3r, aes(x=lab, y=value, hjust=-0.15, label=Label)) +
  scale_y_continuous(limits=c(0,100)) +
#  ggtitle("CLK2 interactors by Gene Ontology") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  labs(y = "-log10(Adjusted p-value)", x="Top GO molecular functions") +
  coord_flip()


pdf(file="SupplFig_CLK2_Interactors_byGO_BP_V5.pdf", useDingbats=FALSE)
q
dev.off()

pdf(file="SupplFig_CLK2_Interactors_byGO_MF_V5.pdf")
p
dev.off()

#library(gridExtra)
#grid.arrange(q, r, nrow=2)




