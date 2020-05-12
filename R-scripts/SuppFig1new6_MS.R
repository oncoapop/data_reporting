# Script to analyse MS data

# ON beast
library(plyr)
library(ggrepel)
library(ggplot2)

# On beast
wd="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
pat="CLK_resMat_PSM_GeneLevel_vsn_Vars_v02.csv"

##################
# 4x N-term FLAG #
#################
pat="CLK_resMat_PSM_GeneLevel_NoC_vsn_Vars_v03.csv"
setwd(wd)
file_names = list.files(pattern = pat)

MS<-as.data.frame(read.csv(file_names[1],header=TRUE))
MS<-read.csv(file_names[1],header=TRUE)
results<-read.csv(file_names[1],header=TRUE)

# Highlight known Uniprot/BioGrid/IntAct UBI database interactors of CLK2
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, sep="\t")

##########################################
### Version for Enrichment lg2 = 0 - V5 ##
##########################################

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

write.csv(allbinding, file="CLK2_interactors_MS_V6.csv", quote=FALSE)
clk2<-as.list(as.character(allbinding$Gene))
system("rm -f CLK2_int_genes_V6.txt")
lapply(clk2, write, "CLK2_int_genes_V6.txt", append=TRUE, ncolumns=nrow(allbinding))

# Download CLK2_interactors_GOList_V6 from PantherGO
system("grep \"splic\" CLK2_interactors_GOList_V6.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_inv_Splicing_V6.txt")

####################

# Highlight known Uniprot/BioGrid/IntAct UBI database interactors of CLK2
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, sep="\t")
# label all those with GO splicing

Splice<-read.table(file="CLK2_inv_Splicing_V5.txt", header=FALSE)
CG.enrich<-read.table(file="/home/dyap/Projects/Takeda_T3/CG_factors/CG_RNA_motifs_full", header=FALSE)

# Nice labelling (non-overlapping)
MS<-read.csv(file_names[1],header=TRUE)

MS<-mutate(MS, Legend=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0, "p<0.05 (n=244)", "Not Significant (n=984)"))
MS<-mutate(MS, Interactors=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Int$V1, "Known Interactors (27/97)", ""))
MS<-mutate(MS, Function=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1, "RNA splicing (n=92)", ""))
MS<-mutate(MS, Both=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1 & MS$Gene.name %in% Int$V1, "Both (n=23)", ""))
MS<-mutate(MS, Enriched=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% CG.enrich$V1, "CG Enriched (n=8)", ""))

table(MS$Legend)
table(MS$Interactors)
table(MS$Function)
table(MS$Both)
table(MS$Enriched)

# Label Uniprot binary interactions

p = ggplot(MS, aes(log2.average.treatment.effect.difference.from.control., 
           -log10(Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+
#  geom_text_repel(data = MS[MS$Interactors == "Known Interactors (27/97)",], aes(label=Gene.name), fontface = "bold") + 
  geom_label_repel(data = MS[MS$Enriched == "CG Enriched (n=8)",], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_y=5.53, nudge_x =-2.16) +
  geom_point(aes(col=Enriched)) + 
  geom_point(aes(col=Function)) + 
  geom_point(aes(col=Interactors)) + 
  scale_color_manual(values=c("transparent","darkgreen","red","grey90", "grey50","blue")) +
  theme(legend.position = c(0.18, 0.85)) + ggtitle("CLK2 interactors by IP-MS (Motif altered in CG, also known interactors)") +
#  theme(legend.position = c(0.18, 0.85)) + ggtitle("CLK2 interactors by IP-MS (known interactors involved in RNA splicing)") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
 # annotate("text", x = -0.5, y = 75, label = "Note the number of CG enriched", fontface = "bold") +
 # annotate("text", x = -0.5, y = 70, label = "factors below (x=1) cutoff", fontface = "bold") +
  annotate("text", x = 2.53, y = 147, label = "Factors with motifs enriched in CG", fontface = "bold") +
  geom_vline(xintercept=0, na.rm = FALSE, show.legend = TRUE, linetype = "dotted") +
  labs(x = "log2(Ave Effect Difference from Control)", y="-log10((Benjamini-Hochberg adjusted Wald p-value)") 

#pdf(file="SupplFig_CLK2_Interactors_RNAknownBinder_V6.pdf")
pdf(file="SupplFig_CLK2_Interactors_CGenriched_V6.pdf")
p
dev.off()


###########################
## FOR PAPER SUPP FIGURE ##
###########################

# Reload fresh MS data
MS<-as.data.frame(read.csv(file_names[1],header=TRUE))

# Download CLK2_interactors_GOList_V5 from PantherGO
system("grep \"splic\" CLK2_interactors_GOList_V6.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_inv_Splicing_V6.txt")

Splice<-read.table(file="CLK2_inv_Splicing_V6.txt", header=FALSE)
CG.enrich<-read.table(file="/home/dyap/Projects/Takeda_T3/CG_factors/CG_RNA_motifs_full", header=FALSE)

# Highlight known interactors of CLK2 (UBI databases)
Int<-read.table(file="all_UBIdb_interactors.txt", header=FALSE, skip=0, sep="\t")
Uniprot<-read.table(file="Uniprot_Interactors.txt", header=FALSE, skip=0, sep="\t")
Biogrid<-read.table(file="BioGrid_Interactors.txt", header=FALSE, skip=0, sep="\t")
IntAct<-read.table(file="IntAct_Interactors.txt", header=FALSE, skip=1, sep="\t")

# grep "RNA binding" CLK2_interactors_GOList.txt | awk -F"\t" '{print $2}' | awk -F";" '{print $2}' > new_CLK2_int_mRNAbinding.txt
system("grep \"RNA binding\" CLK2_interactors_GOList_V6.txt | awk -F\"\t\" '{print $2}' | awk -F\";\" '{print $2}' > CLK2_inv_mRNAbinding_V6.txt")

RBin<-read.table(file="CLK2_inv_mRNAbinding_V6.txt", header=FALSE)

##########
MS<-mutate(MS, Legend=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0, "p<0.05 & Enriched (n=244)", "Filtered (n=984)"))
MS<-mutate(MS, Interactors=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Int$V1, "Known Interactors (27/97)", ""))
MS<-mutate(MS, Function=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% Splice$V1, "RNA splicing (n=95)", ""))
MS<-mutate(MS, Unpi=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Uniprot$V1, 
	"Uniprot Interactors (n=10/26)", ""))
MS<-mutate(MS, Bigi=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% Biogrid$V1, 
	"Biogrid Interactors (n=26/69)", ""))
MS<-mutate(MS, Inai=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% IntAct$V1, 
	"IntAct Interactors (n=13/52)", ""))
MS<-mutate(MS, Binders=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & MS$Gene.name %in% RBin$V1, "RNA binding (n=139)", ""))
MS<-mutate(MS, Enriched=ifelse(MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 & 
	MS$log2.average.treatment.effect.difference.from.control.>0 & 
	MS$Gene.name %in% CG.enrich$V1, "Motif change in CG (n=8)", ""))

print("Check n=, for the following values correspond to the legend")
table(MS$Legend)
table(MS$Interactors)
table(MS$Function)
table(MS$Unpi)
table(MS$Bigi)
table(MS$Inai)
table(MS$Binders)
table(MS$Enriched)

##############################################################################
#For Supplemental Table
# List the selected proteins by cut of and p value - raw data to 2 sig figs (in excel)

sig<-MS[MS$Legend == "p<0.05 & Enriched (n=244)",]
sigcom<-sig[complete.cases(sig),]
sigcom$Binders[sigcom$Binders == "RNA binding (n=139)"]  <- "YES"
sigcom$Function[sigcom$Function == "RNA splicing (n=95)"] <- "YES"
sigcom$Interactors[sigcom$Interactors == "Known Interactors (27/97)"] <- "YES"
sigcom$Unpi[sigcom$Unpi == "Uniprot Interactors (n=10/26)"] <- "YES"
sigcom$Bigi[sigcom$Bigi == "Biogrid Interactors (n=26/69)"] <- "YES"
sigcom$Inai[sigcom$Inai == "IntAct Interactors (n=13/52)"] <- "YES"
sigcom$Enriched[sigcom$Enriched == "Motif change in CG (n=8)"] <- "YES"

################
names(sigcom)
print("check")

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
colnames(sigcom)[34]<-"Binding_Motif_change_in_CG"

# reorder and subset data using [ ]
suptab <- sigcom[c("Gene.name","Accession","Descriptions","Sequence","Annotated.Sequence","Modifications",
	"Absolute.fold.change.direction","Fold.change.for.average.treatment.effect",
	"log2.average.treatment.effect.difference.from.control.","Average.treatment.effect.t.statistic.Wald.based.p.value",
	"Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect",
	"Binding_Motif_change_in_CG","Known_Interactors","Uniprot_Interactors","BioGrid_Interactors","IntAct_Interactors",
	"Involved_in_RNA_Splicing","Involved_in_RNA_binding")]

write.table(suptab, file="SupplTable_IP-MS_CLK2_V6.tsv", col.names = TRUE, quote=FALSE, sep="\t")



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

left=c("SNRNP70","DHX15")
leftlab=c("SRSF7","ELAVL1","SRSF1","SRSF6","KHDRBS1","ZNF638","HNRNPH1","SFPQ")
leftlab1=c("U2AF2")
right=c("DDX46","CPSF7","SRSF3","TRA2B","U2AF1","YTHDC1","CLK1",
	"RNPS1", "TNPO3","TRA2A","SON","NUDT21","LUC7L","CPSF4")
top=c("CLK2","THRAP3","SRRM1","ACIN1","BCLAF1","FIP1L1",
	"LUC7L3","RBM39","LUC7L2","DDX5","CLK3","PNN","SRRM2")
below=c("PRPF38A","CPSF2","CPSF1")

p = ggplot(MS, aes(log2.average.treatment.effect.difference.from.control., 
           -log10(Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% left,], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_x = -0.3, nudge_y=-1.04) +
  theme_light()+geom_label_repel(data = MS[MS$Gene.name %in% leftlab1,], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_x = -0.9) +
  theme_light()+geom_label_repel(data = MS[MS$Gene.name %in% leftlab,], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_y=c(30,-5,15,20,5,30,40), nudge_x =-1.56) +
  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% right,], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_x = 0.56, nudge_y=0.53) +
  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% top,], aes(label=Gene.name),
  fontface = "bold", force = 7, nudge_y=0.23) +
#  theme_light()+geom_label_repel(data = MS[MS$Gene.name %in% toplab,], aes(label=Gene.name),
#  fontface = "bold", force = 7, nudge_y=0.23) +
  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% below,], aes(label=Gene.name),
  fontface = "bold", force = 3, nudge_y=c(-10,-5,-12,-2), nudge_x = c(0.6,0.2,-0.2,0)) +
#  theme_light()+geom_label_repel(data = MS[MS$Gene.name %in% belowlab,], aes(label=Gene.name),
#  fontface = "bold", force = 7, nudge_y=-4.23) +
  scale_y_continuous(limits=c(-10,50)) +
  geom_point(aes(col=Function)) +
  geom_point(aes(col=Interactors)) +
  scale_color_manual(values=c("transparent","grey90","red", "grey50","blue")) +
#  annotate("text", x = 2.53, y = 147, label = "Involved in RNA splicing (selected)", fontface = "bold") +
  theme(legend.position = c(0.17, 0.75)) + ggtitle("CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "log2(Ave Effect Difference from Control)", y="-log10(Benjamini-Hochberg adjusted Wald p-value)") +
  geom_vline(xintercept=0, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")

pdf(file="SupplFig_CLK2_Interactors_Splicing_V6.pdf", useDingbats=FALSE)
p
dev.off()


##############################


##########################################
# plots all interactors to identify binders
o = ggplot(MS, aes(log2.average.treatment.effect.difference.from.control., 
           -log10(Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect)), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+geom_label_repel(data = MS[MS$Benjamini.Hochberg.adjusted.Wald.p.value.for.average.treatment.effect<0.05 &
        MS$log2.average.treatment.effect.difference.from.control.>1,], aes(label=Gene.name), cex=2.5, force=7)+  
#  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% Splice$V1,], aes(label=Gene.name), force=10)+  
#  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% RBin$V1,], aes(label=Gene.name), force=10)+  
#  theme_light()+geom_text_repel(data = MS[MS$Gene.name %in% Int$V1,], aes(label=Gene.name), force=10)+  
  theme(legend.position = c(0.20, 0.85)) + ggtitle("All significant CLK2 interactors by IP-MS") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "log2(Ave Effect Difference from Control)", y="-log10(Benjamini-Hochberg adjusted Wald p-value)") +
  geom_vline(xintercept=1, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")
o
pdf(file="ReferenceFig_CLK2_Interactors_V3all.pdf")
o
dev.off()

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

