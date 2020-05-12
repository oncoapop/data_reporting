# Script to analyse MS data

# ON beast
library(plyr)
library(ggrepel)
library(ggplot2)

# On beast
wd="/home/dyap/Projects/PPP2R2A/ProteoGenomics_Nat2016"
pat="CPTAC_BC_SupplementaryTable09.csv"
setwd(wd)
file_names = list.files(pattern = pat)

# From Nature 2016 (ProteoGenomics May 2016 Paper)

##########################################
### Version for Initial assessment  ##
##########################################

###########################
## FOR PAPER SUPP FIGURE ##
###########################

# Reload fresh data from Suppl table 09
Corr<-as.data.frame(read.csv(file_names[1],header=TRUE, stringsAsFactors=FALSE))
#Corr<-dat[complete.cases(dat),]

highlight=c("PPP2R2A","SPDEF","PDEF","GATA3","MTAP","MAP2K4","ZNF703","PTEN","MYC","CCND1","MDM2","ERBB2","CCNE1","ATM",
	"CLK2","EIF4A3")
#highlight=c("APOBEC1","APOBEC2","APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G","APOBEC3H","CDA","CDADC1",
#	"CTNNB1","ATR","FANCL","FANCD","RAD50","NF1","RPL22","POLE","FOXL2","KMT2B","PER3")
paper=c("CDK12","PAK1","PTK2","RIPK2","TLK2","ESR1","PGR")

##########
Corr$m.lg.Adj_p.value<--log10(as.numeric(gsub("< ","", Corr$FDR.p.value)))
Corr$RP_Corr<-as.numeric(gsub("< ","", Corr$Pearson.corrrelation))

Corr<-mutate(Corr, Legend=ifelse(Corr$m.lg.Adj_p.value<1.30102999566, "Significant (n=3143)", "Not significant (n=6159)"))
Corr<-mutate(Corr, high=ifelse(Corr$Gene.ID %in% highlight, "Interest", ""))
Corr<-mutate(Corr, pap=ifelse(Corr$Gene.ID %in% paper, "Paper", ""))
table(Corr$Legend)
table(Corr$high)
table(Corr$pap)

##############
# FOR Figure #
##############

#############
# Nice labelling (non-overlapping)

p = ggplot(Corr, aes(RP_Corr,m.lg.Adj_p.value), pch=20) +
  geom_point(aes(col=Legend)) +
  theme_light()+geom_label_repel(data = Corr[Corr$Gene.ID %in% highlight,], aes(label=Gene.ID),
  fontface = "bold", force = 5, nudge_x = -0.3, nudge_y=-0.3)+
  geom_point(aes(col=high)) +
  theme_light()+geom_label_repel(data = Corr[Corr$Gene.ID %in% paper,], aes(label=Gene.ID),
  fontface = "bold", force = 5, nudge_x = c(0.1,-0.1,0.1,-0.1,-0.2,-0.15,0.2), nudge_y=c(-0.1,-0.1,-0.1,0.1,-0.2,-0.1,0.2), 
	color = "white", fill = "red")+
	geom_point(aes(col=pap)) +
  scale_color_manual(values=c("transparent","blue","grey40","red","grey80")) +
  theme(legend.position = c(0.17, 0.75)) + ggtitle("mRNA-to-protein correlation analysis (SuppTab09, Nat 2016)") +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "Pearson's corrrelation", y="-log10(FDR p-value)") +
  geom_vline(xintercept=0, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")+
  geom_hline(yintercept=1.30102999566, na.rm = FALSE, show.legend = TRUE, linetype = "dotted")

p


######

pdf(file="SupplFig9_Nat2016_RNA-Prot-Correction.pdf", useDingbats=FALSE)
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


##############################################################################
#For Supplemental Table
# List the selected proteins by cut of and p value - raw data to 2 sig figs (in excel)

sig<-MS[MS$Legend == "Above threshold (n=270)",]
sigcom<-sig[complete.cases(sig),]
sigcom$Binders[sigcom$Binders == "RNA binding (n=149)"]  <- "YES"
sigcom$Function[sigcom$Function == "RNA splicing (n=98)"] <- "YES"
sigcom$Interactors[sigcom$Interactors == "Known Interactors (28/97)"] <- "YES"
sigcom$Unpi[sigcom$Unpi == "Uniprot Interactors (n=10/26)"] <- "YES"
sigcom$Bigi[sigcom$Bigi == "Biogrid Interactors (n=27/69)"] <- "YES"
sigcom$Inai[sigcom$Inai == "IntAct Interactors (n=13/52)"] <- "YES"
sigcom$Enriched[sigcom$Enriched == "Motif change in CG (n=10)"] <- "YES"

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

write.table(suptab, file="SupplTable_V5.tsv", col.names = TRUE, quote=FALSE, sep="\t")




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

