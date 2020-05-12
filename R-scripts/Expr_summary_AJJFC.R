# To plot the 48 hr summary only
# SRRM1/2 data sent for RNA-seq

library(reshape2)
library(ggplot2)

wd="/share/lustre/projects/takeda_splicing_inhibitor/gene_fusions/data/SRRM_KD"
file=paste(wd,"SRRM_KD_readthrough_filtered_ctrl_norm.tsv",sep="/")

QC<-read.table(file = file, stringsAsFactors = FALSE, sep="\t", header=TRUE)

boxplot(QC$ACTB_norm ~ QC$sample_id, outline=FALSE, las=2, cex=0.8)

sub="48 hr post siRNA treatment"
expr <- subset(QC, sample_id=="NT-48hr" | sample_id=="SRRM1-48hr" | 
	sample_id=="SRRM2-48hr" , 
	select=c(ACTB_norm,GAPDH_norm,TUBA1B_norm,sample_id))

# Process the data to be accepted by ggplot multiple plot function
df<-melt(expr, variable.name = "ID", value.name = "IV")
samples <- c("NT-48hr", "SRRM1-48hr", "SRRM2-48hr")
df$sample_id <- factor(df$sample_id, levels=samples)

outdir="/home/dyap/Projects/Takeda_T3/CG"
outfile="SRRM12_AELDU_boxplot_48hr.pdf"
pdffile=paste(outdir,outfile,sep="/")

pdf(file=pdffile, width=15, height=8)
 ggplot(df, aes(sample_id, IV)) +
        geom_boxplot() +
#	geom_jitter(width = 0.2) +
        facet_wrap(~ID, scales="free_y") +
        scale_y_log10() +
        labs(x="Sample", y="log(Normalized expression)") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()



