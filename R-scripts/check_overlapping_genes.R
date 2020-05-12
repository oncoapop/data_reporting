## PURPOSE:
## To check all the genes with monotonic expression changes 
## compare them with the NMD transcripts which are monotonically increase in terms of PSI and expression

# This script reads all the excel files from the eIF4A3 paper GitHUb sync repo
#library(XLConnect)
#rnaseq<-readWorksheetFromFile(file=rnaseqfile,sheet = c("Sheet1"), header=TRUE)

#library(plyr)
#library(dplyr)
#library(RCurl)

# This is where all the relevant excel spreadsheets may lie on MOMAC14 (GitHub sync folder)
dir="/Users/dyap/Documents/Collboration_Projects/EIF4A3/EIF4A3splicing/supplementary/table"
setwd(dir)

# This is the monotonic expression files - increasing and decreasing
expr_dir="Gene_expression"
mono_inc_expr=paste(dir, expr_dir, "gene_pathways_monotonically_increasing_info.csv", sep="/")
mono_dec_expr=paste(dir, expr_dir, "gene_pathways_monotonically_decreasing_info.csv", sep="/")

# To be compared with NMD increasing based on expr and also in ratio
nmd_dir="NMD_expression"
nmd_expr=paste(dir, nmd_dir, "NMD_expression_pathways_monotonically_increasing_info.csv", sep="/")
ratio_expr=paste(dir, nmd_dir, "NMD_ratio_pathways_monotonically_increasing_info.csv", sep="/")

inc_expr<-read.csv(file=mono_inc_expr, header=TRUE, sep=",")
dec_expr<-read.csv(file=mono_dec_expr, header=TRUE, sep=",")
nmd<-read.csv(file=nmd_expr, header=TRUE, sep=",")
ratio<-read.csv(file=ratio_expr, header=TRUE, sep=",")

# There are all the genes names from each of the input files
# add the condition that they need to show monotonic expression changes in both cell lines with both drugs

LL1<-as.character(inc_expr[which(inc_expr$in_Hela_T_202_up=="1" & inc_expr$in_Hela_T_595_up=="1" & inc_expr$in_HCT116_T_202_up=="1" & inc_expr$in_HCT116_T_595_up=="1"), c("gene_name")] )
LL2<-c(LL1,as.character(dec_expr[which(dec_expr$in_Hela_T_202_up=="1" & dec_expr$in_Hela_T_595_up=="1" & dec_expr$in_HCT116_T_202_up=="1" & dec_expr$in_HCT116_T_595_up=="1"), c("gene_name")]))
LL3<-c(LL2,as.character(nmd[which(nmd$in_Hela_T_202_up=="1" & nmd$in_Hela_T_595_up=="1" & nmd$in_HCT116_T_202_up=="1" & nmd$in_HCT116_T_595_up=="1"), c("gene_name")]))
LL4<-c(LL3,as.character(ratio[which(ratio$in_Hela_T_202_up=="1" & ratio$in_Hela_T_595_up=="1" & ratio$in_HCT116_T_202_up=="1" & ratio$in_HCT116_T_595_up=="1"), c("gene_name")]))
LL<-unique(LL4)
univ<-as.data.frame(LL, stringsAsFactors = FALSE)

rec<-nrow(univ)

# defining the data frame with placeholder data
outdf <- data.frame(Gene = rep("", rec),
                  GOTerm = rep("", rec),
                Inc_Expr = rep("", rec),
                Dec_Expr = rep("", rec),
                  NMD_Up = rep("", rec),
                Ratio_Up = rep("", rec),
                     stringsAsFactors = FALSE)

# This is the module where I sum up everything from the various datasets
for (ri in seq(rec)) {

                gene<-univ[[1]][ri]

                if ( length(gene) > 0 )
                                {
                                 outdf$Gene[ri] <- gene

				test <- inc_expr[inc_expr[,6] %in% gene, c(6)]
				if (length(test) > 0 )
					{
      		                        outdf$Inc_Expr[ri] <- unique(as.character(inc_expr[inc_expr[,6] %in% gene, c("gene_name")]))
#					outdf$GOTerm[ri] <-  paste(unique(as.character(inc_expr[inc_expr[,c(2,3,4,5)] < 0.05 & inc_expr[,c(2,3,4,5)] !="NA", c(1)])), collapse=",")
					}

				test <- dec_expr[dec_expr[,6] %in% gene, c(6)]
				if (length(test) > 0 )
					{
      		                        outdf$Dec_Expr[ri] <- unique(as.character(dec_expr[dec_expr[,6] %in% gene, c("gene_name")]))
#					outdf$GOTerm[ri] <-  paste(unique(as.character(dec_expr[dec_expr[,6] %in% gene, c("GO_process")])), collapse=",")
					}

				test <- nmd[nmd[,6] %in% gene, c(6)]
				if (length(test) > 0 )
					{
      		                        outdf$NMD_Up[ri] <- unique(as.character(nmd[nmd[,6] %in% gene, c("gene_name")]))
#					outdf$GOTerm[ri] <-  paste(unique(as.character(nmd[nmd[,6] %in% gene, c("GO_process")])), collapse=",")
					}

				test <- ratio[ratio[,6] %in% gene, c(6)]
				if (length(test) > 0 )
					{
      		                        outdf$Ratio_Up[ri] <- unique(as.character(ratio[ratio[,6] %in% gene, c("gene_name")]))
#					outdf$GOTerm[ri] <-  paste(unique(as.character(ratio[ratio[,6] %in% gene, c("GO_process")])), collapse=",")
					}

                                }

			}

NMD_up_Expr_up <- outdf[which(outdf$Inc_Expr!='' & outdf$NMD_Up !='' & outdf$Ratio_Up!=''),c("Gene")]
NMD_up_Expr_down <- outdf[which(outdf$Dec_Expr!='' & outdf$NMD_Up !='' & outdf$Ratio_Up!=''),c("Gene")]

# This is to cross check with RMPs with motif changes
motifdir="/Users/dyap/Documents/Collboration_Projects/EIF4A3/Experiments/motifs"
motifile="mergedMotifs_proteins"

motifs=paste(motifdir, motifile, sep="/")

RBP<-read.table(file=motifs, header=TRUE, sep=" ", stringsAsFactors = FALSE)

dat<-as.data.frame(RBP)
nameList <- strsplit(dat[,2], split=",")
newdf <- data.frame(names=unlist(nameList))

# Monotonic Increasing RBPs
RBP_up<-unique(as.character(inc_expr[inc_expr[,6] %in% newdf$name, c("gene_name")]))

# Monotonic Decreasing RBPs
RBP_down<-unique(as.character(dec_expr[dec_expr[,6] %in% newdf$name, c("gene_name")]))

# Monotonic Increasing RBPs NMD
RBP_nmd<-unique(as.character(nmd[nmd[,6] %in% newdf$name, c("gene_name")]))

# Monotonic Increasing RBPs Inc ratio
RBP_ratio<-unique(as.character(ratio[ratio[,6] %in% newdf$name, c("gene_name")]))

# Specific NMD isoform (NMD expr + inclusion ratio, both monotonic)
RBP_NMD<-RBP_nmd[RBP_nmd %in% RBP_ratio]

# Specific NMD isoform (NMD expr + inclusion ratio, both monotonic) + MONOtonic UP
RBP_NMD_UP<-RBP_up[RBP_up %in% RBP_NMD]

# Specific NMD isoform (NMD expr + inclusion ratio, both monotonic) + MONOtonic DOWN
RBP_NMD_DOWN<-RBP_down[RBP_down %in% RBP_NMD]


# Output NMD
NMD_up_Expr_up
NMD_up_Expr_down



# RBP
RBP_up
RBP_down
RBP_NMD
RBP_NMD_UP
RBP_NMD_DOWN
