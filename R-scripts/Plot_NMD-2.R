# Script to plot data from NMD expression correlation

# Load R libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
#library(plotly) # Pleione only

##########################################################
### CHANGE THESE PARAMETERS FOR EACH RUN OF THE SCRIPT ###
##########################################################
 
# Inputs
file="HCT116_T_202_normalized_highest_expression_HCT116_T_228_normalized.res"
#file="HCT116_T_595_normalized_highest_expression_HCT116_T_598_normalized.res"
#file="HCT116_T_202_normalized_longest_HCT116_T_228_normalized.res"
#file="HCT116_T_595_normalized_longest_HCT116_T_598_normalized.res"
#file="Hela_T_202_normalized_longest_Hela_T_228_normalized.res"
#file="Hela_T_595_normalized_longest_Hela_T_598_normalized.res"
dir="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation/other_compounds"

# Outputs
design="Normalized_NMD_Expr_vs_Highest_Isoform"
#design="Normalized_NMD_Expr_vs_Longest_Isoform"
#cell_line="HeLa"
cell_line="HCT-116"
#drug="T-202"
#drug="T-228"
drug="T-595"
#drug="T-598"
series="0, 0.5, 2.0, 5.0, 10, 20"
Expt<-paste(paste(design,cell_line,sep="_"),drug,sep="_")

##########################################################

outdir="/home/dyap/Projects/eIF4A3/Plots"
fname=paste(outdir,Expt,sep="/")
pdfname=paste(fname,"pdf",sep=".")
graphtitle=gsub("_"," ",design)
makepdf = TRUE

setwd(dir)
data<-read.table(file=file, header=TRUE)
#api<-"afuTYWB3jSd30tdqJmPY"
#user<-"oncoapop"
#py <- plotly(username=user, key=api)

NMD <- data[c(1,6:11)]
Con <- data[c(1,12:17)]
Expr <- data[c(1,18:23)]

names(NMD)[2:7] <- as.numeric(gsub("NMD_exp_", "", names(NMD)[2:7] ))

plotNMD <- melt(NMD, id.vars="Gene_name")
plotCon <- melt(Con, id.vars="Gene_name")
plotExpr <- melt(Expr, id.vars="Gene_name")

###############
# QC QC QC
file
cell_line
drug
design
head(NMD)
head(Con)
head(Expr)
####################################################
########## FOR checking slope only #################

design.mat <- cbind(1,1:6)
response.mat <- t(NMD[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
NMD <- cbind(NMD, t(reg))
names(NMD)[-(1:7)] <- c("NMD_Intercept","NMD_Slope")
NMD_ord<-NMD[order(-NMD$NMD_Slope),]

response.mat <- t(Con[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Con <- cbind(Con, t(reg))
names(Con)[-(1:7)] <- c("Con_Intercept","Con_Slope")
Con_ord<-Con[order(-Con$Con_Slope),]

response.mat <- t(Expr[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Expr <- cbind(Expr, t(reg))
names(Expr)[-(1:7)] <- c("Expr_Intercept","Expr_Slope")
Expr_ord<-Expr[order(-Expr$Expr_Slope),]


#############################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Added a title using one row at the bottom of all graphs
# title needs 'graphtitle' variable
#############################################

# variables needed
# might need to change the heights ratio
graphspace = 29
titlespace = 1
# title font size
fs=20

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    intlayout <- matrix(seq(from = 1, by = 1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
    layout <- rbind(intlayout,(numPlots+1))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights = c(rep_len(graphspace, nrow(layout)-1), titlespace) ) ))
    grid.text(graphtitle, gp = gpar(fontsize = fs), vp = viewport(layout.pos.row = ceiling(numPlots/cols)+1, layout.pos.col = 1:cols))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Separate plots
p1<-ggplot(plotNMD, aes(y=value,x=variable)) + 
  geom_point() + 
  stat_smooth() +
  xlab(bquote(.(series) ~ mu * 'M' ~ .(drug) ~ 'in' ~ .(cell_line) )) +
  ylab(paste("Normalized Expression of NMD isoform", cell_line, sep=" in ")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  facet_wrap(~Gene_name, scales = "free") +
  theme(legend.position="none", axis.text.x=element_blank(),
        strip.text.x = element_text(size = 6, colour = "red"))

p2<-ggplot(plotCon, aes(y=value,x=variable)) + 
  geom_point() + 
  stat_smooth() +
  xlab(bquote(.(series) ~ mu * 'M' ~ .(drug) ~ 'in' ~ .(cell_line) )) +
  ylab(paste("Normalized Expression of Consensus isoform", cell_line, sep=" in ")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  facet_wrap(~Gene_name, scales = "free") +
  theme(legend.position="none", axis.text.x=element_blank(),
         strip.text.x = element_text(size = 6, colour = "blue"))

p3<-ggplot(plotExpr, aes(y=value,x=variable)) + 
  geom_point() + 
  stat_smooth() +
  xlab(bquote(.(series) ~ mu * 'M' ~ .(drug) ~ 'in' ~ .(cell_line) )) +
  ylab(paste("Normalized Gene Expression", cell_line, sep=" in ")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  facet_wrap(~Gene_name, scales = "free") +
  theme(legend.position="none", axis.text.x=element_blank(),
        strip.text.x = element_text(size = 6, colour = "darkgreen"))

if (makepdf == TRUE) { pdf(pdfname,height=15, width=30) }

multiplot(p1, p2, p3, cols=3)

if (makepdf == TRUE) { dev.off() }


