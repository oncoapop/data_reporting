# Script to plot data from NMD expression correlation

# Load R libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
library(scales)

##########################################################
### CHANGE THESE PARAMETERS FOR EACH RUN OF THE SCRIPT ###
##########################################################
 
# Inputs
#file1="HCT116_T_202_normalized_highest_expression.res"
file1="HCT116_T_595_normalized_highest_expression.res"
#file1="HCT116_T_202_normalized_longest.res"
#file1="HCT116_T_595_normalized_longest.res"
#file1="Hela_T_202_normalized_longest.res"
#file1="Hela_T_595_normalized_longest.res"
dir1="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation"
infile1=paste(dir1,file1,sep="/")

# To make sure do not mix and match, use $base
base=strsplit(file1,split="[.]")[[1]][1]
#file2=paste(base,"HCT116_T_228_normalized.res", sep="_")
file2=paste(base,"HCT116_T_598_normalized.res", sep="_")
#file2=paste(base,"HCT116_T_228_normalized.res", sep="_")
#file2=paste(base,"HCT116_T_598_normalized.res", sep="_")
#file2=paste(base,"Hela_T_228_normalized.res", sep="_")
#file2=paste(base,"Hela_T_598_normalized.res", sep="_")
dir2="/share/scratch/amazloomian_temp/EIF4A3_STAR/validation/other_compounds"
infile2=paste(dir2,file2,sep="/")

##########################################################

# Outputs
design="Normalized_NMD_Expr_vs_Highest_Isoform"
# design="Normalized_NMD_Expr_vs_Longest_Isoform"
# cell_line="HeLa"
cell_line="HCT-116"
drug1="T-595"
drug2="T-598"
series="0, 0.5, 2.0, 5.0, 10, 20"
Expt<-paste(paste(paste(design,cell_line,sep="_"),drug1,sep="_"),drug2,sep="-")

outdir="/home/dyap/Projects/eIF4A3/Plots"
fname=paste(outdir,Expt,sep="/")
pdfname=paste(fname,"pdf",sep=".")
graphtitle=gsub("_"," ",design)
makepdf = TRUE

##########################################################


data1<-read.table(file=infile1, header=TRUE)
data2<-read.table(file=infile2, header=TRUE)


NMD1 <- data1[c(1,6:11)]
Con1 <- data1[c(1,12:17)]
Expr1 <- data1[c(1,18:23)]

NMD2 <- data2[c(1,6:11)]
Con2 <- data2[c(1,12:17)]
Expr2 <- data2[c(1,18:23)]

names(NMD1)[2:7] <- as.numeric(gsub("NMD_exp_", "", names(NMD1)[2:7] ))
names(NMD2)[2:7] <- as.numeric(gsub("NMD_exp_", "", names(NMD2)[2:7] ))

plotNMD1 <- melt(NMD1, id.vars="Gene_name")
plotCon1 <- melt(Con1, id.vars="Gene_name")
plotExpr1 <- melt(Expr1, id.vars="Gene_name")

plotNMD2 <- melt(NMD2, id.vars="Gene_name")
plotCon2 <- melt(Con2, id.vars="Gene_name")
plotExpr2 <- melt(Expr2, id.vars="Gene_name")

###############
# QC QC QC

file1
file2
cell_line
drug1
drug2
design

####################################################
########## FOR checking slope only #################

design.mat <- cbind(1,1:6)

response.mat <- t(NMD1[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
NMD1 <- cbind(NMD1, t(reg))
names(NMD1)[-(1:7)] <- c("NMD1_Intercept","NMD1_Slope")
NMD1_ord<-NMD1[order(-NMD1$NMD1_Slope),]

response.mat <- t(Con1[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Con1 <- cbind(Con1, t(reg))
names(Con1)[-(1:7)] <- c("Con1_Intercept","Con1_Slope")
Con1_ord<-Con1[order(-Con1$Con1_Slope),]

response.mat <- t(Expr1[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Expr1 <- cbind(Expr1, t(reg))
names(Expr1)[-(1:7)] <- c("Expr1_Intercept","Expr1_Slope")
Expr1_ord<-Expr1[order(-Expr1$Expr1_Slope),]

response.mat <- t(NMD2[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
NMD2 <- cbind(NMD2, t(reg))
names(NMD2)[-(1:7)] <- c("NMD2_Intercept","NMD2_Slope")
NMD2_ord<-NMD2[order(-NMD2$NMD2_Slope),]

response.mat <- t(Con2[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Con2 <- cbind(Con2, t(reg))
names(Con2)[-(1:7)] <- c("Con2_Intercept","Con2_Slope")
Con2_ord<-Con2[order(-Con2$Con2_Slope),]

response.mat <- t(Expr2[,2:7])
reg <- lm.fit(design.mat, response.mat)$coefficients
Expr2 <- cbind(Expr2, t(reg))
names(Expr2)[-(1:7)] <- c("Expr2_Intercept","Expr2_Slope")
Expr2_ord<-Expr2[order(-Expr2$Expr2_Slope),]


###################################################
# testing ########## for all genes

trans<-NULL
as.data.frame(trans)

################# run this for manual selection step ##########

for (i in 1:length(NMD1$Gene_name)) {

	gene<-(as.character(NMD1$Gene_name[i]))
	y1<-as.numeric(NMD1[i,][2:7])
	y2<-as.numeric(NMD2[i,][2:7])
	x<-as.numeric(names(NMD1)[2:7])

	plot(x,y1,type="l",col="red")
	lines(x,y2,col="green")

	fit1 <- lm(y1 ~ x)
	fit2 <- lm(y2 ~ x)
	abline(fit1)
	abline(fit2)
	summary(fit1)
	summary(fit2)

	c1<-fit1$coefficient[1]
	m1<-fit1$coefficient[2]
	r1<-summary(fit1)$r.squared
	ypos1<-3/4 * ceiling(max(y1) * 2) / 2
	ypos21<-1/4 * ceiling(max(y1) * 2) / 2
	eq1<-paste("y =", round(m1,3), "x + ", round(c1,3), " r^2 = ", round(r1,3))
	text(2,ypos1, eq1)

	c2<-fit2$coefficient[1]
	m2<-fit2$coefficient[2]
	r2<-summary(fit2)$r.squared
	ypos2<-3/4 * ceiling(max(y2) * 2) / 2
	ypos22<-1/4 * ceiling(max(y2) * 2) / 2
	eq2<-paste("y =", round(m2,3), "x + ", round(c2,3), " r^2 = ", round(r2,3))
	text(2,ypos2, eq2)

	title(main = gene, sub = NULL, xlab = "Drug conc", ylab = "Relative Norm. Con Gene Expression",
        line = NA, outer = FALSE)

	# manual selection processed based on the viewing of graphs
	#save<-readline(prompt="Press [s]+[ENTER] to save or [ENTER] alone to skip - Enter your response: ")

	#if (save == "s") 
	#	{ 
		pdf3=paste(outdir, "/", "Graphs_", gene, "-", i, ".pdf", sep="" )
		pdf(pdf3,height=10, width=7) 
		plot(x,y1,type="l",col="red")
		lines(x,y2,col="green")
		abline(fit1)
		abline(fit2)
		text(2,ypos1, eq1)
		text(2,ypos2, eq2)
		title(main = gene, sub = "Dose Response : T-595 (red) and T-598 (green)", xlab = "Drug conc", ylab = "Relative Norm. Con Gene Expression")
		dev.off() 

		trans$Gene_name[i]<-as.character(data1[i,]$Gene_name)
		trans$NMD_transcript_id[i]<-as.character(data1[i,]$NMD_transcript_id)
		trans$consensus_transcript[i]<-as.character(data1[i,]$consensus_transcript)
	#	}
	}


################# run this for manual selection step ##########

trans<-as.data.frame(trans)
seltrans<-trans[complete.cases(trans),]

outfile=paste(outdir, "all_transcripts.csv", sep="/" )
write.csv(seltrans, file = outfile)

############################################
# Expression of Consensus (highest expressing isoform)

Trend <- NULL

Trend$Gene <-  Con1[,c("Gene_name")]
Trend$T.595 <-  Con1[,c("Con1_Slope")]
Trend$T.598 <-  Con2[,c("Con2_Slope")]
Trend$T.595.intercept <-  Con1[,c("Con1_Intercept")]
Trend$T.598.intercept <-  Con2[,c("Con2_Intercept")]
df<-as.data.frame(Trend)

# Plot histogram
pdf1="/home/dyap/Projects/eIF4A3/Plots/Con_Expr_histo.pdf"

if (makepdf == TRUE) { pdf(pdf1,height=15, width=10) }

bins <- seq(-20, 65, by=0.4)
hist(df$T.595, breaks=bins, col=rgb(1,0,0,0.5), main="T-595 causes greater variation in Norm. Con. Expr.", xlab="Effect of Drug treatment, slope=(drug conc/norm consensus expr)", ylim=c(0,90))
hist(df$T.598, breaks=bins, col=rgb(0,0,1,0.5), add=T)
legend("topright", c("T-595", "T-598"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

if (makepdf == TRUE) { dev.off() }

############################################
# Expression of NMD isoform

Trend_NMD <- NULL

Trend_NMD$Gene <-  NMD1[,c("Gene_name")]
Trend_NMD$T.595 <-  NMD1[,c("NMD1_Slope")]
Trend_NMD$T.598 <-  NMD2[,c("NMD2_Slope")]
Trend_NMD$T.595.intercept <-  NMD1[,c("NMD1_Intercept")]
Trend_NMD$T.598.intercept <-  NMD2[,c("NMD2_Intercept")]
df<-as.data.frame(Trend_NMD)

# Plot histogram
pdf2="/home/dyap/Projects/eIF4A3/Plots/NMD_Expr_histo.pdf"

if (makepdf == TRUE) { pdf(pdf2,height=15, width=10) }

bins <- seq(-20, 65, by=1)
hist(df$T.595, breaks=bins, col=rgb(1,0,0,0.5), main="T-595 causes greater increase in NMD isoform Expr.", xlab="Effect of Drug treatment, slope=(drug conc/norm NMD expr)", ylim=c(0,150))
hist(df$T.598, breaks=bins, col=rgb(0,0,1,0.5), add=T)
legend("topright", c("T-595", "T-598"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

if (makepdf == TRUE) { dev.off() }


# Order the original dataframe by decreasing NMD1_Slope (effect of T-595)
# Top increase in expression of NMD isoforms by T-595 (data1)
data1_ord<-data1[order(match(data1$Gene_name,NMD1_ord$Gene_name)),]

NMD1 <- data1_ord[c(1,6:11)]
Con1 <- data1_ord[c(1,12:17)]
Expr1 <- data1_ord[c(1,18:23)]

names(NMD1)[2:7] <- as.numeric(gsub("NMD_exp_", "", names(NMD1)[2:7] ))

plotNMD <- melt(NMD1, id.vars="Gene_name")
plotCon <- melt(Con1, id.vars="Gene_name")
plotExpr <- melt(Expr1, id.vars="Gene_name")

drug<-drug1
pdfname1<-"/home/dyap/Projects/eIF4A3/Plots/Re-plot_by_NMD_slope.pdf"
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


fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = FALSE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     parse(text=l)
}

# 
plotNMD$value <- factor(plotNMD$value, levels = sort(unique(plotNMD$value)))


# Separate plots
p1<-ggplot(plotNMD, aes(y=value,x=variable)) + 
  geom_point() + 
  stat_smooth() +
  xlab(bquote(.(series) ~ mu * 'M' ~ .(drug) ~ 'in' ~ .(cell_line) )) +
  ylab(paste("Normalized Expression of NMD isoform", cell_line, sep=" in ")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  facet_wrap(~Gene_name, scales = "free") +
  scale_y_discrete(breaks= pretty_breaks()) +
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

if (makepdf == TRUE) { pdf(pdfname1,height=15, width=30) }

multiplot(p1, p2, p3, cols=3)

if (makepdf == TRUE) { dev.off() }


