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
#file1="HCT116_T_202_normalized.fpkm"
#file1="HCT116_T_595_normalized.fpkm"
#file1="Hela_T_202_normalized.fpkm"
file1="Hela_T_595_normalized.fpkm"
dir1="/home/amazloomian/Projects/SplicingProject/EIF4A3_STAR/expressionClustering_NMD_plots/data_normalized"
infile1=paste(dir1,file1,sep="/")

# To make sure do not mix and match, use $base
base=""
#file2=paste(base, "HCT116_T_228_normalized.fpkm", sep="")
#file2=paste(base,"HCT116_T_598_normalized.fpkm", sep="")
#file2=paste(base,"Hela_T_228_normalized.fpkm", sep="")
file2=paste(base,"Hela_T_598_normalized.fpkm", sep="")
dir2="/home/amazloomian/Projects/SplicingProject/EIF4A3_STAR/expressionClustering_NMD_plots/data_normalized"
infile2=paste(dir2,file2,sep="/")

##########################################################

# Outputs
cell_line="HeLa"
#cell_line="HCT-116"
drug1="T-595"
#drug1="T-202"
drug2="T-598"
#drug2="T-298"
series="0, 0.5, 2.0, 5.0, 10, 20"
Expt<-paste(paste(paste(design,cell_line,sep="_"),drug1,sep="_"),drug2,sep="-")
#Design
design=paste("Normalized_",drug1,"vs",drug2,sep="")


outdir="/home/dyap/Projects/eIF4A3_NMD/Plots"
fname=paste(outdir,Expt,sep="/")

##########################################################


data1<-read.table(file=infile1, header=TRUE)
data2<-read.table(file=infile2, header=TRUE)

#Expr1 <- data1[grep("MAPK", data1$gene_short_name), ]
Expr1<-data1[data1[,2] %in% c("MAPK4","MAPK2","MAPK1", "MAPK3","MAPK5","MAPK6","MAPK7","MAPK8","MAPK9","MAPK10","MAPK11","MAPK12","MAPK13","MAPK14","MAPK15"),]

#Expr2 <- data2[grep("MAPK", data2$gene_short_name), ]
Expr2<-data2[data2[,2] %in% c("MAPK4","MAPK2","MAPK1", "MAPK3","MAPK5","MAPK6","MAPK7","MAPK8","MAPK9","MAPK10","MAPK11","MAPK12","MAPK13","MAPK14","MAPK15"),]
plotExpr1 <- melt(Expr1, id.vars="gene_short_name")

plotExpr2 <- melt(Expr2, id.vars="gene_short_name")

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

response.mat <- t(Expr1[,3:8])
reg <- lm.fit(design.mat, response.mat)$coefficients
Expr1 <- cbind(Expr1, t(reg))
names(Expr1)[-(1:8)] <- c("Expr1_Intercept","Expr1_Slope")
Expr1_ord<-Expr1[order(-Expr1$Expr1_Slope),]

response.mat <- t(Expr2[,3:8])
reg <- lm.fit(design.mat, response.mat)$coefficients
Expr2 <- cbind(Expr2, t(reg))
names(Expr2)[-(1:8)] <- c("Expr2_Intercept","Expr2_Slope")
Expr2_ord<-Expr2[order(-Expr2$Expr2_Slope),]


###################################################
# testing ########## for all genes

conc=c("0.0","0.5","2.0","5.0","10.0","20.0")
names(Expr1)[3:8]<-conc
names(Expr2)[3:8]<-conc

		pdf3=paste(outdir, "/", "Chart_", "MAPKs_", drug1, "_", cell_line, ".pdf", sep="" )
		pdf(pdf3,height=15, width=10)
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		row=4
		coln=2
		par(mfrow=c(row,coln))

for (i in 1:length(Expr1$gene_short_name)) {

	gene<-(as.character(Expr1$gene_short_name[i]))
	y1<-as.numeric(Expr1[i,][3:8])
	y2<-as.numeric(Expr2[i,][3:8])
	x<-as.numeric(names(Expr1)[3:8])

	low=min(y1,y2)
	high=max(y1,y2)

		col1="red"
		col2="dark green"
		des=paste(drug1," (",col1, ") vs ", drug2," (",col2,")")
  
		if ( low == high )
			{
			print("skip!")
			} else {
				plot(x,y1,type="l",col=col1, ylim=range(low,high), xlab="", ylab="")
				par(new=TRUE)
				plot(x,y2,type="l",col=col2, ylim=range(low,high), xlab="", ylab="")
				title(main = gene, sub = des, xlab = "Drug conc", ylab = "Relative Norm. Gene Expr", cex=0.5)

		fit1 <- lm(y1 ~ x)
		fit2 <- lm(y2 ~ x)
#		abline(fit1, lty=1)
#		abline(fit2, lty=5)
		summary(fit1)
		summary(fit2)

		c1<-fit1$coefficient[1]
		m1<-fit1$coefficient[2]
		r1<-summary(fit1)$r.squared
		ypos1<-3/4 * ceiling(max(y1) * 2) / 2
		eq1<-paste("y =", round(m1,3), "x + ", round(c1,3), " r^2 = ", round(r1,3))

		c2<-fit2$coefficient[1]
		m2<-fit2$coefficient[2]
		r2<-summary(fit2)$r.squared
		ypos2<-3/4 * ceiling(max(y2) * 2) / 2
		eq2<-paste("y =", round(m2,3), "x + ", round(c2,3), " r^2 = ", round(r2,3))

		# Add Legend
		legend("topleft",legend=c(drug1,drug2),
  		text.col=c(col1,col2), lty=1:1, col=c(col1,col2), cex=0.8)
#		legend("topright",legend=c(eq1,eq2),
#  		text.col=c(col1,col2),lty=1:5, col=c(col1,col2), cex=0.8)
				}

	}

#%%%%%%%%%%%%%%%%%%%%%%%%
		dev.off() 



