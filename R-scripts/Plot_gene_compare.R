# Script to plot data from NMD expression correlation

# Load R libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
library(scales)

##########################################################
### CHANGE THESE PARAMETERS FOR EACH RUN OF THE SCRIPT ###
##########################################################

# Common variables
dir="/home/amazloomian/Projects/SplicingProject/EIF4A3_STAR/expressionClustering_NMD_plots/data_normalized"
outdir="/home/dyap/Projects/eIF4A3_NMD/Plots"

#Comparisons
cl<-c("HCT116","HCT116","Hela","Hela")
d1<-c("T-595","T-202")
d2<-c("T-598","T-598")
com<-data.frame(cl,d1,d2)

#Set-up
	col1="red"
	col2="dark green"
	row=2
	coln=2

#name<-"combined"
#pdf3=paste(outdir, "/", "Chart_",name,".pdf", sep="" )

gene=c("MAPK13")
#gene=c("AURKB","SRSF6","RBMX","SRSF2","RBM38","CDK7","CDC27","AURKB","BCL2L1","TIA1","MAPK6","MAPK7","MAPK13","G3BP1")
pdf3=paste(outdir, "/", "Chart_",gene,".pdf", sep="" )

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pdf(pdf3,height=10, width=7)
	par(mfrow=c(row,coln))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (j in 1:4 )
	{
	print(j)
		
cell_line<-as.character(com$cl[j])
drug1<-as.character(com$d1[j])
drug2<-as.character(com$d2[j])

fd1<-gsub("-","_", drug1)
fd2<-gsub("-","_", drug2)

# Inputs
infile1=paste(dir,"/",cell_line,"_",fd1,"_normalized.fpkm",sep="")
infile2=paste(dir,"/",cell_line,"_",fd2,"_normalized.fpkm",sep="")

# Outputs
design=paste(cell_line,"_",drug1,"vs",drug2,sep="")
print(design)

data1<-read.table(file=infile1, header=TRUE)
data2<-read.table(file=infile2, header=TRUE)

Expr1<-data1[data1[,2] %in% gene,]
Expr2<-data2[data2[,2] %in% gene,]

conc=c("0.0","0.5","2.0","5.0","10.0","20.0")
names(Expr1)[3:8]<-conc
names(Expr2)[3:8]<-conc

for (i in 1:length(Expr1$gene_short_name)) {

	trans<-(as.character(Expr1$tracking_id[i]))
	gene<-(as.character(Expr1$gene_short_name[i]))
	y1<-as.numeric(Expr1[i,][3:8])
	y2<-as.numeric(Expr2[i,][3:8])
	x<-as.numeric(names(Expr1)[3:8])

	low=min(y1,y2)
	high=max(y1,y2)

		des=paste("Dose Response of",trans)
		ylab=paste("Relative Normalized Expr (",gene,")",sep="")
  
		if ( low == high )
			{
				plot(x,y1,type="l",col=col1, xlab="", ylab="")
				par(new=TRUE)
				plot(x,y2,type="l",col=col2, xlab="", ylab="")
			} else {
				plot(x,y1,type="l",col=col1, ylim=range(low,high), xlab="", ylab="")
				par(new=TRUE)
				plot(x,y2,type="l",col=col2, ylim=range(low,high), xlab="", ylab="")
				title(main = cell_line, sub = des , xlab = "Drug conc (uM)", ylab = ylab, cex=1.0, cex.lab=1.0)
				legend("topleft",legend=c(drug1,drug2),
  				text.col=c(col1,col2), lty=1:1, col=c(col1,col2), cex=1.2)
				}
	}
}

#%%%%%%%%%%%%%%%%%%%%%%%%
		dev.off() 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
