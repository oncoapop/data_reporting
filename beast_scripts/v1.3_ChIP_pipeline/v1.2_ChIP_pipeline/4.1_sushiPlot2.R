#source("http://bioconductor.org/biocLite.R")
#biocLite("Sushi")

# Modified for the original script from Sam Aparicio (Jun 2015)
# Modified by Damian Yap 9 Jun 2015

# Script to plot multiple ChIP-seq on the same graph and region for comparison
# Regions are specified in the data.frame bglist (this corresponds to the .broadpeak file
# usually output from the MACS)

# Comparisons 1. Drug-input-signal - H2AX-input-10-6CX_219_treat_pileup.bedgraph
#	      2. Drug-input-background - H2AX-input-10-6CX_219_control_lambda.bedgraph
#	      3. Drug-NODRUG-signal - H2AX-0vs10-6CX_219_treat_pileup.bedgraph
#             4. Drug-NODRUG-background - H2AX-0vs10-6CX_219_control_lambda.bedgraph

# Obtaining the list of relevant files in the directory

# This checks to see all the drug conc in the directory
drugconc<-system("ls | grep CX | grep -v bam | awk -F- '{print $4}' | awk -F_ '{print $1}' | grep CX | sort -u", intern=TRUE)

makepdf = TRUE

if (makepdf == TRUE)
{
   system("rm -f *.pdf")
}

for (m in 1:length(drugconc)) {

ctrltreat<-paste(paste("ls",paste("grep",drugconc[m],sep=" "),sep="|"),"| grep input | grep treat_pileup.bedgraph")
ctrlcontrol<-paste(paste("ls",paste("grep",drugconc[m],sep=" "),sep="|"),"| grep input | grep control_lambda.bedgraph")

expttreat<-paste(paste("ls",paste("grep",drugconc[m],sep=" "),sep="|"),"| grep 0vs10 | grep treat_pileup.bedgraph")
exptcontrol<-paste(paste("ls",paste("grep",drugconc[m],sep=" "),sep="|"),"| grep 0vs10 | grep control_lambda.bedgraph")

exptpeak<-paste(paste("ls",paste("grep",drugconc[m],sep=" "),sep="|"),"| grep 0vs10 | grep peaks.broadPeak")

# This is the list of the files to read in
Drug_input_treat<-system(ctrltreat, intern=TRUE)
Drug_input_control<-system(ctrlcontrol, intern=TRUE)
Drug_ChIP_nodrug_treat<-system(expttreat, intern=TRUE)
Drug_ChIP_nodrug_control<-system(exptcontrol, intern=TRUE)

peaks<-system(exptpeak, intern=TRUE)
fname<-peaks

########
# CODE
########

library("Sushi")
pdfname=paste(fname,"pdf",sep=".")

if (makepdf == TRUE)
{
   pdf(pdfname,height=10, width=12)
}

#FUNCTION FOR READING IN BEDGRAPHS VERY QUICKLY
read.bedgraph <- function(file) {
 dat <- scan(file=file, 
         what=list(character(),integer(),integer(),numeric()), 
         sep="\t", skip=1)
 dat <- data.frame(chr=dat[[1]], start=dat[[2]], 
         end=dat[[3]], val=dat[[4]])
 return(dat)
}


#READ IN THE BEDGRAPH FILES

Drug_input_Treat<-read.bedgraph(Drug_input_treat)
Drug_input_Control<-read.bedgraph(Drug_input_control)

Drug_ChIP_nodrug_Treat<-read.bedgraph(Drug_ChIP_nodrug_treat)
Drug_ChIP_nodrug_Control<-read.bedgraph(Drug_ChIP_nodrug_control)

#READ IN REGIONS - 
bedpeaks<-data.frame(read.table(peaks,col.names=c("chrom","chromstart","chromend","name","score","strand","value","pval","qval")))


#GENERIC FUNCTION FOR PLOTTING FOUR LIBRARIES AT A TIME
bed.region.plot <- function(d1,d2,d3,d4,chrom,chromstart,chromend,range,name,delta) 
{

# Place holder graph in event of failure to plot previous plots
# plot(1, type="n", axes=F, xlab="", ylab="")
plot.new()

try(plotBedgraph(d1,chrom,chromstart-delta,chromend+delta,transparency=.25,color=SushiColors(4)(4)[1],range=range,main=name))
try(plotBedgraph(d2,chrom,chromstart-delta,chromend+delta,transparency=.25,color=SushiColors(4)(4)[2],overlay=TRUE))
try(plotBedgraph(d3,chrom,chromstart-delta,chromend+delta,transparency=.25,color=SushiColors(4)(4)[3],overlay=TRUE))
try(plotBedgraph(d4,chrom,chromstart-delta,chromend+delta,transparency=.25,color=SushiColors(4)(4)[4],overlay=TRUE))


try(axis(side=2,las=2,tcl=.2))
try(mtext("Read Depth",side=2,line=1.75,cex=.75,font=2))
try(labelgenome(chrom,chromstart-delta,chromend+delta,n=10,scale="bp"))
try(legend("topright",inset=0.025,legend=c("Drug_input_Treat","Drug_input_Control","Drug_ChIP_nodrug_Treat","Drug_ChIP_nodrug_Control"),fill=opaque(SushiColors(4)(4)),border=SushiColors(4)(4),text.font=2,cex=1.0))

}

#FUNCTION FOR EXECUTING PLOTS OVER SEVERAL REGIONS
execute.plot<-function(bedpeaks) {
for (n in 1:nrow(bedpeaks)) {
bed.region.plot(Drug_input_Treat,Drug_input_Control,Drug_ChIP_nodrug_Treat,Drug_ChIP_nodrug_Control,bedpeaks$chrom[n],bedpeaks$chromstart[n],bedpeaks$chromend[n],c(0,100),bedpeaks$name[n],1500)
}
}

execute.plot(bedpeaks)


if (makepdf == TRUE)
{
	dev.off()
}

}

