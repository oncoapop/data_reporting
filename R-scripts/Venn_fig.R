# Script to generate Venn diagram
# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# install.packages("gplots")
library(gplots)

# For MS data
msdir="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
ms=paste(msdir,"CLK2_allIP.txt",sep="/")

# For CG factors
cgdir="/home/dyap/Projects/Takeda_T3/CG_factors"
cg=paste(cgdir,"CG-factors",sep="/")

# For TiO2 phosphopeptide expt
tidir="/home/dyap/Projects/Takeda_T3/TiO2"
phos=paste(tidir,"Filtered_top42_splice.txt",sep="/")

msf<-read.table(file=ms)
cgf<-read.table(file=cg)
phosf<-read.table(file=phos)

input<-list(msf,cgf,phosf)

out="/home/dyap/Projects/Takeda_T3"
pdf(file=paste(out,"Venn_Intersect.pdf",sep="/"))
par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
#box("figure",lty="solid", col="green")

#plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
# all4<-attr(tmp, "intersections")$'1111'
# The interestion of all conditions
overlap<-attr(tmp, "intersections")

# This section prints out the genes in the intersection
        createCounter <- function(value) { function(j) { value <<- value+j} }
        createRowCounter <- function(value) { function(j) { value <<- value-r} }
        counter <- createCounter(1)
        rowc <- createRowCounter(30)
	name <- createCounter(1)

#        mtext("IP-MS ectopic expressed FLAG-tagged CLKs", side=3)
#        mtext("Venn diagram", side=3)
#	text(40,350, "CLK1")
	text(200,380, "Decreased phospho-peptide in TiO2")
     	text(55,80, "High Conf interactor of CLK2", srt=-45)
     	text(335,80, "Factors w/motif enrich in CG", srt=45)
#	text(365,350, "CLK3_cs")

dev.off()

outfile=paste(out,"Intersects.txt",sep="/")

# works but no header
#writeLines(unlist(lapply(overlap, paste, collapse=" ")))

sink(outfile) 
lapply(overlap, print) 
sink() 



########################## 
rows=ceiling((length(all4))/10)

	for (r in seq(rows)){
		print(r)
		y <- ((rows-r)*8)-20
              for (j in 1:10) {x <- counter(40)+r; text(x, y, all4[name(1)], srt=35, cex=0.6); print(x)}
	        counter <- createCounter(1)
	        name <- createCounter(r*10)
		print("hello world")
		print(y)
		}


