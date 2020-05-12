# Script to analyse MS data
# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# install.packages("gplots")
library(gplots)

wd="/home/dyap/Projects/Takeda_T3/MS data/Expt2_ectopicCLK2_Tag-unTag_IP"
pat="*ch_CLK-IP-MS_TMT6_proteins.csv"
setwd(wd)
file_names = list.files(pattern = pat)

all_MS<-read.csv(file_names[1],header=TRUE)

cut=1
p=0.05
pep=2


N1CLK2<-subset(all_MS, (rN1 > cut & N_p.adj < p & pepNum > pep), select=c(Gene))
N2CLK2<-subset(all_MS, (rN2 > cut & N_p.adj < p & pepNum > pep), select=c(Gene))
C1CLK2<-subset(all_MS, (rC1 > cut & C_p.adj < p & pepNum > pep), select=c(Gene))
C2CLK2<-subset(all_MS, (rC2 > cut & C_p.adj < p & pepNum > pep), select=c(Gene))
CLK2<-subset(all_MS, (rU > cut & All_p.adj < p & pepNum > pep), select=c(Gene))

inputN<-list(N1CLK2,N2CLK2)
inputC<-list(C1CLK2,C2CLK2)
inputU<-list(CLK2)

tmpN <- venn(inputN, show.plot=FALSE)
Noverlap<-attr(tmpN, "intersections")$'11'
tmpC <- venn(inputC, show.plot=FALSE)
Coverlap<-attr(tmpC, "intersections")$'11'
tmpU <- venn(inputU, show.plot=FALSE)
Uoverlap<-attr(tmpU, "intersections")$'11'

input<-list(Noverlap,Coverlap,Uoverlap)

pdf(file="CLK2_allIP.pdf")
par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
#box("figure",lty="solid", col="green")

#plot.new()
#tmp <- venn(input, showSetLogicLabel=TRUE)
tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
# all4<-attr(tmp, "intersections")$'1111'
# The interestion of all 3 conditions
overlap<-attr(tmp, "intersections")$'111'

# This section prints out the genes in the intersection
        createCounter <- function(value) { function(j) { value <<- value+j} }
        createRowCounter <- function(value) { function(j) { value <<- value-r} }
        counter <- createCounter(1)
        rowc <- createRowCounter(30)
	name <- createCounter(1)

#        mtext("IP-MS ectopic expressed FLAG-tagged CLKs", side=3)
        mtext("IP-MS ectopic expressed untagged or FLAG-tagged CLK2 (n=6)", side=3)
#	text(40,350, "CLK1")
	text(200,380, "Untagged CLK2 IP (n=2)")
     	text(55,80, "N-FLAG tagged CLK2 IP (n=2)", srt=-45)
     	text(335,80, "C-FLAG tagged CLK2 IP (n=2)", srt=45)
#	text(365,350, "CLK3_cs")

dev.off()

write.table(overlap, file = "CLK2_allIP.txt", append = FALSE, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = FALSE )



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


