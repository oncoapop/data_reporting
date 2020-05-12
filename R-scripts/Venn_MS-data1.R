	y<-counter(5)
# Script to analyse MS data
# myR
# R version 3.0.2 Patched (2014-01-16 r64804) -- "Frisbee Sailing"
# install.packages("gplots")
library(gplots)

wd="/home/dyap/Projects/Takeda_T3/MS data/Expt1_ectopicCLK2"
pat="*-genes"
setwd(wd)
file_names = list.files(pattern = pat)

# Reads the file name into Q (in a list)
QC_list = lapply(file_names, read.table, header = FALSE)


for (i in seq(length(QC_list)))	{
	assign(strsplit(file_names[i],split="-")[[1]][1], QC_list[[i]])
	print(strsplit(file_names[i],split="-")[[1]][1])	
	}

#input<-list(CLK1_sol,CLK2_sol,CLK3_cs,CLK3_sc)
input<-list(CLK1_sol,CLK2_sol)
#input<-list(CLK3_cs,CLK3_sc)

pdf(file="CLK1_2.pdf")
par(oma=c(0,1,2,1))
par(mar=c(6, 4, 2, 2) + 0.1)
par(pin=c(5,5))
#box("figure",lty="solid", col="green")

#plot.new()
tmp <- venn(input, showSetLogicLabel=TRUE)
#tmp <- venn(input, show.plot=TRUE)

# The interestion of all 4 conditions
# all4<-attr(tmp, "intersections")$'1111'
# The interestion of all 2 conditions
overlap<-attr(tmp, "intersections")$'11'

# This section prints out the genes in the intersection
        createCounter <- function(value) { function(j) { value <<- value+j} }
        createRowCounter <- function(value) { function(j) { value <<- value-r} }
        counter <- createCounter(1)
        rowc <- createRowCounter(30)
	name <- createCounter(1)
	y <- rowc(30)

#        mtext("IP-MS ectopic expressed FLAG-tagged CLKs", side=3)
        mtext("IP-MS ectopic expressed FLAG-tagged CLK1 & 2", side=3)
#	text(40,350, "CLK1")
#	text(137,375, "CLK2")
#	text(260,375, "CLK3_sc")
#	text(365,350, "CLK3_cs")

dev.off()

write.table(overlap, file = "CLK1_2.txt", append = FALSE, quote = FALSE, sep = ",",
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


