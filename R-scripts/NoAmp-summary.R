# This script takes that output in summary *.read.txt.sum.txt
# Using the command line command as follows:
# grep "0,0,0,0,0" *.sum.txt | awk -F"-" '{print $2,$3}' | awk -F: '{print $1","$2","$3}' | awk -F, '{print $1,$2,$3}' | awk -F" " '{print $1,$3":"$4}' > noamp.txt
# This give the positions and nuclei in which the primers give no reads

noamp <- 
read.table(file="/share/lustre/backup/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq/Blastfasta/summary/noamp.txt", 
stringsAsFactors = FALSE)

names(noamp)[1] <- "Nuclei#"
names(noamp)[2] <- "Position"

count <- table(noamp$"Position")

sorted <- sort(count)

pdf("/share/lustre/backup/dyap/Projects/Single_Cell/Karn-VBA0038-run3-MiSeq/Blastfasta/summary/noamp.pdf", width=15, height=20)

par(mar=c(8,4,4,1)+.1)

plot(sorted, main="Positions with no reads in Karn's dataset (28 nuclei)", ylab="No of Nuclei", las=2, xaxt = "n")

midp<-seq(1,21,by=1)
axis(1, at=c(midp), labels=names(sorted), col.axis="black", las=2, cex.lab=0.5)

dev.off()


