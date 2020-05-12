# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name in three (3) place in this script to make 
it run correctly

# R script to read in the calculated and obtained Tms of amplicons
# Change the sample name in three (3) place in this script to make 
it run correctly

library(calibrate)

#for (ri in "SA429-N" "SA429-X2")
#{
df<-read.csv(file="/home/dyap/dyap_temp/Temp_SA496-X2.txt", header=TRUE)

#df<-read.csv(file="/home/dyap/dyap_temp/Temp_"ri".txt", header=TRUE)

# This lists the number of cases with missing values (these are to be momitted)
df[!complete.cases(df),]

ff <- df[!complete.cases(df),]

# Include only complete cases with no NA values
ef <- na.omit(df)

pdf("/home/dyap/Projects/Tumour_Evol/Primer3_outputs/Temp_Corr_SA496-X2.pdf", width=6, height=6)

plot(ef$Actual_Temp, ef$Cal_Temp, xlab="Measured Amplicon melting Temperature", ylab="Calculated Amplicon Melting Temp", main="Correlation of Amplicon Melting Temperatures for SA496-X2", pch=19, ylim=c(60,100), xlim=c(60,100))

# text(seq(65,90, by=5), 65, label = ff$PCR_well, cex=0.8)

createCounter <- function(value) { function(i) { value <<- value+i} }
counter <- createCounter(60)

for (i in 1:length(ff$PCR_well)) {x <- counter(3); text(x, 65, ff$PCR_well[i], cex=0.9); print(x)}


textxy(ef$Actual_Temp, ef$Cal_Temp, ef$PCR_well, cx=0.8, m=c(mean(ef$Actual_Temp),mean(ef$Cal_Temp)))

dev.off()

#res=lm(ef$Cal_Temp ~ ef$Actual_Temp)

#abline(res)
#}

