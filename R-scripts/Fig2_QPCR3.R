# start from scratch

path="/home/dyap/Projects/Takeda_T3/siRNA"
setwd(path)
file="151104_SRRM1_2_siKD_qPCR_Result_Data_Expt1"

qpcr<-read.table(file=file, sep="\t", header=TRUE)

df<-qpcr[c("Sample","Detector","Task","Ct")]

boxplot(qpcr$Ct ~ qpcr$Sample)
