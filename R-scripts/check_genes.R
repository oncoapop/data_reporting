# R-script to check all the synomyms of genes 
# then use them to check against the list of genes in a list

hgnc <- 
read.delim(url("http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&hgnc_dbtag=on&col=gd_app_sym&col=gd_aliases&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit&.cgifields=&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag") )

genes<-c("CDC16","CD27","CENPE")

for (ri in 1:length(genes)) {

input<-genes[ri]
syn<-hgnc[hgnc[,1] %in% genes[ri],c(1)]
sym<-hgnc[hgnc[,2] %in% genes[ri],c(2)]

input
syn
sym

print("****")

}
