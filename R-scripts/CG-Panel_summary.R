# R Script to summaries the results of the samples for each CG panel primer pair

projdir <- "/home/dyap/Projects/Takeda_T3/CG"

# This is the input file which can be in ang format
# This particular script takes in the SDS2.4 RQ format and recommends that you manually edit out blank lines and headers
# USE /home/dyap/Projects/Takeda_T3/CG/pre-process_SDS.sh to pre-process the SDS file

exptname="qPCR_QC_T3_treatment_CG_Panel"

infile="T3_CGPanel.csv"
sam="T3_conc_samples"
input <- paste(projdir,infile, sep = "/")
samfile <-paste(projdir, sam, sep="/")

pdffile <- paste(paste(prodir,exptname,sep="/"),"pdf",sep=".")

# This is the whole set of primers ordered
pri <- paste(projdir, "/CG_primers_ordered.txt", sep = "")

pridf <- read.table(file = pri, header = FALSE, stringsAsFactors = FALSE)
names(pridf)[1] <- "Primers"
head(pridf)

ppdf <- read.table(file = input, sep="\t", stringsAsFactors = FALSE, header = TRUE)

head(ppdf)

#####################
# samples <- read.table(file=samfile, sep="\n", stringsAsFactors = FALSE, header = FALSE)
# importing using read tables does NOT work
samples <- c("1_Old_0", "2_Old_Med", "3_Old_Hi", "4_New_0", "5_New_Med", "6_New_Hi")
sdf <- data.frame(Samples = samples, stringsAsFactors = FALSE)
#colnames(sdf)[1]<-"Samples"

### Ordering of primer gene names in output is arbitrary - need to make Primers variable same for either order of gene names.
ppdf$g1nm <- paste(ppdf$gene1, "@", ppdf$breakpt1, sep = "")
ppdf$g2nm <- paste(ppdf$gene2, "@", ppdf$breakpt2, sep = "")

g1priml <- lapply( ppdf$g1nm, function(x) which( grepl(x, pridf$Primers) ) )
g2priml <- lapply( ppdf$g2nm, function(x) which( grepl(x, pridf$Primers) ) )
intprim <- sapply( seq(nrow(ppdf) ), function(x) c(intersect(g1priml[[x]], g2priml[[x]] ), NA_integer_)[1] )
if (length(intprim) == nrow(ppdf)) { ppdf$PrimNo <- intprim } else { stop("Primer match error") }
ppdf$PrimersN <- sapply( seq(nrow(ppdf) ), function(x) c(intersect(g1priml[[x]], g2priml[[x]] ), NA_integer_)[1] )
ppdf$Primers <-  pridf$Primers[ppdf$PrimersN]

ppdf$g1prim <- sapply(ppdf$g1nm, function(x) which(grepl(x, pridf$Primers)))
ppdf$g2prim <- sapply(ppdf$g2nm, function(x) which(grepl(x, pridf$Primers)))

### ?? Need to make dataframe for all samples and all primer pairs and merge with pipeline data

prsadf <- expand.grid(as.character(samples), as.character(pridf$Primers), stringsAsFactors = FALSE)
names(prsadf) <- c("sample_id", "Primers")

ppdf$g1nm <- NULL
ppdf$g2nm <- NULL
ppdf$g1prim <- NULL
ppdf$g2prim <- NULL
ppdf$PrimersN <- NULL

presdf <- merge(prsadf, ppdf, all.x = TRUE, all.y = FALSE, stringsAsFactors = FALSE)
table(prsadf$sample_id, useNA = "always")
table(prsadf$Primers, useNA = "always")
table(table(prsadf$sample_id, useNA = "always"))
table(table(prsadf$Primers, useNA = "always"))
table(presdf$sample_id, useNA = "always")
table(presdf$Primers, useNA = "always")
table(table(presdf$sample_id, useNA = "always"))
table(table(presdf$Primers, useNA = "always"))
presdfpso <- order(presdf$Primers, presdf$sample_id)
presdf[presdfpso, ]
head(presdf[presdfpso, ])
presdf[presdfpso, c("sample_id", "Primers",  "RQ")]
head(presdf[presdfpso, c("sample_id", "RQ")])
     
presdf[presdfpso, c("sample_id",  "RQ")]


## Why two records here?
## SA467	20	-	TP53RK	SLC13A3	45317771	45242364	0.3893030794165316	1.0049608511146972	0.522612128182016
## SA467	20	-	TP53RK	SLC13A3	45317771	45242364	0.38981245658717295	1.0062757755065448	0.5232959313710254

table(apply(presdf[, c("RQ")], 1, function(x) sum(is.na(x))))
###   0   3 
### 417 367 

### All are missing, or none are missing.

X11()

require("nlme")

presdf$sidf <- factor(presdf$sample_id)
presdf$sidn <- as.numeric(presdf$sidf)
presdf$Primerssp <- gsub(":", "\n", presdf$Primers)
trgd <- groupedData(RQ ~ sidn | Primerssp, data = presdf, order.groups = FALSE)

pdf(file = pdffile, width = 8, height = 10)
plot(trgd, aspect = "fill", par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table = TRUE)
dev.off()

