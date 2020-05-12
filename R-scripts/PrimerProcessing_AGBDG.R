# This R script is modified from Steve Mckinney's script to take into
# Miseq runs which have been normalized with ACTB, GAPDH or TUB

projdir <- "/home/dyap/Projects/Takeda_T3/CG"
wdir <- "/home/dyap/Projects/Takeda_T3/CG_factors/screen"

# From Miseq run /share/lustre/archive/MiSeq/MiSeq_Analysis_Files/160209_M02348_0108_000000000-AGBDG
infile="KD_CG_factors_readthrough_filtered_ctrl_norm.tsv"
miseq="AGBDG"
outfile=paste("MiSeqValidation_CG_Panel",miseq,sep="-")

input <- paste(wdir, infile, sep = "/")
output <- paste(wdir, outfile, sep = "/")

pri <- paste(projdir, "/CG_primers_ordered.txt", sep = "")

pridf <- read.table(file = pri, header = FALSE, stringsAsFactors = FALSE)
names(pridf)[1] <- "Primers"
head(pridf)

ppdf <- read.table(file = input, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

head(ppdf)

# This gets the samples from the data itself
samples <- unique(ppdf$sample_id)
sdf <- data.frame(Samples = samples, stringsAsFactors = FALSE)

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
presdf[presdfpso, c("sample_id", "Primers",  "ACTB_norm", "GAPDH_norm", "TUBA1B_norm")]
head(presdf[presdfpso, c("sample_id", "Primers", "ACTB_norm", "GAPDH_norm", "TUBA1B_norm")])
     
presdf[presdfpso, c("sample_id",  "ACTB_norm", "GAPDH_norm", "TUBA1B_norm")]

table(apply(presdf[, c("ACTB_norm", "GAPDH_norm", "TUBA1B_norm")], 1, function(x) sum(is.na(x))))

### All are missing, or none are missing.

X11()

require("nlme")

presdf$sidf <- factor(presdf$sample_id)
presdf$sidn <- as.numeric(presdf$sidf)
presdf$Primerssp <- gsub(":", "\n", presdf$Primers)
trgd <- groupedData(ACTB_norm ~ sidn | Primerssp, data = presdf, order.groups = FALSE)
pdf(file = paste(output,"ACTB_v01.pdf",sep="_"), width = 8, height = 10)
plot(trgd, aspect = "fill", par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table = TRUE)
dev.off()

trgd <- groupedData(GAPDH_norm ~ sidn | Primerssp, data = presdf, order.groups = FALSE)
pdf(file =  paste(output, "GAPDH_v01.pdf", sep="_"), width = 8, height = 10)
plot(trgd, aspect = "fill", par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table = TRUE)
dev.off()

trgd <- groupedData(TUBA1B_norm ~ sidn | Primerssp, data = presdf, order.groups = FALSE)
pdf(file =  paste(output, "TUBA1B_v01.pdf", sep="_"), width = 8, height = 10)
plot(trgd, aspect = "fill", par.strip.text=list(cex=0.7, lines = 3), layout = c(4, 5), as.table = TRUE)
dev.off()
