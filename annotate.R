if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# BiocManager::available()
pkgs <- rownames(installed.packages())
BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)

source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("GenomicFeatures", "AnnotationDbi"))
BiocInstaller::biocLite(c("BSgenome","BSgenome.Hsapiens.UCSC.hg19"))

