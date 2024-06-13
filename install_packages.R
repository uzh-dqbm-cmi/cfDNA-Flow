# Title     : Install R package dependencies
# Objective : prepare environment to run the R scripts in this project

# Set the CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_dependencies <- function(packages) for (package in packages) if (!requireNamespace(package, quietly = TRUE)) install.packages(package, quietly = TRUE, dependencies = TRUE, force = TRUE) else print(paste("Package", package, " already installed."))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quietly = TRUE)
install_Bioc_dependencies <- function(packages) for (package in packages) if (!requireNamespace(package, quietly = TRUE)) BiocManager::install(package, quietly = TRUE) else print(paste("Bioconductor package", package, " already installed."))

install_jypiter_dependencies <- function (packages) {
  if (!requireNamespace("IRkernel", quietly = TRUE)){
    install.packages("IRkernel", quietly = TRUE)
  }
  IRkernel::installspec(user = FALSE)
  install_dependencies(packages)
}

#install.packages("xml2", quietly = TRUE, dependencies = TRUE)
#install.packages("stringi", quietly = TRUE, dependencies = TRUE)
#update.packages()
# install_dependencies(packages = c("xml2", "stringi"))

install_dependencies(packages = c("dplyr", "Seurat", "reshape", "data.table", "readr", "hash", "ggplot2","scales",
                                   "RColorBrewer", "gridExtra", "grid", "gtable", "plotly", "gridExtra","optparse",
                                   "devtools", "this.path", "shiny", "openssl", "plotmm", "mixtools", "EMCluster",
                                   "flexmix", "tidyverse", "icesTAF", "parallel",  "doParallel", "doMC", "doFuture", "plyr")) # , "Rdsm"
install_Bioc_dependencies(packages = c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene",
                                        "TxDb.Hsapiens.UCSC.hg38.knownGene", "QDNAseq.hg19", "QDNAseq.hg38", "BSgenome.Hsapiens.UCSC.hg38",
                                        "BSgenome.Hsapiens.UCSC.hg19"))

library(devtools)
devtools::install_github("asntech/QDNAseq.hg38@main")
#
# # CNAclinic; tMAD
#devtools::install_github("sdchandra/CNAclinic", build_vignettes = TRUE, dependencies=TRUE)
#library(CNAclinic)



# Install annotation packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "QDNAseq.hg19"))
# Installing CNAclinic
install.packages("devtools")
library(devtools)
install_github("sdchandra/CNAclinic", build_vignettes = TRUE, dependencies=TRUE)
library(CNAclinic)


# # ichorCNA
devtools::install_github("broadinstitute/ichorCNA")

## install_jypiter_dependencies(packages = c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'uuid', 'digest'))
