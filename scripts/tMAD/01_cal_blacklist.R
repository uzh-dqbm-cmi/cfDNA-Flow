####################################################################
####### Create empirical blacklists from all control samples #######
####################################################################

library(CNAclinic)
library(QDNAseq.hg38)
library(this.path)

source(paste0(trimws(this.dir()), "/00_bedprocessing.R"))

###########################################
# Arguments to change
###########################################


args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=4) {
  ctrl_sample_dir <- args[1]   # input directory with BED files
  ctrl_sample_name <- readLines(args[2])  # List of healthy samples
  binSizes <- c(as.integer(args[3]))  # Reads binned into 1000 Kbp windows
  RCNormalization <- as.logical(args[4])
} else {
  stop("Unexpected number of arguments. usage: \n Rscript --vanilla 01_cal_blacklist.R <ctrl_sample_dir> <ctrl_sample_names.txt> <binZize> <RCNormalization:true,false>",
       call.=FALSE)
}

# binSizes <- c(1000)

# ctrl_sample_dir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# input directory with BED files
#ctrl_sample_dir <- "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"
# List of healthy samples
#ctrl_sample_name <- c("HV01","HV03","HV04","HV05","HV06","HV07")
# intput from
FragmentFileSuffix <- "_point.csv"

# RCNormalization <- TRUE
methodOfRCNormalization <- "mean"
path_to_blacklist <- ctrl_sample_dir

###########################################
###########################################

# For each bin size value calculate the corresponding blacklist:

for(b in 1:length(binSizes)){
  
  binSize <- binSizes[b]
  
  FragmentFiles <- paste0(ctrl_sample_dir,ctrl_sample_name,
                          FragmentFileSuffix)
  
  cat(paste0("Create empirical blacklists from samples:\n",paste(FragmentFiles, collapse = '\n'),'\n'))
  
  FragmentFileNames <- unlist(lapply(strsplit(FragmentFiles, "/"), function(x){ x[[length(x)]]}))
  
  ################################
  # obtain bin annotations including binwise GC%, mappability, blacklist
  # available & pre-calculated for genome build
  # hg19/38 and bin sizes 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp in the
  # QDNAseq.hg19/38 package
  
  
  userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize, genome="hg38")
  
  ################################
  # create QDNAseqReadCounts object from BED file
  
  cat("Read files.\n")
  
  readCounts <- ReadFragmentsAsQDNAseqReadCounts(bins=userMadeBins,
                                   FragmentFiles=FragmentFiles,
                                   FragmentFileNames=FragmentFileNames)
  
  ################################
  # Normalize read counts so that samples with different coverage can be used together
  # for each sample, divide the read count values by the average read count
  cat("Normalization.\n")
  
  if(RCNormalization){
    # readCounts <- estimateCorrection(readCounts)
    readCopyNumber <- correctBins(readCounts,
                              method="none")
    
    readCopyNumber <- QDNAseq::normalizeBins(readCopyNumber,
                                        method=methodOfRCNormalization)
    
    readCounts <- new('QDNAseqReadCounts', 
                      bins=Biobase::featureData(readCopyNumber),
                      phenodata=Biobase::phenoData(readCopyNumber),
                      counts=Biobase::assayDataElement(readCopyNumber, "copynumber"))
  }
  
  ctrl <- readCounts
  
  ################################
  # filter out some chromosomes
  
  readCounts <- QDNAseq::applyFilters(readCounts, residual=FALSE, 
                                      blacklist=FALSE,
                                      mappability=FALSE, 
                                      bases=FALSE,
                                      chromosomes = c("X", "Y", "MT", "M"))
  
  ################################
  # calculate median residuals of the LOESS from the control dataset
  # (pre-calculated residuals from the bin annotations in  QDNAseq.hg19/38 package are overwrited. 
  # Only the precalculated shape and coordinates of the bins dataframe are used.)
  cat("Calculate LOESS fit.\n")
  
  userMadeBins$residual <- QDNAseq::iterateResiduals(readCounts,cutoff=10)
  
  chromosomes = c("X", "Y", "MT", "M")
  
  ################################
  # Create a residual filter from cfDNA controls. filtering conditions:
  ## not in the specified chromosomes
  ## gc% not NA
  ## residuals not NA
  cat("Create filter.\n")
  
  condition <- rep(TRUE, times=nrow(readCounts))
  condition <- !(Biobase::fData(readCounts)$chromosome %in% chromosomes)
  condition <- condition & !is.na(Biobase::fData(readCounts)$gc)
  
  residuals <- userMadeBins$residual
  cutoff <- TRUE * matrixStats::madDiff(residuals, na.rm=TRUE) 
  residualsMissing <- aggregate(residuals,
                                by=list(chromosome=Biobase::fData(readCounts)$chromosome),
                                function(x) all(is.na(x)))
  chromosomesWithResidualsMissing <-
    residualsMissing$chromosome[residualsMissing$x]
  chromosomesToInclude <-
    setdiff(chromosomesWithResidualsMissing, chromosomes)
  if (length(chromosomesToInclude) > 0) {
    message("Note: Residual filter missing for chromosomes: ",
            paste(chromosomesToInclude, collapse=", "))
    residuals[Biobase::fData(readCounts)$chromosome %in% chromosomesToInclude] <- 0
  }
  
  # If FALSE, filter the bin from analysis, if TRUE keep the bin.
  condition <- condition & !is.na(residuals)
  
  saveRDS(condition, 
          file=paste0(path_to_blacklist,"control_blacklist_", 
                      binSize, "Kbp","_norm_",RCNormalization,".RDS"))   
}
