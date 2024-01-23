####################################################################
# 1. Process fragment files for segmentation. Take     
# fragment files and do the same things as        
# CNAclinic::processForSegmentation does           
####################################################################
# 2. Generate reference sample using a set of     
# control samples and store in .csv file.         
####################################################################

library(QDNAseq.hg38)
library(this.path)
library(openssl)

source(paste0(trimws(this.dir()), "/00_bedprocessing.R"))


genRefName <- function(fileNames, prefix="Ref_", hashit=FALSE){
  # sort fileNames to keep it independant from the order

  ref_str <- paste(sort(fileNames), collapse = '_')
  if (hashit || length(fileNames)>7) {  # could choose other criteria for deciding whether to hash it or not i.e. length
    ref_ID <- md5(ref_str)
  } else {
    ref_ID <- ref_str
  }
  RefName <- make.names(paste0(prefix, ref_ID))
  cat("genRefName refName: ", RefName, "\n\t for: ", ref_str, "\n")
  return(RefName)
}

processFragmentFilesForSegmentation <- function(
  FragmentFiles=NULL, FragmentFileNames=NULL, refSamples=NULL,
  CountFiles=NULL,CountFileNames=NULL,
  binSize=NULL, genome="hg38", outputType="CNAclinicData", RCNormalization=TRUE,
  typeOfPreMadeBins="SR50",
  userMadeBins=NULL,
  residualFilter=TRUE, blacklistFilter=TRUE, mappabilityFilter=15,
  chromosomesFilter=c("X", "Y", "M", "MT"),
  spanForLoess=0.65, familyForLoess="symmetric",
  maxIterForCorrection=1, cutoffForCorrection=4.0,
  variablesForCorrection=c("gc", "mappability"), methodOfCorrection="ratio",
  methodOfNormalization="median", methodOfRCNormalization="mean",
  logTransformForSmoothing=TRUE,
  skipMedianNormalization=FALSE, skipOutlierSmoothing=FALSE,
  saveCountData=FALSE, filename="corrected_QDNAseqCopyNumbers"){
  
  InputFiles = c(FragmentFiles,CountFiles)
  InputFileNames = c(FragmentFileNames,CountFileNames)
  ############################################################################
  # Check user specified arguments
  ############################################################################
  
  # if FragmentFileNames are present but not distinct complain!
  if(!is.null(FragmentFileNames)){
    if(!(length(unique(FragmentFileNames)) == length(FragmentFiles))){
      errorMsg <- paste("Please provide unique FragmentFileNames:[", paste(unique(FragmentFileNames), collapse = ','), "] with length ",
                        length(unique(FragmentFileNames)), " that match the order of the FragmentFiles with length ", length(FragmentFiles))
      stop(errorMsg)
    }
  }
  # if CountFileNames are present but not distinct complain!
  if(!is.null(CountFileNames)){
    if(!(length(unique(CountFileNames)) == length(CountFiles))){
      errorMsg <- "Please provide unique CountFileNames
      that match the order of the CountFiles"
      stop(errorMsg)
    }
  }
  
  if(!is.null(refSamples)){
    if(!all(refSamples %in% c(FragmentFileNames,CountFileNames, NA, 'NA', 'drop'))
       || length(refSamples) != (length(InputFileNames))){
      errorMsg <-
        "refSamples should contain reference sample names
      that are to be used in normalizing each sample in FragmentFileNames and CountFileNames.
      If not NULL, refSamples must be the same length as FragmentFileNames + CountFileNames
      and should only include sample names contained in FragmentFileNames, CountFileNames, NA or
      'drop'.
      \n
      When 'drop', the corresponding sample from FragmentFileNames/CountFileNames will be removed
      from the output. When NA, the corresponding sample will be
      normalized by its median.
      \n
      As an example, if FragmentFileNames = c('tumour1', 'tumour2', 'normal2') and
      refSamples = c(NA, 'normal2', 'drop'),
      tumour1 will be kept as is, since it does not have a matched normal
      tumour2 will be divided by its matched reference: normal2 and
      normal2 will be dropped from further analysis."
      
      stop(errorMsg)
    }else{
      refSampleIndex = vector("numeric", length = (length(InputFileNames)))
      for(i in 1:length(refSamples)){
        ref <- refSamples[i]
        if(ref %in% c(NA, 'NA')){
          refSampleIndex[i] <- NA
        }else if(ref %in% c('drop')){
          refSampleIndex[i] <- FALSE
        }
        else{
          refSampleIndex[i] <- which(ref == InputFileNames)
        }
      }
    }
  }
  
  if(!is.null(userMadeBins)){
      if(!(binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000))){
        errorMsg <- paste(errorMsg, "bin size. Available bin sizes are:
                          hg38: 1, 5, 10, 15, 30, 50, 100, 500, 1000 Kbp.\n")
      
      stop(paste(errorMsg,
                 "CNAclinic does not generate the required annotation as this is
                 a time consuming step.
                 Please generate bin annotations using the QDNAseq package
                 and pass in the output via the userMadeBins argument."))
    }
    }else{
      
      if(is.null(binSize)){
        errorMsg <- "Please specify binSize argument"
        stop(errorMsg)
      }
    }
  
  outputType <- match.arg(outputType,
                          choices=c("CNAclinicData", "QDNAseqCopyNumbers"),
                          several.ok=FALSE)
  
  ############################################################################
  # Rename arguments to fit QDNAseq parameters
  ############################################################################
  binSize <- binSize[1]
  type <- typeOfPreMadeBins
  
  residual <- residualFilter
  blacklist <- blacklistFilter
  
  chromosomes <- chromosomesFilter
  chromosomes <- unique(c(chromosomes, "MT", "M"))
  
  span <- spanForLoess
  
  maxIter <- maxIterForCorrection
  cutoff <- cutoffForCorrection
  variables <- variablesForCorrection
  fit <- NULL
  logTransform <- logTransformForSmoothing
  ############################################################################
  
  if(is.null(userMadeBins)){
    
    # Read in pre-made QDNAseq bin annotation
    userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize,
                                               genome=genome,
                                               type=type, path = NULL)
  }
  
  readCountsFromF <- ReadFragmentsAsQDNAseqReadCounts(bins=userMadeBins,
                                                 FragmentFiles=FragmentFiles,
                                                 FragmentFileNames=FragmentFileNames)
  
  readCountsFromC <- ReadCountsAsQDNAseqReadCounts(bins=userMadeBins, 
                                CountFiles=CountFiles, 
                                CountFileNames=CountFileNames)
  
  readCounts <- CombineQDNAseqReadCounts(readCountsFromF,readCountsFromC)

  # Normalize read counts so that samples with different coverage can be used together
  # for each sample, divide the read count values by the average read count
  
  if(RCNormalization){
    readCopyNumber <- correctBins(readCounts,
                                  method="none")
    
    readCopyNumber <- QDNAseq::normalizeBins(readCopyNumber,
                                             method=methodOfRCNormalization)
    
    readCounts <- new('QDNAseqReadCounts', 
                      bins=Biobase::featureData(readCopyNumber),
                      phenodata=Biobase::phenoData(readCopyNumber),
                      counts=Biobase::assayDataElement(readCopyNumber, "copynumber"))
  }
  
  # Apply Blacklist filters
  readCountsFiltered <- QDNAseq::applyFilters(readCounts,
                                              residual=residual,
                                              blacklist=blacklist,
                                              mappability=mappabilityFilter,
                                              chromosomes=unique(c(chromosomes, "X", "Y")))
  
  # Estimate the correction for GC content and mappability
  readCountsFiltered <- QDNAseq::estimateCorrection(readCountsFiltered,
                                                    span=span,
                                                    family=familyForLoess,
                                                    adjustIncompletes=TRUE,
                                                    maxIter=maxIter,
                                                    cutoff=cutoff,
                                                    variables=variables)
  
  rm(list=c('readCounts'))
  gc(verbose=FALSE)
  
  if(!all(c("X", "Y") %in% chromosomes)){
    # Reverse the filtering for sex-chromosomes
    # to avoid confounding due to gender when estimating loess correction
    readCountsFiltered <- QDNAseq::applyFilters(readCountsFiltered,
                                                residual=residual,
                                                blacklist=blacklist,
                                                mappability=mappabilityFilter,
                                                chromosomes=NA)
  }
  
  
  # Apply the correction for GC content and mappability.
  # This returns a QDNAseqCopyNumbers object which is used to
  # normalise and smooth outliers
  copyNumbers <- QDNAseq::correctBins(readCountsFiltered,
                                      fit=fit,
                                      method=methodOfCorrection,
                                      adjustIncompletes=TRUE)
  
  if(saveCountData){
    
    saveRDS(copyNumbers, file=filename)
    
  }
  
  rm(list=c('readCountsFiltered'))
  gc(verbose=FALSE)
  
  #Apply median normalisation (Optional)
  
  if(!skipMedianNormalization)
    copyNumbers <- QDNAseq::normalizeBins(copyNumbers,
                                          method=methodOfNormalization)
  
  gc(verbose=FALSE)
  
  #Smooth outliers (Optional)
  if(!skipOutlierSmoothing)
    copyNumbers <- QDNAseq::smoothOutlierBins(copyNumbers,
                                              logTransform=logTransform)
  
  gc(verbose=FALSE)
  
  # Normalize by reference samples
  if(!is.null(refSamples)){
    copyNumbers <- QDNAseq::compareToReference(
      copyNumbers,
      references=refSampleIndex)
  }
  
  message("Pre-processing steps are completed.")
  
  if(outputType == "CNAclinicData"){
    
    bins <- data.frame(
      chromosome=QDNAseq::chromosomes(copyNumbers),
      start=QDNAseq::bpstart(copyNumbers),
      end=QDNAseq::bpend(copyNumbers),
      usebin=QDNAseq:::binsToUse(copyNumbers))
    
    sampleNames <- Biobase::sampleNames(copyNumbers)
    
    copyNumber<-as.data.frame(
      QDNAseq:::log2adhoc(
        Biobase::assayDataElement(copyNumbers, "copynumber")))
    
    names(copyNumber) <- sampleNames
    
    return(CNAclinic::CNAclinicData(bins=bins, copyNumber=copyNumber))
    
  }else{
    return(copyNumbers)
  }
}


GenerateRefSample <- function(ControlFiles, ControlFileNames, RefName,
                              OutputDir, genome="hg38",
                              RCNormalization=TRUE,
                              methodOfRCNormalization="mean",
                              bins=NULL,
                              binSize=NULL){
  print(paste("GenerateRefSample",RefName, OutputDir, genome, RCNormalization, bins, binSize))
  if(!is.null(bins)){
    if(!(binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000))){
      errorMsg <- paste(errorMsg, "bin size. Available bin sizes are:
                        hg38: 1, 5, 10, 15, 30, 50, 100, 500, 1000 Kbp.\n")
      
      stop(paste(errorMsg,
                 "CNAclinic does not generate the required annotation as this is
                 a time consuming step.
                 Please generate bin annotations using the QDNAseq package
                 and pass in the output via the userMadeBins argument."))
    }
    }else{
      
      if(is.null(binSize)){
        errorMsg <- "Please specify binSize argument"
        stop(errorMsg)
      }
    }
  
  if(is.null(bins)){
    
    # Read in pre-made QDNAseq bin annotation
    userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize,
                                               genome=genome)
  }
  
  readCounts <- ReadFragmentsAsQDNAseqReadCounts(bins=userMadeBins,
                                                 FragmentFiles=ControlFiles,
                                                 FragmentFileNames=ControlFileNames)
  
  # Normalize read counts so that samples with different coverage can be used together
  # for each sample, divide the read count values by the average read count
  
  if(RCNormalization){
    readCopyNumber <- correctBins(readCounts,
                                  method="none")
    
    readCopyNumber <- QDNAseq::normalizeBins(readCopyNumber,
                                             method=methodOfRCNormalization)
    
    readCounts <- new('QDNAseqReadCounts', 
                      bins=Biobase::featureData(readCopyNumber),
                      phenodata=Biobase::phenoData(readCopyNumber),
                      counts=Biobase::assayDataElement(readCopyNumber, "copynumber"))
  }
  
  counts <- Biobase::assayDataElement(readCounts, "counts")
  condition <- readCounts@featureData@data[["use"]]
  
  # take the median
  median.counts <- matrix(NA_integer_, nrow=nrow(counts), ncol=4,
         dimnames=list(Biobase::featureNames(readCounts), c("chromosome","start","end",RefName)))
  median.counts[,4] <- apply(counts,MARGIN=1, 
                         FUN = median, na.rm=TRUE)
  median.counts[!condition,] <- NA
  
  # save median.counts as _count.csv
  
  bins <- strsplit(rownames(median.counts),":|-")
  bins <- matrix(unlist(bins), nrow=length(bins), byrow=T)
  
  filename <- paste0(OutputDir,RefName,"_norm",RCNormalization,
                     "_",binSize,"kbp_",genome,"_count.csv")
  
  median.counts[,1:3] <- bins[,1:3]
  
  write.table(median.counts,filename,
              sep="\t",row.names = FALSE,
              col.names = TRUE,quote = FALSE)
}

######################################################
## test
#############
# 
# ControlFiles <- c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV01_point.csv",
#                   "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV02_point.csv")
# ControlFileNames <- c("HV01","HV02")
# RefName <- "HV_test"
# OutputDir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# binSize <- 1000
# 
# # GenerateRefSample(ControlFiles=ControlFiles,
# #                   ControlFileNames=ControlFileNames,
# #                   RefName=RefName,
# #                   OutputDir=OutputDir,
# #                   binSize=binSize)
# 
# FragmentFiles <- c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV01_point.csv",
#                   "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV02_point.csv")
# FragmentFileNames <- c("HV01","HV02")
# refSamples <- c("HV_test","HV_test",NA)
# CountFiles=c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/HV_test_normTRUE_1000kbp_hg38_count.csv")
# CountFileNames=c("HV_test")
# 
# test <- processFragmentFilesForSegmentation(
#   FragmentFiles=FragmentFiles, FragmentFileNames=FragmentFileNames, refSamples=refSamples,
#   CountFiles=CountFiles,CountFileNames=CountFileNames,
#   binSize=binSize)
