####################################################################
# 1. Do segmentation
####################################################################

library(this.path)

source(paste0(trimws(this.dir()), "/03_BeforeSegmentation.R"))

DoSegmentation <- function(FragmentFiles=NULL,FragmentFileNames=NULL,
                           CountFiles=NULL,CountFileNames=NULL,
                           binSize,segType="CH",control_name="none",
                           genome="hg38",
                           segmentType=c("CBS"),
                           summaryMethod="mean",OutDir="./",
                           ...){
  
  InputFiles <- c(FragmentFiles,CountFiles)
  InputFileNames <- c(FragmentFileNames,CountFileNames)
  
  cat('DoSegmentation for samples: ', paste(FragmentFileNames,collapse=','),'\n')
  prefix <- if (length(FragmentFileNames)>1) "ALL_" else ""
  sampleNamesRef <- genRefName(FragmentFileNames, prefix)
  outfile_suffix <- paste0(binSize, "_control_",control_name,"_samples_", sampleNamesRef)

  if(control_name == "none"){

    processedData <- processFragmentFilesForSegmentation(
      FragmentFiles=FragmentFiles,
      FragmentFileNames=FragmentFileNames,
      CountFiles=CountFiles,
      CountFileNames=CountFileNames,
      binSize=binSize,
      ...)
    
    # Run segmentation
    if(segType == "CH"){
      # CNAData <- runSegmentationNew(processedData, genome=genome,
      #                            segmentType=segmentType,
      #                            summaryMethod=summaryMethod)
      
      CNAData <- CNAclinic::runSegmentation(processedData, genome=genome,
                                    segmentType=segmentType,
                                    summaryMethod=summaryMethod)
      RDSFile <- paste0(OutDir,"CNAData_", outfile_suffix,".RDS")
      cat("saveRDS CNAData for control ", control_name ,"in file", RDSFile, "\n")
      saveRDS(CNAData, file=RDSFile)
    }
  }else if(control_name != "none"){
    
    if(length(control_name) != 1)
      stop("Check control specified")
    
    which_control <- which(InputFileNames == control_name)
    refSamples <- rep(control_name, length(InputFileNames))
    
    # drop the control sample after normalizing test samples 
    # as we will have log2R=0 if it is normalized by its own bin counts
    
    refSamples[which_control] <- 'drop'

    cat('DoSegmentation control_name: ', control_name, '; FragmentFiles=', paste(FragmentFiles,collapse=','),';\n FragmentFileNames=', paste(FragmentFileNames,collapse=','), '\n' )

    processedData <- processFragmentFilesForSegmentation(FragmentFiles=FragmentFiles,
                                                         FragmentFileNames=FragmentFileNames,
                                                         binSize=binSize,
                                                         refSamples=refSamples,
                                                         CountFiles=CountFiles,
                                                         CountFileNames=CountFileNames,
                                                         ...)
    
    # Run segmentation
    if(segType == "CH"){
      # CNAData <- runSegmentationNew(processedData, genome=genome,
      #                            segmentType=segmentType,
      #                            summaryMethod=summaryMethod)
      # 
      CNAData <- CNAclinic::runSegmentation(processedData, genome=genome,
                                    segmentType=segmentType,
                                    summaryMethod=summaryMethod)
      RDSFile <- paste0(OutDir,"CNAData_", outfile_suffix,".RDS")
      cat("saveRDS CNAData for control ", control_name ,"in file", RDSFile, "\n")
      saveRDS(CNAData, file=RDSFile)
    }
  }
  CNAData
}


#######################################################
## Test
################
# 
# FragmentFiles <- c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV01_point.csv",
#                   "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV02_point.csv")
# FragmentFileNames <- c("HV01","HV02")
# refSamples <- c("HV_test","HV_test",NA)
# CountFiles=c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/HV_test_normTRUE_1000kbp_hg38_count.csv")
# CountFileNames=c("HV_test")
# binSize <- 1000
# OutDir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# segmentType=c("CBS") # "HMM" doesn't work, due to the bug of CNAclinic
# test <- DoSegmentation(
#   FragmentFiles=FragmentFiles,FragmentFileNames=FragmentFileNames,
#   CountFiles=CountFiles,CountFileNames=CountFileNames,
#   binSize=binSize,
#   segmentType=segmentType,
#   OutDir = OutDir)
