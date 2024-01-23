####################################################################
# Calculate t-MAD for one sample or for multiple samples 
# using one reference file
####################################################################

library(this.path)

source(paste0(trimws(this.dir()), "/04_DoSegmentation.R"))

cal_t_MAD <- function(path_to_blacklist,blacklist_Normalization,
                      binSize, OutputDir,
                      the_control,control_name,
                      FragmentFiles=NULL,FragmentFileNames=NULL,
                      isControl = FALSE,
                      ...){

  controlSample <- c(control_name, "none")
  blacklist_file <- paste0(path_to_blacklist, 
                           "control_blacklist_", 
                           binSize, "Kbp","_norm_",blacklist_Normalization,".RDS")
  
  cfDNA_blacklist <- readRDS(blacklist_file)

  tMAD_results <- NULL

  for(j in 1:length(controlSample)){
    
    # Doing only median normalization for these samples        
    control <- controlSample[j]
    
    processedData <- NULL   # TDDO confusing asigment: check if needed

    prefix <- if (length(FragmentFileNames)>1) "ALL_" else ""
    sampleNamesRef <- genRefName(FragmentFileNames, prefix)
    # column_suffix <- paste0("bs", binSize,"_control_",control)  # ,"_samples_",sampleNamesRef)

    if(control == "none"){
      column_suffix <- "none"

      CNAData <- DoSegmentation(FragmentFiles=FragmentFiles,
                                FragmentFileNames=FragmentFileNames,
                                binSize=binSize,control_name="none",
                                CountFiles=the_control,
                                CountFileNames=control_name,
                                OutDir = OutputDir,
                                ...)

    }else if(control == control_name){
      column_suffix <- "control"

      CNAData <- DoSegmentation(FragmentFiles=FragmentFiles,
                                FragmentFileNames=FragmentFileNames,
                                binSize=binSize,control_name = control_name,
                                CountFiles=the_control,
                                CountFileNames=control_name,
                                skipMedianNormalization=TRUE,
                                OutDir = OutputDir,
                                ...)
    }
    
    sampleNames <- unlist(lapply(CNAclinic::sampleNames(CNAData), 
                                 function(x){strsplit(x, ".vs.")[[1]][1]}))
    
    colnames(CNAData@copyNumber) <- sampleNames
    
    stopifnot(length(cfDNA_blacklist) == length(CNAclinic::usebin(CNAData)))
    
    segData <-  CNAclinic::segSummary(CNAData)
    
    segData <- segData[cfDNA_blacklist & CNAclinic::usebin(CNAData), ]
    
    segData[abs(segData) >=5 ] <- NA
    if (is.null(dim(segData))){
      segData <- matrix(segData,ncol=length(CNAclinic::sampleNames(CNAData)))
    }
    tMAD <- apply(segData, 2,  function(x){ 
      abs(mad(x = x, center=0, na.rm = TRUE))
    })
    
    outData <- data.frame(sample=CNAclinic::sampleNames(CNAData),
                          tMAD)
    colnames(outData) <- c("sample", paste0("tmad_vs_", column_suffix))

    # write.csv(outData, file=tMADfile, quote=F, row.names=F)
    if (is.null(tMAD_results)){
      tMAD_results <- outData
    } else {
      tMAD_results <- merge(tMAD_results, outData, by.y = "sample")
    }
  }
  tMAD_results$tmad_is_control <- isControl

  tMADfile <- paste0(OutputDir,"tMAD_results.tsv")
  cat("save tMADs for control ", control , "and sampleNamesRef", sampleNamesRef,"in file", tMADfile, "\n")

  append <- file.exists(tMADfile)
  if (!dir.exists(dirname(tMADfile))){
    mkdir(dirname(tMADfile))
  }
  write.table(tMAD_results, file = tMADfile, row.names=FALSE, append = append, col.names = !append, sep="\t")
}

#######################################################
## Test
################

# FragmentFiles <- c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV01_point.csv",
#                    "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/data/HV02_point.csv")
# FragmentFileNames <- c("HV01","HV02")
# refSamples <- c("HV_test","HV_test",NA)
# CountFiles=c("D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/HV_test_normTRUE_1000kbp_hg38_count.csv")
# CountFileNames=c("HV_test")
# binSize <- 1000
# sample_dir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# FragmentFileSuffix <- "*point.csv"
# selected_binSize <- c(1000)
# the_control <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/HV_test_normTRUE_1000kbp_hg38_count.csv"
# control_name <- "HV_test"
# path_to_blacklist <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# blacklist_Normalization <- FALSE
# RCNormalization <- TRUE
# 
# CNADatas <- cal_t_MAD(path_to_blacklist=path_to_blacklist,
#           blacklist_Normalization=blacklist_Normalization,
#           binSize=binSize,
#           FragmentFiles=FragmentFiles,
#           FragmentFileNames=FragmentFileNames,
#           the_control=the_control,control_name=control_name,
#           OutputDir = path_to_blacklist)