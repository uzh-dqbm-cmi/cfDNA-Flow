# process -> segmentation -> blacklist -> calculate t-MAD
# creates two .csv files with t-MAD values per sample
# The _control_none.csv file is the analyses run without a control (sample is median-normalized)
# The second file is the analyses run with an external control or panel
###########################################
# Arguments to change
###########################################

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=6) {
  sample_dir <- args[1]   # input directory with BED files
  HVFragmentFileNames <- readLines(args[2])  # List of healthy samples: ctrl_sample_name
  selected_binSize <- c(as.integer(args[3])) # Reads binned into 1000 Kbp windows
  RCNormalization <- as.logical(args[4])
  blacklist_Normalization <- RCNormalization
  genome <- args[5]
  do_control_refs <- as.logical(args[6])  #
  samples <- read.table(file = args[7], sep = '\t', header = TRUE)
  if (length(args)>=8) {
    OutputDir <- args[8]
  } else {
    OutputDir <- sample_dir
  }
} else {
  # /Users/todor/data/cfdna/cfDNA-KiSpi/szsel_90_150/ /Users/todor/data/cfdna/cfDNA-KiSpi/control_samples_KiSpi.txt 1000 TRUE hg38 /Users/todor/data/cfdna/cfDNA-KiSpi/samples_KiSpi.tsv /Users/todor/data/cfdna/cfDNA-KiSpi/szsel_90_150/tMAD/
  ## test
  sample_dir <- "/Users/todor/data/cfdna/cfDNA-KiSpi/szsel_90_150/"
  HVFragmentFileNames <- readLines("/Users/todor/data/cfdna/cfDNA-KiSpi/control_samples_KiSpi.txt")
  samples <- read.table(file = "/Users/todor/data/cfdna/cfDNA-KiSpi/samples_KiSpi.tsv", sep = '\t', header = TRUE)
  # sample_dir <- "test/data/BED/"
  # HVFragmentFileNames <- readLines("scripts/data/control_samples_RAH.txt")
  selected_binSize <- 1000
  RCNormalization <- TRUE
  blacklist_Normalization <- RCNormalization
  genome <- "hg38"
  do_control_refs <- TRUE
  #samples <- read.table(file = "scripts/data/samples_RAH.tsv", sep = '\t', header = TRUE)
  OutputDir <- paste0(sample_dir,"tMAD/")
  #stop("Unexpected number of arguments. usage: \n Rscript --vanilla 05_cal_t_MAD_forall.R <ctrl_sample_dir> <ctrl_sample_names.txt> <binZize> <RCNormalization:true,false> <genome:hg38> <sample.tsv>",
  #     call.=FALSE)
}

library(CNAclinic)
library(this.path)
library(hash)

source(paste0(trimws(this.dir()), "/05_cal_t_MAD.R"))

# 1) Path to all your sample files

# sample_dir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
#sample_dir <- "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"
FragmentFileSuffix <- "_point.csv"
# names of all samples, both tumor and healthy
# take it from the samples.tsv file
FragmentFileNames <- as.array(as.list(samples["sample_name"])[[1]])

#FragmentFileNames <- c("HV01","HV03","HV04","HV05",
#                       "HV06","HV07","HV08","PV01","PV03",
#                       "OMD1_PreRT","OMD001_d2","OMD001_d5")

# 2) Which genomic bin size should be used as resolution, scale or vector
# selected_binSize <- c(1000)
# genome = "hg38"
# 3) The path to blacklist created in the previous step
path_to_blacklist <- sample_dir # "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"

# 4) Whether the control samples are normalized or not 
# when calculating blacklist
# blacklist_Normalization <- TRUE # same as other script 01

# 5) If TRUE, normalize read counts before any processing 
# so that samples with different coverage can be used together

# RCNormalization <- FALSE  # should be same as in  other scripts == blacklist_Normalization
# directory with healthy samples, BED files
RefDir <- sample_dir  # "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"
HVRefNames <- hash()

if (do_control_refs && length(HVFragmentFileNames) > 1) {
  for (h in 1:length(HVFragmentFileNames)){
    ControlFileNames <- HVFragmentFileNames[-h]
    RefName <- genRefName(ControlFileNames, "Ref_")  # paste0("Ref_",paste(ControlFileNames, collapse = '_'))
    HVRefNames[[HVFragmentFileNames[h]]] <- RefName
    cat('control HVFragmentFileNames:',HVFragmentFileNames[h],' RefName:', RefName)
  }
}
PVRefName <- genRefName(HVFragmentFileNames, "Ref_")   # paste0("Ref_",paste(HVFragmentFileNames, collapse = '_'))
#PVRefName <- "Ref_HV01_HV03_HV04_HV05_HV06_HV07" # all healthy
###########################################
# Check the code below for necessary changes
###########################################


segType <- "CH"

for(t in 1:length(selected_binSize)){
  
  binSize <- selected_binSize[t]
  
  cat('Calculate tMAD for bin size ',binSize, 'kbp.\n')
  
  FragmentFiles <- c()
  for (i in FragmentFileNames) {
    FragmentFiles <- c(FragmentFiles, paste0(sample_dir,i,FragmentFileSuffix))
  }
  controls_pattern<-paste0(paste0("^", HVFragmentFileNames), collapse = '|')
  HVIndex <- grepl(controls_pattern, FragmentFileNames)
  HVFragmentFileNames <- FragmentFileNames[HVIndex]  # should be == previous HVFragmentFileNames
  HVFragmentFiles <- FragmentFiles[HVIndex]
  
  PVIndex <- !HVIndex
  PVFragmentFileNames <- FragmentFileNames[PVIndex]
  PVFragmentFiles <- FragmentFiles[PVIndex]

  if (do_control_refs && length(HVFragmentFileNames) > 1) {
    for (h in 1:length(HVFragmentFileNames)){
      sample_name <- HVFragmentFileNames[h]
      sample_file <- HVFragmentFiles[h]
      RefName <- HVRefNames[[sample_name]]
      HVRefFile <- paste0(RefDir,RefName,"_norm",RCNormalization,
                          "_",binSize,"kbp_",genome,"_count.csv")

      FragmentFiles <- c(sample_file)
      # now all balck list controls must be in the sample.tsv file must be in the same order as in control_sample.txt / best sorded
      CountFile <- c(HVRefFile)
      FragmentFileNames <- c(sample_name)
      CountFileName <- c(RefName)

      cat('Calculate tMAD for sample(s) ',
          paste(FragmentFileNames,collapse=', '),
          'with and without control. \n')
      cat('Control file: ',CountFileName,'\n')
      if (is.na(CountFileName)){
        next # skip not used in the black list
      }
      cal_t_MAD(path_to_blacklist = path_to_blacklist,
                blacklist_Normalization = blacklist_Normalization,
                binSize = binSize,OutputDir = OutputDir,
                FragmentFiles = FragmentFiles,
                FragmentFileNames = FragmentFileNames,
                isControl = TRUE,
                the_control = CountFile,
                control_name = CountFileName,
                genome=genome,
                segmentType=c("CBS"),
                summaryMethod="mean",
                chromosomesFilter=c("X", "Y", "MT", "M"),
                saveCountData=FALSE,
                RCNormalization=RCNormalization,
                segType=segType)
    }
  }

  PVRefFile <- paste0(RefDir,PVRefName,"_norm",RCNormalization,
                      "_",binSize,"kbp_",genome,"_count.csv")

  FragmentFiles <- c(PVFragmentFiles)
  FragmentFileNames <- c(PVFragmentFileNames)

  CountFile <- c(PVRefFile)
  CountFileName <- c(PVRefName)
  
  cat('Calculate tMAD for sample(s) ',
      paste(FragmentFileNames,collapse=', '),
      ' with and without control. \n')
  cat('Control file: ',CountFileName,'\n')

  cal_t_MAD(path_to_blacklist = path_to_blacklist,
            blacklist_Normalization = blacklist_Normalization,
            binSize = binSize,OutputDir = OutputDir,
            FragmentFiles = FragmentFiles,
            FragmentFileNames = FragmentFileNames,
            the_control = CountFile,
            control_name = CountFileName,
            genome=genome,
            segmentType=c("CBS"),
            summaryMethod="mean",
            chromosomesFilter=c("X", "Y", "MT", "M"),
            saveCountData=FALSE,
            RCNormalization=RCNormalization,
            segType=segType)
}
