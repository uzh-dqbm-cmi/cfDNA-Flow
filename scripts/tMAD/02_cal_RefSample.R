####################################################################
####### Create reference from control samples for each sample ######
####################################################################

library(this.path)

source(paste0(trimws(this.dir()), "/03_BeforeSegmentation.R"))

###########################################
# Arguments to change
###########################################

# RefDir <- "D:/Courses/2020_fall/cfDNA_analysis_sp2/cfDNA-pregnancy/results/"
# input directory with BED files
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=5) {
  RefDir <- args[1]   # input directory with BED files: ctrl_sample_dir
  HVFragmentFileNames <- readLines(args[2])  # List of healthy samples: ctrl_sample_name
  binSize <- as.integer(args[3])  # Reads binned into 1000 Kbp windows
  RCNormalization <- as.logical(args[4])
  genome <- args[5]
  do_control_refs <- as.logical(args[6])  #
} else {
  # for local test
  RefDir <- "/Users/todor/data/cfdna/cfDNA-KiSpi/szsel_90_150/"
  HVFragmentFileNames <- readLines("/Users/todor/data/cfdna/cfDNA-KiSpi/control_samples_KiSpi.txt")
  # sample_dir <- "test/data/BED/"
  # HVFragmentFileNames <- readLines("scripts/data/control_samples_RAH.txt")
  binSize <- 1000
  RCNormalization <- TRUE
  genome <- "hg38"
  do_control_refs <- TRUE
  #stop("Unexpected number of arguments. usage: \n Rscript --vanilla 02_cal_RefSample.R <ctrl_sample_dir> <ctrl_sample_names.txt> <binZize> <RCNormalization:true,false> <genome:hg38>",
  #    call.=FALSE)
}

print(paste("02_cal_RefSample.R arguments:", RefDir, HVFragmentFileNames, binSize, RCNormalization, genome, do_control_refs))
# RefDir <- "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"
# binSize <- 1000
FragmentFileDir <- RefDir
# healthy samples
# HVFragmentFileNames <- c("HV01","HV03","HV04","HV05","HV06","HV07")

# RCNormalization <- FALSE
cat(paste0("Create reference for sample:\n",paste(HVFragmentFileNames, collapse = ', '),'\n'))
HVFragmentFiles <- c()
for (i in HVFragmentFileNames) {
  HVFragmentFiles <- c(HVFragmentFiles,
                       paste0(FragmentFileDir,i,"_point.csv"))
}


if (do_control_refs && length(HVFragmentFileNames) > 1) {
  for (h in 1:length(HVFragmentFileNames)){

    ControlFileNames <- HVFragmentFileNames[-h]
    ControlFiles <- HVFragmentFiles[-h]
    RefName <- genRefName(ControlFileNames, "Ref_")  # paste0("Ref_",paste(ControlFileNames, collapse = '_'))
    cat(paste0("Reference name: ",RefName,"\n"))

    GenerateRefSample(ControlFiles=ControlFiles,
                      ControlFileNames=ControlFileNames,
                      RefName=RefName,
                      OutputDir=RefDir, genome=genome,
                      RCNormalization=RCNormalization,
                      methodOfRCNormalization="mean",
                      binSize=binSize)
  }
}

cat(paste0("Create reference for cancer samples from:\n",paste(HVFragmentFileNames, collapse = ', '),'\n'))
RefName <- genRefName(HVFragmentFileNames, "Ref_")  # paste0("Ref_",paste(HVFragmentFileNames, collapse = '_'))
cat(paste0("Reference name: ",RefName,"\n"))

GenerateRefSample(ControlFiles=HVFragmentFiles, 
                  ControlFileNames=HVFragmentFileNames, 
                  RefName=RefName,
                  OutputDir=RefDir, genome=genome,
                  RCNormalization=RCNormalization,
                  methodOfRCNormalization="mean",
                  binSize=binSize)
    