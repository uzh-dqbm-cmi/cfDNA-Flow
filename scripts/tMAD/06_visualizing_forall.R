# Plot the CBS segment values for all samples read from 
# a CNAData saved in RDS files

# constants: default
ref_hashit <- FALSE

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=3) {
  sample_dir <- args[1]   # input directory with BED files
  HVFragmentFileNames <- readLines(args[2])  # List of healthy samples: ctrl_sample_name
  binSize <- as.integer(args[3]) # Reads binned into binSize Kbp windows
  do_control_refs <- as.integer(args[4])  #
  samples <- read.table(file = args[5], sep = '\t', header = TRUE)
  if (length(args)>=6) {
    OutputDir <- args[6]
  } else {
    OutputDir <- sample_dir
  }
} else {
  ## test
  sample_dir <-                         "/Users/todor/data/cfdna/cfDNA-KiSpi/szsel_90_150/tMAD/";
  HVFragmentFileNames <- readLines("/Users/todor/data/cfdna/cfDNA-KiSpi/control_samples_radiooncology.txt")
  samples <- read.table(file =          "/Users/todor/data/cfdna/cfDNA-KiSpi/samples_KiSpiRH.tsv", sep = '\t', header = TRUE)
  #sample_dir <- "/Users/todor/data/cfdna/radiooncology/BED/szsel_90_150/tMAD/"
  #sample_dir <- "/Users/todor/data/cfdna/radiooncology/BED/szsel_200_250/tMAD/"; ref_hashit <- FALSE
  #sample_dir <- "/Users/todor/data/cfdna/radiooncology/BED/RAH_tMAD/"; ref_hashit <- TRUE
  #HVFragmentFileNames <- readLines("scripts/data/control_samples_radiooncology.txt")
  #sample_dir <- "/Users/todor/data/cfdna/AmsterdamUMC/results/BED/FalseD25630/szsel_90_150/tMAD/"
  #HVFragmentFileNames <- readLines("scripts/data/control_samples_RAH.txt")
  binSize <- 1000
  do_control_refs <- TRUE
  #samples <- read.table(file = "scripts/data/samples_HN_HV.tsv", sep = '\t', header = TRUE)
  #samples <- read.table(file = "scripts/data/samples_radiooncology.tsv", sep = '\t', header = TRUE)
  #samples <- read.table(file = "scripts/data/samples_RAH.tsv", sep = '\t', header = TRUE)
  #samples <- read.table(file = "scripts/data/samples_AmsterdamUMC.tsv", sep = '\t', header = TRUE)
  OutputDir <- sample_dir
  #stop("Unexpected number of arguments. usage: \n Rscript --vanilla 05_cal_t_MAD_forall.R <ctrl_sample_dir> <ctrl_sample_names.txt> <binZize> <RCNormalization:true,false> <genome:hg38> <sample.tsv>",
  #     call.=FALSE)
}

# ref_hashit <- FALSE # for newer results set to TRUE as in previous steps

library(this.path)

source(paste0(trimws(this.dir()), "/06_Visualizing.R"))

# HVFragmentFileNames <- c("HV01","HV03","HV04","HV05","HV06","HV07","HV08")


# HVRefNames <- c("Ref_HV03_HV04_HV05_HV06_HV07",
#                 "Ref_HV01_HV04_HV05_HV06_HV07",
#                 "Ref_HV01_HV03_HV05_HV06_HV07",
#                 "Ref_HV01_HV03_HV04_HV06_HV07",
#                 "Ref_HV01_HV03_HV04_HV05_HV07",
#                 "Ref_HV01_HV03_HV04_HV05_HV06",
#                 "Ref_HV01_HV03_HV04_HV05_HV06_HV07")

# binSize <- 1000
# OutputDir <- "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/results_new_reads/analysis/"

# targets
FragmentFileNames <- as.array(as.list(samples["sample_name"])[[1]])

PVRefName <- genRefName(HVFragmentFileNames, "Ref_", ref_hashit)    # "Ref_HV01_HV03_HV04_HV05_HV06_HV07"

controls_pattern<-paste0(paste0("^", HVFragmentFileNames), collapse = '|')
HVIndex <- grepl(controls_pattern, FragmentFileNames)

# controls
if (do_control_refs && length(HVFragmentFileNames) > 1) {
  for (h in 1:length(HVFragmentFileNames)) {
    control_refName <- genRefName(HVFragmentFileNames[-h], "Ref_", ref_hashit)  # same as in 05_cal_t_MAD_forall.R
    controls <- c("none", control_refName)
    plot_CNV(controls,binSize,OutputDir,HVFragmentFileNames[h])
  }
}
PVIndex <- !HVIndex
PVFragmentFileNames <- FragmentFileNames[PVIndex]  #  c("PV01","PV03","OMD1_PreRT","OMD001_d2","OMD001_d5")
# selected only
#PVFragmentFileNames <- c("OMD004")
#PVFragmentFileNames <- c("OMD009-BL", "OMD009-d2", "OMD009-d5", "OMD009-3m")

controls <- c("none", PVRefName) # for szsel_90_150
#controls <- c("none") #  for szsel_380_1000
#controls <- c("Ref") #  for Amsterdam
plot_CNV(controls,binSize,OutputDir,PVFragmentFileNames)
