
library(this.path)

source(paste0(trimws(this.dir()), "/03_BeforeSegmentation.R"))

rStartsWith <- function(str, text){
  return(startsWith(text, str))
}

plot_CNV <- function(controls,binSize,OutputDir,FragmentFileNames){
  
  for (control_name in controls){
    genRefName(HVFragmentFileNames, "Ref_", ref_hashit)
    cat('plot_CNV for all samples: ', paste(FragmentFileNames,collapse=','),'\n')
    prefix <- if (length(FragmentFileNames)>1) "ALL_" else ""
    sampleNamesRef <- genRefName(FragmentFileNames, prefix)
    outfile_suffix <- paste0(binSize,"_control_",control_name,"_samples_",sampleNamesRef)
    
    CNAData <- readRDS(file=paste0(OutputDir,"CNAData_", outfile_suffix,".RDS"))
    sample_names <- CNAclinic::sampleNames(CNAData)
    # samples_pattern<-paste0(paste0("^", FragmentFileNames), collapse = '|')
    for (selected_sample in sample_names){
      # only for chery pick selection code:
      #if (!startsWith(selected_sample, "OMD")
      #  && !startsWith(selected_sample, "HN")
      #  && !startsWith(selected_sample, "PV")
      #) {
      #  next
      #}
      P8_SF_Bx <- CNAclinic::subsetData(CNAData, selected_sample)
      
      # Get log2R plot with summary segments and calls
      
      genomewide_plot <- CNAclinic::plotSampleData(P8_SF_Bx, 
                                                   showCalls=TRUE, 
                                                   segmentType="summary",
                                                   xaxSize=7,
                                                   mainSize=12)
      
      sample_plot <- gridExtra::grid.arrange(genomewide_plot[[1]], 
                                                  ncol=1)
      
      ggplot2::ggsave(sample_plot, width=10, height=5, filename=paste0(OutputDir,"CNV_sample_",selected_sample,".pdf"))
                            #_control_",control_name,
    }
  }
}
