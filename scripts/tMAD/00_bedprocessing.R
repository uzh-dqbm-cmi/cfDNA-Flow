####################################################################
# 1. Bin fragments by chromosome. The FragmentFile     
# are supposed to have been preprocessed such that 
# unwanted fragments are not included      
####################################################################
# 2. Load Fragment information stored in .csv file into 
# QDNAseqReadCounts objects       
####################################################################
# 3. Load read count information stored in .csv file into 
# QDNAseqReadCounts objects       
####################################################################

binFragmentsPerSample <- function(bins, FragmentFile) {
  # Read one fragment file (fields: chromosome, position, ...) and bin by chromosome.

  cat('binFragmentsPerSample: read.table(',str(FragmentFile),')\n')
  fragments <- read.table(FragmentFile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  chrs <- unique(fragments[["V1"]])
  
  # if (any(grepl("chr",chrs_fragments))){
  #   chrs <- split(unlist(strsplit(chrs_fragments, split='chr', fixed=TRUE)),1:2)[[2]]
  # }
  
  readCounts <- integer(length=nrow(bins))
  for (chromosome in chrs) {
    keep_bins <- which(bins$chromosome == chromosome);
    
    ## No bins for this chromosome?
    if (length(keep_bins) == 0L)
      next
    
    chromosomeBreaks <- c(bins$start[keep_bins], max(bins$end[keep_bins]) + 1)
    keep_fragments <- which(fragments[["V1"]] == chromosome)
    counts <- matrixStats::binCounts(fragments[keep_fragments,"V2"], bx=chromosomeBreaks)
    readCounts[keep_bins] <- readCounts[keep_bins] + counts
    
    ## Not needed anymore
    chromosomeBreaks <- keep_bins <- counts <- keep_fragments <- NULL
  }
  ## Not needed anymore
  rm(list=c("fragments"))
  gc(FALSE)
readCounts
}

ReadFragmentsAsQDNAseqReadCounts <- function(bins, FragmentFiles, FragmentFileNames) {
  # Read a group of fragment files (fields: chromosome, position, ...), 
  # bin by chromosome, and store the read count information as QDNAseqReadCounts object.
  
  if (class(bins) == 'data.frame'){
    bins <- AnnotatedDataFrame(bins)
  }
  
  counts <- matrix(NA_integer_, nrow=nrow(bins), ncol=length(FragmentFileNames),
                   dimnames=list(Biobase::featureNames(bins), FragmentFileNames))
  
  for (i in seq_along(FragmentFiles)) {
    counts[, i] <- binFragmentsPerSample(bins,FragmentFiles[i])
  }
  condition <- bins@data[["use"]]
  phenodata <- Biobase::AnnotatedDataFrame(
    data.frame(
      sampleNames=colnames(counts),
      total.reads=colSums(counts),
      used.reads=matrixStats::colSums2(counts, rows=condition)
    ))
  
  object <- new('QDNAseqReadCounts', bins=bins, counts=counts,
                phenodata=phenodata)
  object
}

# write a new function to read _count.csv as the object

ReadCountsAsQDNAseqReadCounts <- function(bins, CountFiles, CountFileNames) {
  # Read a group of count files (fields: chromosome, start, end, count ...), 
  # and store as QDNAseqReadCounts object.
  
  if (class(bins) == 'data.frame'){
    bins <- AnnotatedDataFrame(bins)
  }
  
  counts <- matrix(NA_integer_, nrow=nrow(bins), ncol=length(CountFileNames),
                   dimnames=list(Biobase::featureNames(bins), CountFileNames))
  
  for (i in seq_along(CountFiles)) {
    count <- read.table(CountFiles[i],header = TRUE, 
                        sep="\t",stringsAsFactors=FALSE, 
                        quote="")
    counts[, i] <- count[,CountFileNames[i]]
  }
  condition <- bins@data[["use"]]
  phenodata <- Biobase::AnnotatedDataFrame(
    data.frame(
      sampleNames=colnames(counts),
      total.reads=colSums(counts),
      used.reads=matrixStats::colSums2(counts, rows=condition)
    ))
  
  object <- new('QDNAseqReadCounts', bins=bins, counts=counts,
                phenodata=phenodata)
  object
}

CombineQDNAseqReadCounts <- function(readCounts1,readCounts2){
  
  counts1 <- Biobase::assayDataElement(readCounts1, "counts")
  condition1 <- readCounts1@featureData@data[["use"]]
  
  counts2 <- Biobase::assayDataElement(readCounts2, "counts")
  condition2 <- readCounts2@featureData@data[["use"]]
  
  counts <- cbind(counts1,counts2)
  condition <- condition1 & condition2
  
  phenodata <- Biobase::AnnotatedDataFrame(
    data.frame(
      sampleNames=colnames(counts),
      total.reads=colSums(counts),
      used.reads=matrixStats::colSums2(counts, rows=condition)
    ))
  
  bins1 <- Biobase::featureData(readCounts1)
  bins2 <- Biobase::featureData(readCounts2)
  if(!all(Biobase::featureNames(bins1) == Biobase::featureNames(bins2))){
    errorMsg <-
      "Two objects have different bin annotations. Please check."
      stop(errorMsg)
  }else{
    bins <- bins1
  }
  object <- new('QDNAseqReadCounts', bins=bins, counts=counts,
                phenodata=phenodata)
  object
}