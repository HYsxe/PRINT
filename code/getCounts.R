# All the functions needed for getting the group-by-region-by-position counts tensor from a fragment file

# Read from a fragments file and retrieve GRanges of fragments
fragsToRanges <- function (fragFile, # Path to the fragments file
                           barcodeList = NULL, # Whether to use a specific list of barcodes
                           startsAre0based = TRUE, # Whether the coordinates are 0-based 
                           nrows=Inf # Number of rows to read in
                           ) 
{
  cat("Reading in fragment file ..\n")
  if (is.null(barcodeList)) {
    
    cat("Reading all barcodes found within file ..\n")
    frags <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, nrows=nrows) %>% 
      data.frame() %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
                                                               start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
                                                               starts.in.df.are.0based = startsAre0based)
  }
  else {
    
    cat("Reading only select barcodes specified within list from file ..\n")
    frags <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, nrows=nrows) %>% 
      data.frame() 
    
    frags <- frags %>% filter(V4 %in% barcodeList) %>% 
      GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
                                              start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
                                              starts.in.df.are.0based = startsAre0based)
  }
  return(frags)
}

# Method to get countTensor from the footprintingProject object
setGeneric("countTensor", function(x) standardGeneric("countTensor"))

setMethod("countTensor", "footprintingProject", function(x) x@countTensor)

# Method to set countTensor from the footprintingProject object
setGeneric("countTensor<-", function(x, value) standardGeneric("countTensor<-"))

setMethod("countTensor<-", "footprintingProject", function(x, value) {
  x@countTensor <- value
  x
})

# Method to get the group-by-region-by-position counts tensor
setGeneric("getCountTensor",
           function(project, pathToFrags, barcodeGroups,...) standardGeneric("getCountTensor"))

setMethod("getCountTensor", 
          signature = c(project = "footprintingProject"),
          function(project, # footprintingProject object
                   pathToFrags, # Path to fragments file
                   barcodeGroups, # Data.frame specifying membership of barcodes in pseudobulks. First column is barcodes and second is groupID 
                   maxFragLength = NULL, # Fragment length upper limit
                   nrows = Inf, # Number of rows to read in
                   chunkSize = NULL, # Chunk size for parallel processing of regions
                   nCores = 16, # Number of cores to use
                   pairedEnd = T, # True/False. Whether data is paired-end or single-end
                   returnCombined = T # Whether to return the combined result for all chunks. Set it to False when data is too big
          ) {
            
            # Set project slots
            fragFile(project) <- pathToFrags
            barcodeGrouping(project) <- barcodeGroups
            
            # Retrieve genomic ranges of the regions to footprint
            regions <- regionRanges(project)
            
            # We chunk the regions into batches to reduce memory usage
            if(is.null(chunkSize)){
              chunkSize <- regionChunkSize(project)
            }
            print(paste0("Using chunk size = ", chunkSize))
            
            # Get countTensor
            countTensor(project) <- computeCountTensor(pathToFrags = pathToFrags, 
                                                       regions = regions, 
                                                       barcodeGroups = barcodeGroups, 
                                                       maxFragLength = maxFragLength,
                                                       tmpDir = dataDir(project),
                                                       nrows = nrows,
                                                       chunkSize = min(chunkSize, length(regions)),
                                                       nCores = nCores,
                                                       pairedEnd = pairedEnd,
                                                       returnCombined = returnCombined)
            
            project
          })

# Helper function for getCountTensor()
computeCountTensor <- function(pathToFrags, # Path to fragments file
                               regions, # Genomic ranges of the regions to footprint
                               barcodeGroups, # Data.frame specifying membership of barcodes in pseudobulks. First column is barcodes and second is groupID 
                               maxFragLength = NULL, # Fragment length upper limit
                               tmpDir = "./", # Directory to store results
                               nrows = Inf, # Max number of rows when reading from fragments file
                               chunkSize = 2000, # Chunk size for parallel processing of regions
                               nCores = 16, # Number of cores to use
                               pairedEnd = T, # True/False. Whether data is paired-end or single-end
                               returnCombined = T # Whether to return the combined result for all chunks. Set it to False when data is too big
) {
  
  start_time <- Sys.time()
  if (is.null(barcodeGroups)){
    stop("Need to supply barcode grouping ..\n")
  }
  
  # Rename columns in barcodeGrouping
  colnames(barcodeGroups) <- c("barcodeID", "group")
  
  #################### Load fragments data, convert to GRanges #######################
  
  cat("Make 1bp step .. \n")
  nRegions <- length(regions)
  
  # Retrieve genomic ranges of fragments
  frags <- fragsToRanges(pathToFrags, 
                         barcodeList = barcodeGroups$barcodeID, 
                         startsAre0based = TRUE, 
                         nrows=nrows)
  
  # Rename extra columns
  if (ncol(mcols(frags)) > 1) {
    colnames(mcols(frags)) <- c("barcodeID", "pcrDup")
  } else {
    colnames(mcols(frags)) <- "barcodeID"
  }
  
  # Filter fragments by size (optional)
  if (!is.null(maxFragLength)) {
    cat("Removing frags with length > ", maxFragLength, " bp ..\n")
    frags <- frags[width(frags) <= maxFragLength]
    if (length(frags) == 0) 
      stop("Fragment filtering resulting in 0 aligned fragments. Please check / change the provided filter size ..\n")
  }
  
  # Get group ID for each cell
  groupInds <- gtools::mixedsort(unique(barcodeGroups$group))
  gc()
  
  #################### Get region-by-position-by-pseudobulk Tn5 insertion count tensor #######################
  
  # To reduce memory usage, we chunk the data in to smaller chunks
  cat("Reformating counts data into a list (each element is data for a region) ..\n")
  chunkIntervals <- getChunkInterval(regions, chunkSize = chunkSize)
  starts <- chunkIntervals[["starts"]]
  ends <- chunkIntervals[["ends"]]
  
  # Create a folder for saving intermediate results
  chunkTmpDir <- paste0(tmpDir, "chunkedCountTensor/")
  if(!dir.exists(chunkTmpDir)){
    system(paste("mkdir -p", chunkTmpDir))
  }
  
  # For each chunk, we extract data for individual regions
  # For each region, we store the data in a 3-column data.frame (columns are group, position and counts)
  # Re-organize data into lists
  for(i in 1:length(starts)){
    
    print("Re-organizing data into lists")
    print(paste0(Sys.time()," Processing chunk ", i, " out of ", length(starts), " chunks"))
    
    # Skip current chunk if result already exists
    if(file.exists(paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))){
      next
    }
    
    # Get fragments within the current chunk
    chunkRegions <- starts[i]:ends[i]
    chunkFrags <- subsetByOverlaps(frags, regions[chunkRegions])
    
    # Go through each group (pseudobulk) and retrive cutsite data
    # Should result in a data table with 4 columns: Region index, position in this region, group ID, count
    groupedCountTensor <- pbmcapply::pbmclapply(
      groupInds,
      function(groupInd){
        
        # Get fragmenst belonging to the current cell group (i.e. pseudobulk)
        groupBarcodes <- barcodeGroups$barcodeID[barcodeGroups$group %in% groupInd]
        groupFrags <- chunkFrags[chunkFrags$barcodeID %in% groupBarcodes]
        
        if(pairedEnd){
          # Get all Tn5 insertion sites (single base pair resolution)
          # Note: The input fragments file should be +4/-5 shifted to accommodate common practice
          # However, +4/-5 actually points to the base immediately to the left of the center of the 9bp staggered end
          # The 1bp cut position should actually by +5/-4. Therefore we need to shift both start and end by +1
          # The function fragsToRanges already shifts start position by +1 when we specify "startsAre0based = T"
          # Therefore here we only need to further shift the end position by + 1
          cutsites <- c(resize(groupFrags, width = 1, fix = "start"),
                        IRanges::shift(resize(groupFrags, width = 1, fix = "end"), 1))
        }else{
          # If we were provided single-end data, the fragments file should be in the format of
          # (chr name, cut position, cut position, barcode, number of insertions at this pos)
          cutsites <- c(resize(groupFrags, width = 1, fix = "start"))
        }
        
        # Get position of cutsites within regions
        ovRegions <- findOverlaps(query = regions, 
                                subject = cutsites)
        if(length(ovRegions) == 0){
          groupCountTensor <- NULL
        }else{
          positions <- start(cutsites)[ovRegions@to] - start(regions)[ovRegions@from] + 1
          
          # Generate a data frame with 4 columns: Region index, position in this region, group ID, count
          cat(paste0("Generating matrix of counts for group ", groupInd, "..\n"))
          groupCountTensor <- Matrix::sparseMatrix(i = ovRegions@from,
                                                   j = positions,
                                                   x = 1)
          groupCountTensor <- summary(t(groupCountTensor))
          colnames(groupCountTensor) <- c("position", "region", "count")
          groupCountTensor$group <- groupInd
          groupCountTensor <- groupCountTensor[,c("region", "position", "group", "count")]
          groupCountTensor <- as_tibble(groupCountTensor)
        }
        
        groupCountTensor
      },
      mc.cores = nCores
    )
    groupedCountTensor <- data.table::rbindlist(groupedCountTensor)
    
    # Re-organize the data into a list where each element is the data for a single region.
    if(dim(groupedCountTensor)[1] > 0){
      cluster <- prep_cluster(length(chunkRegions), n_cores = nCores)
      opts <- cluster[["opts"]]
      cl <- cluster[["cl"]]
      countTensorChunk <- foreach(regionInd = chunkRegions,
                                  .options.snow = opts, 
                                  .packages = c("dplyr","Matrix")) %dopar%   {
                                    regionCountTensor <- groupedCountTensor %>% filter(region %in% regionInd)
                                    return(regionCountTensor)
                                  }
      stopCluster(cl)
    }else{
      countTensorChunk <- lapply(chunkRegions, function(x){data.table::data.table()})
    }
    
    # Save results
    saveRDS(countTensorChunk,
            paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))
    
    # Release unused memory
    if((i %% 10) == 0) gc()
  }
  
  cat("Done!\n")
  end_time <- Sys.time()
  cat("Time elapsed: ", end_time - start_time, units(end_time - start_time), " \n\n")
  
  if(returnCombined){
    countTensorAll <- pbapply::pblapply(
      1:length(starts),
      function(i){
        readRDS(paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))
      }
    )
    countTensorAll <- Reduce(c, countTensorAll)
  }else{
    countTensorAll <- list()
  }
  return(countTensorAll)
  
}

# One-hot encoding of a single DNA sequence
onehotEncode <- function(seq){
  onehotMapping <- list("A" = 1, "C" = 2, "G" = 3, "T" = 4, "N" = 0)
  len <- nchar(seq)
  onehotInd <- sapply(
    substring(seq, 1:len , 1:len),
    function(x){onehotMapping[[x]]})
  onehot <- matrix(0, nrow = 4, ncol = len)
  undeterminedBases <- onehotInd == 0
  onehot[cbind(unname(onehotInd), 1:len)[!undeterminedBases,]] <- 1
  onehot
}

# Use PWM to score a sequence
PWMScoring <- function(seq, PWM){
  onehotSeq <- onehotEncode(seq)
  sum(PWM * onehotSeq)
}

# Method to get the Tn5 bias PWM from the current dataset
setGeneric("getTn5PWM",
           function(project,...) standardGeneric("getTn5PWM"))

setMethod("getTn5PWM", 
          signature = c(project = "footprintingProject"),
          function(project, # footprintingProject object
                   PWMRadius = 10, # Radius of the PWM to generate
                   nSample = 1000 # Number of regions to sample
          ){
            # Check if the region-by-position matrix has been computed
            if(dim(ATACTracks(project))[1] == 0){
              project <- getATACTracks(project)
            }
            
            # Get reference genome
            referenceGenome <- refGenome(project)
            if(referenceGenome == "hg19"){
              genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
            }else if(referenceGenome == "hg38"){
              genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
            }else if(referenceGenome == "mm10"){
              genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
            }else{
              stop("Specified reference genome is not available!")
            }
            
            regions <- regionRanges(project)
            
            regionPFM <- pbmcapply::pbmclapply(
              sample(1:length(regions),1000),
              function(regionInd){
                regionRange <- regions[regionInd]
                regionPositions <- IRanges::tile(regionRange, width = 1)[[1]]
                contextRanges <- IRanges::resize(regionPositions, fix = "center", width = PWMRadius * 2 + 1)
                contextSeq <- Biostrings::getSeq(genome, contextRanges, as.character = T)
                ATACTrack <- ATACTracks(project)[regionInd, ]
                
                # One-hot encoding of each sequence context that appeared in the BAC
                # This is used to quantify the background sequence frequencies
                backgroundOnehot <- lapply(
                  1:length(contextSeq), 
                  function(i){onehotEncode(contextSeq[i])})
                backgroundPFM <- Reduce("+", backgroundOnehot)
                
                # Since each context is cut with different frequencies, we calculate
                # the foreground frequencies by weighing the background with cutting density 
                foregroundOnehot <- lapply(
                  1:length(contextSeq), 
                  function(i){backgroundOnehot[[i]] * ATACTrack[i]})
                foregroundPFM <- Reduce("+", foregroundOnehot)
                
                rownames(foregroundPFM) <- c("A", "C", "G", "T")
                rownames(backgroundPFM) <- c("A", "C", "G", "T")
                
                list(foregroundPFM, backgroundPFM)
              }
            )
            fgPFM <- Reduce("+", lapply(regionPFM, function(x){x[[1]]}))
            bgPFM <- Reduce("+", lapply(regionPFM, function(x){x[[2]]}))
            
            # Get background base frequencies
            background <- rowSums(bgPFM)
            background <- background / sum(background)
            
            # Get PPM 
            PPM <- t(t(fgPFM) / colSums(fgPFM))
            adjustedPPM <- PPM / background
            PWM <- log2(adjustedPPM)
            
            list(
              PPM = PPM,
              background = background,
              adjustedPPM = adjustedPPM,
              PWM = PWM
            )
          })

