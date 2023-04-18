# All the functions needed for getting the group-by-region-by-position counts tensor from a fragment file

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
                                                       returnCombined = T)
            
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
        
        # Get all Tn5 insertion sites (single base pair resolution)
        # Note: The input fragments file should be +4/-5 shifted to accommodate common practice (0-based indexing)
        # However, +4/-5 actually points to the base immediately to the left of the center of the 9bp staggered end
        # The 1bp cut position should actually by +5/-4 (for 0-based indexing, if using 1-based it should be +4/-4). 
        # Therefore we need to shift both start and end by +1
        # The function fragsToRanges already shifts start position by +1 when we specify "startsAre0based = T"
        # Therefore here we only need to further shift the end position by + 1
        cutsites <- c(resize(groupFrags, width = 1, fix = "start"),
                      IRanges::shift(resize(groupFrags, width = 1, fix = "end"), 1))
        
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
    return(countTensorAll)
  }else{
    countTensorAll <- list()
    return("Finished")
  }
  
}

