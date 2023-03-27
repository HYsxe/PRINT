# Method to get footprint tracks from the footprintingProject object
setGeneric("footprintTracks", function(x, ...) standardGeneric("footprintTracks"))

setMethod("footprintTracks", "footprintingProject", function(x, mode) x@footprintTracks[[mode]])

# Method to set footprint tracks from the footprintingProject object
setGeneric("footprintTracks<-", function(x, value, ...) standardGeneric("footprintTracks<-"))

setMethod("footprintTracks<-", "footprintingProject", function(x, value, mode) {
  x@footprintTracks[[mode]] <- value
  x
})

# Tracking footprint dynamics across pseudo-time
setGeneric("footprintTracking",
           function(project, # footprintingProject object
                    lineageGroups, # Indices of pseudo-bulks in the lineage we want to run tracking on
                    threshold = 2, # Identify position of object by finding footprint signal summit above this threshold value
                    pseudoTimeWindowSize = 10, # Size of the sliding window across pseudo-time. 
                    chunkSize = NULL, # Number of regions per chunk
                    footprintRadius = 50, # Radius of the footprint region
                    nCores = 12, # Number of cores to use
                    ...) standardGeneric("footprintTracking"))

setMethod("footprintTracking", 
          signature = c(project = "footprintingProject"),
          function(project,
                   lineageGroups,
                   threshold = 1,
                   pseudoTimeWindowSize = 10, 
                   chunkSize = NULL,
                   footprintRadius = 50,
                   nCores = 12){
            
            # Determine chunk size
            if(is.null(chunkSize)){
              chunkSize <- regionChunkSize(project)
            }
            print(paste0("Using chunk size = ", chunkSize))
            
            # Create a folder for saving intermediate results
            chunkTmpDir <- paste(dataDir(project), "chunkedTrackingResults/", footprintRadius, "bp/", sep = "")
            if(!dir.exists(chunkTmpDir)){
              system(paste("mkdir -p", chunkTmpDir))
            }
            
            # Get the ATAC count tensor from project
            projectCountTensor <- countTensor(project)
            
            # Get sequence bias
            seqBias <- regionBias(project)
            rownames(seqBias) <- 1:length(regionRanges(project))
            width <- regionWidth(project)
            dispersionModel <- dispModel(project, as.character(footprintRadius))
            smoothRadius <- as.integer(footprintRadius / 2)
            
            # To reduce memory usage, we chunk the region list in to smaller chunks
            chunkIntervals <- getChunkInterval(regionRanges(project))
            starts <- chunkIntervals[["starts"]]
            ends <- chunkIntervals[["ends"]]
            
            # Remove ties from pseudo-time values and re-order pseudo-bulks by pseudo-time
            lineagePseudoTime <- groupPseudoTime(project)[lineageGroups]
            lineagePseudoTime <- rank(lineagePseudoTime, ties.method = "random")
            lineageGroups <- lineageGroups[order(lineagePseudoTime)]
            
            # Process each chunk
            for(i in 1:length(starts)){
              
              # Select regions in the current chunk
              print(paste("Processing region chunk",i))
              print(Sys.time())
              
              # Get ATAC insertion data for the current chunk
              chunkRegions <- starts[i]:ends[i]
              if(length(projectCountTensor) == 0){
                chunkTensorDir <- paste0(dataDir(project), "chunkedCountTensor/")
                chunkCountTensor <- readRDS(paste(chunkTensorDir, "chunk_",i, ".rds", sep = ""))
              }else{
                chunkCountTensor <- projectCountTensor[chunkRegions]
              }
              names(chunkCountTensor) <- chunkRegions
              
              chunkBias <- seqBias[chunkRegions,]
              
              # Skip current chunk if result already exists
              if(file.exists(paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))){
                next
              }
              
              # Tracking footprint positions in chunks
              print("Tracking footprint positions in chunks")
              cluster <- prep_cluster(chunkSize, n_cores = nCores)
              opts <- cluster[["opts"]]
              cl <- cluster[["cl"]]
              chunkTrackList <- 
                foreach(regionInd = chunkRegions,
                        .options.snow = opts, 
                        .packages = c("dplyr","Matrix"),
                        .export = c("regionBias", "getRegionATAC","regionWidth",
                                    "footprintScoring", "findSummits", "conv",
                                    "getRegionFootprintTracks", "countTensor",
                                    "footprintWindowSum", "predictDispersion")) %dopar%   {
                                      
                                      # Get Tn5 bias for every bp in the current region
                                      Tn5Bias <- chunkBias[as.character(regionInd), ]
                                      
                                      # Create a sliding window along pseudo-time
                                      peudoTimeWindows <- lapply(1:(length(lineageGroups) - pseudoTimeWindowSize + 1),
                                                                 function(i){i:(i + pseudoTimeWindowSize - 1)})
                                      
                                      # Get ATAC signal tracks for each pseudo-bulk
                                      regionATAC <- getRegionATAC(chunkCountTensor, as.character(regionInd), lineageGroups, width[regionInd])
                                      
                                      # Aggregate data from all pseudo-bulks and identify rough positioning of DNA-bound objects
                                      # These positions will be used to determine a local window for tracking the movement of each object later
                                      aggregateATAC <- rowSums(regionATAC)
                                      aggregateFootprintPvals <- footprintScoring(aggregateATAC, 
                                                                                  Tn5Bias,
                                                                                  dispersionModel,
                                                                                  footprintRadius = footprintRadius,
                                                                                  flankRadius = footprintRadius)
                                      aggregateFootprintScores <- -log10(aggregateFootprintPvals)
                                      aggregateFootprintScores <- caTools::runmax(aggregateFootprintScores, 2 * smoothRadius)
                                      aggregateFootprintScores <- conv(aggregateFootprintScores, smoothRadius) / (2 * smoothRadius)
                                      aggregateSummits <- findSummits(aggregateFootprintScores, 
                                                                      r = footprintRadius,
                                                                      threshold = threshold)
                                      
                                      # For every pseudo-time window, we pool the data of pseudo-bulks in the same window
                                      # Then we calculate the position-by-pseudoTimeWindow matrix of footprint scores
                                      regionFootprintScores <- sapply(
                                        peudoTimeWindows,
                                        function(peudoTimeWindow){
                                          regionWindowATAC <- rowSums(regionATAC[,peudoTimeWindow])
                                          pvals <- footprintScoring(regionWindowATAC, 
                                                                    Tn5Bias,
                                                                    dispersionModel,
                                                                    footprintRadius = footprintRadius,
                                                                    flankRadius = footprintRadius)
                                          scores <- -log10(pvals)
                                          scores <- caTools::runmax(scores, 2 * smoothRadius)
                                          scores <- conv(scores, smoothRadius) / (2 * smoothRadius)
                                        }
                                      )
                                      
                                      # For each footprint summit position we found in the aggregate data,
                                      # we go through the sliding pseudotime windows to locate the associated object in
                                      # each pseudotime window. We end up with a track of positions for each object across pseudo-time
                                      regionFootprintTracks <- getRegionFootprintTracks(aggregateSummits, 
                                                                                    regionFootprintScores, 
                                                                                    footprintRadius,
                                                                                    width[regionInd])
                                      
                                      return(regionFootprintTracks)
                                    }
              stopCluster(cl)
              print("Finished!")
              
              # Save results
              saveRDS(chunkTrackList,
                      paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))
            }
            
            # Integrate results for all chunks
            chunkFiles <- gtools::mixedsort(list.files(chunkTmpDir))
            chunkFiles <- chunkFiles[stringr::str_detect(chunkFiles, "chunk_")]
            chunkResults <- lapply(chunkFiles, function(f){readRDS(paste(chunkTmpDir, f, sep = ""))})
            names(chunkResults) <- sapply(chunkFiles, function(f){strsplit(f, "\\.")[[1]][1]})
            trackingResults <- Reduce(c, chunkResults)
            rm(chunkResults)
            names(trackingResults) <- as.character(regionRanges(project))
            
            # Assign results to the corresponding slot in the project object
            trackingResults
            
          }
)

# In a specific region, find the rough position of each DNA-binding object using aggregate footprint signal summits.
# Then within a certain window around each summit, go through pseudo-time to identify the precise position of the object
# at each time point by finding local max position of signal.
getRegionFootprintTracks <- function(aggregateSummits, # Signal summits obtained by aggregating all data across pseudo-time
                                   regionFootprintScores,  # Position-by-pseudoTimeWindow matrix of footprint scores
                                   footprintRadius, # Radius of the footprint region. Used to determine size of local search window
                                   width # Width of the whole CRE region
                                   ){
  regionFootprintTracks <- lapply(
    aggregateSummits,
    function(summit){
      # For each summit, extract footprint scores within a local window
      summitNeighborhood <- regionFootprintScores[max(summit - footprintRadius, 1):min(summit + footprintRadius, width), ]
      t(sapply(
        1:dim(regionFootprintScores)[2],
        function(ptime){
          # For each pseudo-time point, find position of max signal within the local window
          # and use it as the position of the object at that moment
          ptimeSummitPosition <- which(summitNeighborhood[, ptime] == max(summitNeighborhood[, ptime])) 
          if(length(ptimeSummitPosition) > 1){# If there is more than one local maximum, we can't determine the position of the nucleosome
            ptimeSummitPosition <- sample(ptimeSummitPosition, 1)
          }
          ptimeSummitScore <- summitNeighborhood[ptimeSummitPosition, ptime]
          c(ptimeSummitPosition + summit - footprintRadius - 1, ptimeSummitScore)
        }
      ))
    }
  )
}
