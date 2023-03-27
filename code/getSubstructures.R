# Iterate through each region we footprinted and segment the CRE regions into substructures
setGeneric("getSubstructureSE",
           function(project, # footprintingProject object
                    chunkSize = 2000, # Chunk size for parallel processing of regions
                    nCores = 16, # Number of cores to use
                    windowSize = 10 # Window size for the segmentation algorithm
                    ) standardGeneric("getSubstructureSE"))

setMethod("getSubstructureSE", 
          signature = c(project = "footprintingProject"),
          function(project, 
                   chunkSize = 2000,
                   nCores = 16,
                   windowSize = 10){
            
            chunkDir <- paste0(dataDir(project), "chunkedTFBSResults/")
            chunkFiles <- gtools::mixedsort(list.files(chunkDir))
            if(length(chunkFiles) == 0){stop("Need to compute TFBS scores first!")}
            
            substructureScores <- list()
            substructureRanges <- list()
            for(chunkInd in 1:length(chunkFiles)){
              
              # Load precomputed TFBS data for the current chunk of regions
              print(paste("Getting substructures in chunk", chunkInd))
              chunkTFBSPath <- paste0(chunkDir, chunkFiles[chunkInd])
              if(!file.exists(chunkTFBSPath)){
                stop(paste0("For chunk ", chunkInd, ", need to run getTFBS and compute TFBS scores first!"))
              }
              chunkTFBS <- readRDS(chunkTFBSPath)
              
              chunkSubstructures <- pbmcapply::pbmclapply(
                chunkTFBS,
                function(regionTFBS){
                  
                  if(length(regionTFBS$position) == 0){
                    list("substructureRanges" = NULL,
                         "substructureScores" = NULL)
                  }else if(length(regionTFBS$position) == 1){
                    list("substructureRanges" = paste0(as.character(seqnames(regionTFBS$region)), ":",
                                                start(regionTFBS$sites), "-",
                                                end(regionTFBS$sites)),
                         "substructureScores" = as.data.frame(t(regionTFBS$TFBSScores)))
                  }else{
                    
                    # Calculate correlation of TF binding scores among matched motif positions
                    TFBSScores <- regionTFBS$TFBSScores
                    corMat <- cor(t(TFBSScores))
                    
                    # Segment the enhancer region
                    segmentation <- CRESegmentation(corMat, windowSize = windowSize)
                    boundaries <- segmentation[["Boundaries"]]
                    BI <- segmentation[["boundaryIndex"]]
                    
                    # Get substructure-by-pseudobulk matrix of substructure activity
                    substructureRanges <- NULL
                    substructureScoreMat <- NULL
                    for(i in 1:(length(boundaries) - 1)){
                      if(i == 1){
                        substructureStart <- 1
                      }else{
                        substructureStart <- boundaries[i] + 1
                      }
                      substructureEnd <- boundaries[i + 1]
                      substructureScores <- colMeans(TFBSScores[substructureStart:substructureEnd, ,drop = F])
                      substructureScoreMat <- rbind(substructureScoreMat, substructureScores)
                      substructureRange <- paste0(as.character(seqnames(regionTFBS$region)), ":",
                                           start(regionTFBS$sites[substructureStart]), "-",
                                           end(regionTFBS$sites[substructureEnd]))
                      substructureRanges <- c(substructureRanges, substructureRange)
                    }
                    
                    list("substructureRanges" = substructureRanges,
                         "substructureScores" = as.data.frame(substructureScoreMat))
                  }
                },
                mc.cores = 16
              )
              
              # Combine results for the current chunk
              chunkSubstructureRanges <- Reduce(c, lapply(chunkSubstructures, function(x){x$substructureRanges}))
              chunkSubstructureScores <- data.table::rbindlist(lapply(chunkSubstructures, function(x){x$substructureScores}))
              
              substructureScores[[chunkInd]] <- chunkSubstructureScores
              substructureRanges[[chunkInd]] <- chunkSubstructureRanges
            }
            
            substructureScores <- data.table::rbindlist(substructureScores)
            substructureRanges <- Reduce(c, substructureRanges)
            
            substructureSE <- SummarizedExperiment(assays = list("substructures" = as.matrix(substructureScores)),
                                            rowRanges = GRanges(substructureRanges))
            
            substructureSE
          }
        )

# Given a TFBS-by-TFBS matrix for all TFBS within a CRE, segment the CRE into substructures
CRESegmentation <- function(data, # TFBS-by-TFBS correlation matrix 
                            windowSize = 10 # Local window radius for substructure boundary detection
){
  
  width <- dim(data)[1]
  r <- as.integer(windowSize/2)
  
  # Calculate boundary index
  boundaryIndex <- sapply(
    1:(width - 1),
    function(x){
      scoreA <- mean(data[max(0, x - windowSize):x, max(0, x - windowSize):x])
      scoreB <- mean(data[max(0, x - windowSize):x, (x + 1):min(width, x + windowSize)])
      scoreC <- mean(data[(x + 1):min(width, x + windowSize), (x + 1):min(width, x + windowSize)])
      max(scoreA, scoreC) - scoreB
    }
  )
  
  # Smooth boundary index
  boundaryIndex <- sapply(
    1:(width - 1),
    function(x){
      mean(boundaryIndex[max(0, x - r):min(width - 1, x + r)])
    }
  )
  
  # Detect boundaries
  boundaryIndex <- sapply(
    1:(width - 1),
    function(x){
      boundaryIndex[x] - max(boundaryIndex[max(0, x - r):min(width - 1, x + r)]) 
    }
  )
  boundaries <- which(abs(boundaryIndex - 0) < 1e-5)
  boundaries <- boundaries[(boundaries > r) & (boundaries <= (width - r))]
  
  list("boundaryIndex" = boundaryIndex, "Boundaries" = c(1, boundaries, width))
  
}

