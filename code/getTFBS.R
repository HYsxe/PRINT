filterMotifs <- function(regions,
                         genome,
                         motifs,
                         p.cutoff = 1e-5,
                         nCores = 16,
                         fdrThreshold = 0.1){
  
  # Get chromosome ranges
  if(genome == "hg38"){
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    chrRanges <- as(seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene), "GRanges")
  }else if(genome == "hg19"){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    chrRanges <- as(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene), "GRanges")
  }else if(genome == "mm10"){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    chrRanges <- as(seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene), "GRanges")
  }
  
  # Get footprint and background ranges
  print("Getting footprint and background ranges")
  regions <- regionRanges(project)
  backgroundRanges <- c(shift(regions, 5000),
                        shift(regions, -5000)) # Background is selected by sampling nearby regions
  
  # Filter out background ranges that are outside of chromosome ranges
  backgroundRanges <- subsetByOverlaps(backgroundRanges, chrRanges, type = "within")
  backgroundRanges <- backgroundRanges[-findOverlaps(backgroundRanges, regions)@from]
  backgroundRanges <- backgroundRanges[sample(1:length(backgroundRanges), length(regions))]
  
  # Scan through each foreground and background region to find motif matches
  print("Scanning through each foreground and background region to find motif matches")
  
  # Find motif matches in the CRE regions
  foregroundMatches <- pbmcapply::pbmcmapply(
    function(TF){
      as.matrix(assay(motifmatchr::matchMotifs(pwms = motifs[TF],
                                               subject = regions, 
                                               genome = genome,
                                               p.cutoff = p.cutoff)))
    },
    names(motifs),
    mc.cores = 16
  )
  
  # Find motif matches in background regions
  backgroundMatches <- pbmcapply::pbmcmapply(
    function(TF){
      as.matrix(assay(motifmatchr::matchMotifs(pwms = motifs[TF],
                                               subject = backgroundRanges, 
                                               genome = genome,
                                               p.cutoff = p.cutoff)))
    },
    names(motifs),
    mc.cores = 16
  )
  
  # Calculate motif enrichment pvals in foreground regions compared to background regions
  print("Calculating motif enrichment in foreground regions compared to background regions")
  enrichmentPvals <- pbmcapply::pbmcmapply(
    function(motifInd){
      contingency_table <- c(sum(foregroundMatches[, motifInd]),
                             sum(!foregroundMatches[, motifInd]),
                             sum(backgroundMatches[, motifInd]),
                             sum(!backgroundMatches[, motifInd]))
      contingency_table <- array(contingency_table, dim = c(2,2))
      enrichmentPval <- fisher.test(contingency_table, alternative = "greater")$p.value
      enrichmentPval
    },
    1:length(motifs),
    mc.cores = nCores
  )
  enrichmentFDRs <- p.adjust(enrichmentPvals, method = "fdr")
  enrichedInds <- which(enrichmentFDRs < fdrThreshold)
  enrichedTFs <- names(motifs)[enrichedInds]
  
  motifs[enrichedTFs]
}

# Use multi-scale footprints and some additional metadata to predict TF binding
# at a specific site
predictTFBS <- function(TFBSData, # Model input data
                        TFBSModel # Weights of the neural net model
){
  
  # First layer
  x <- TFBSData %*% TFBSModel[[1]] # Linear transform
  x <- t(t(x) + as.numeric(TFBSModel[[2]])) # Bias
  x[x < 0] <- 0 # ReLU
  
  # Second layer
  x <- x %*% TFBSModel[[3]] # Linear transform
  x <- t(t(x) + as.numeric(TFBSModel[[4]])) # Bias
  x[x < 0] <- 0 # ReLU
  
  # Third layer
  x <- x %*% TFBSModel[[5]] # Linear transform
  x <- t(t(x) + as.numeric(TFBSModel[[6]])) # Bias
  x <- 1 / (1 + exp(-x)) # Sigmoid
}

# Method to get TFBS prediction model from the footprintingProject object
setGeneric("TFBindingModel", function(x) standardGeneric("TFBindingModel"))

setMethod("TFBindingModel", "footprintingProject", function(x) x@TFBindingModel)

# Method to set TFBS prediction model in the footprintingProject object
setGeneric("TFBindingModel<-", function(x, value) standardGeneric("TFBindingModel<-"))

setMethod("TFBindingModel<-", "footprintingProject", function(x, value) {
  x@TFBindingModel <- value
  x
})

# Method to get predicted TFBS from the footprintingProject object
setGeneric("TFBS", function(x) standardGeneric("TFBS"))

setMethod("TFBS", "footprintingProject", function(x) x@TFBS)

# Method to set predicted TFBS in the footprintingProject object
setGeneric("TFBS<-", function(x, value) standardGeneric("TFBS<-"))

setMethod("TFBS<-", "footprintingProject", function(x, value) {
  x@TFBS <- value
  x
})

# Get predicted TF binding scores for a specific region
getRegionTFBS <- function(regionATAC, # Position-by-pseudobulk matrix of ATAC data for the current region
                          Tn5Bias, # Numeric vector of predicted Tn5 bias for the current region
                          region, # GRanges object for the current region
                          dispModels, # Background dispersion model for center-vs-(center + flank) ratio insertion ratio
                          TFBSModel, # Model for predicting TF binding at any motif site
                          sites = NULL, # GRanges object. Genomic ranges of motif matches in the current region
                                               # Must have $score attribute which indicates how well the motif is matched (between 0 and 1)
                                               # If NULL, the region will be divided into equally sized tiles for TF binding scoring
                          tileSize = 10, # Size of tiles if sites is NULL
                          contextRadius = 100 # Local radius of model input (in bp)
){
  
  if(!is.null(sites)){
    # Only keep motif matches within the region
    # Although we specified that the input sites should already be matches in the current region
    # this step prevents erros in case the user for got to take the overlap
    sites <- subsetByOverlaps(sites, region, type = "within")
  }else{
    sites <- GenomicRanges::slidingWindows(region, width = tileSize, step = tileSize)[[1]]
    sites$score <- 1
    sites$TF <- ""
  }

  width <- length(Tn5Bias)
  
  # These are the footprint radii used by the model. 
  # For our first version, we use 10bp, 20bp, 30bp, 50bp, 80bp, 100bp footprints
  scales <- as.character(TFBSModel$scales)
  
  # Get the pseudobulk-by-position matrix of footprint scores for the current region for each scale
  multiScaleFootprints <- lapply(
    scales,
    function(scale){
      
      smoothRadius <- as.integer(as.integer(scale) / 2)
      
      # Get the pseudobulk-by-position matrix of footprint scores for the current region with the current scale (footprint size)
      footprintScores <- t(sapply(
        1:dim(regionATAC)[2],
        function(groupInd){
          
          # Get Tn5 insertion track for one pseudobulk
          Tn5Insertion <- regionATAC[, groupInd]
          
          # Calculate footprint pval for each position
          pvals <- footprintScoring(
            Tn5Insertion = Tn5Insertion,
            Tn5Bias = Tn5Bias,
            dispersionModel = dispModels[[scale]],
            footprintRadius = as.integer(scale),
            flankRadius = as.integer(scale)
          )
          
          # Convert pvals to scores
          scores <- -log10(pvals)      
          scores[!is.finite(scores)] <- 0
          scores <- caTools::runmax(scores, 2 * smoothRadius)
          scores <- conv(scores, smoothRadius) / (2 * smoothRadius)
          scores
        }
      ))
      
      footprintScores
    }
  )
  names(multiScaleFootprints) <- scales
  
  # Calculate positions of TF sites relative to the start of the CRE region
  relativePos <- start(resize(sites, 1, fix = "center")) - start(region) + 1
  
  # Only keep sites with distance to CRE edge >= contextRadius
  siteFilter <- (relativePos > contextRadius) & (relativePos <= (width - contextRadius))
  sites <- sites[siteFilter]
  relativePos <- relativePos[siteFilter]
  
  # Go through each site and calculate predicted TF binding score for each pseudobulk
  if(length(sites) > 0){
    TFBSScores <- t(sapply(
      1:length(sites),
      function(siteInd){
        siteFootprints <- do.call("rbind", lapply(
          scales,
          function(scale){
            contextInds <- (relativePos[siteInd] - contextRadius) : (relativePos[siteInd] + contextRadius)
            scaleFootprints <- collapse::fsubset(t(multiScaleFootprints[[scale]]), 1:width %in% contextInds)
            #Flip the pattern if the site is on the negative strand
            if(as.character(strand(sites[siteInd])) == "-"){
              scaleFootprints <- scaleFootprints[rev(1:dim(scaleFootprints)[1]), ,drop = F]
            }
            scaleFootprints
          }
        ))
        siteFootprints <- t((siteFootprints - TFBSModel$footprintMean) / TFBSModel$footprintSd)
        
        # Retrieve TF site score for the current site
        TFBSData <- cbind(siteFootprints, sites[siteInd]$score)
        
        # Calculate TF binding prediction scores
        TFBSScore <- predictTFBS(TFBSData, TFBSModel)
        
      }
    ))
    
    # If there is only 1 pseudobulk, flip the axis to make sure TFBSScores is TFBS-by-pseudobulk
    if(dim(regionATAC)[2] == 1){
      TFBSScores <- t(TFBSScores)
    }
    
    colnames(TFBSScores) <- colnames(regionATAC)
    
  }else{
    TFBSScores = NULL
  }
  
  list("position" = relativePos,
       "region" = region,
       "sites" = sites,
       "TFBSScores" = TFBSScores)
  
}

# Find genomic ranges of motfi matches for all regions
setGeneric("getMotifPositions",
           function(project, # footprintingProject object
                    motifs, # PWMatrixList object of TF motifs. 
                    nCores = 16, # Number of cores to use
                    combineTFs = T # Whether to combine results for different TFs
           ) standardGeneric("getMotifPositions"))

setMethod("getMotifPositions", 
          signature = c(project = "footprintingProject"),
          function(project,
                   motifs,
                   nCores = 16,
                   combineTFs = T){
            
            if(combineTFs){
              path <- paste0(dataDir(project), "motifPositions.rds")
            }else{
              path <- paste0(dataDir(project), "motifPositionsList.rds")
            }
            
            if(file.exists(path)){
              motifPositions <- readRDS(path)
            }else{
              regions <- regionRanges(project)
              # Find motif matches across all regions
              print("Getting TF motif matches within CREs")
              motifPositions <- pbmcapply::pbmclapply(
                names(motifs),
                function(TF){
                  TFMotifPositions <- motifmatchr::matchMotifs(pwms = motifs[[TF]], 
                                                               subject = regions, 
                                                               genome = refGenome(project),
                                                               out = "positions")[[1]]
                  if(length(TFMotifPositions) > 0){
                    TFMotifPositions$TF <- TF
                    TFMotifPositions$score <- rank(TFMotifPositions$score) / length(TFMotifPositions$score)
                    TFMotifPositions <- mergeRegions(TFMotifPositions)
                    TFMotifPositions
                  }
                },
                mc.cores = nCores
              )
              names(motifPositions) <- names(motifs)
              if(combineTFs){motifPositions <- Reduce(c, motifPositions)}
              saveRDS(motifPositions, path)
            }
            
            motifPositions
          })

# Calculate predicted TF binding scores for all candidate TF binding sites in all regions of interest
setGeneric("getTFBS",
           function(project, # footprintingProject object
                    motifs = NULL, # PWMatrixList object of TF motifs. 
                    # If NULL, the region will be divided into equally sized tiles for TF binding scoring
                    contextRadius = 100, # Local radius of model input (in bp)
                    chunkSize = 2000, # Number of regions per chunk
                    chunkInds = NULL, # Whether to only run the model on a specific subset of chunks. If so, provide chunk indices here
                    innerChunkSize = 10, # Number of regions per inner chunk
                    nCores = 16 # Number of cores to use
           ) standardGeneric("getTFBS"))

setMethod("getTFBS", 
          signature = c(project = "footprintingProject"),
          function(project,
                   motifs = NULL,
                   contextRadius = 100,
                   chunkSize = 2000,
                   chunkInds = NULL,
                   innerChunkSize = 10,
                   nCores = 16) {
            
            # Directory for storing intermediate results
            tmpDir <- dataDir(project)
            
            # Determine chunk size
            if(is.null(chunkSize)){
              chunkSize <- regionChunkSize(project)
            }
            print(paste0("Using chunk size = ", chunkSize))
            
            if(length(countTensor(project)) != 0){
              chunkSize = min(chunkSize, length(countTensor(project)))
            }
            
            # Create a folder for saving intermediate results
            chunkTmpDir <- paste(tmpDir, "chunkedTFBSResults/", sep = "")
            if(!dir.exists(chunkTmpDir)){
              system(paste("mkdir -p", chunkTmpDir))
            }
            
            # Get region data we need to use later
            width <- regionWidth(project)
            seqBias <- regionBias(project)
            regions <- regionRanges(project)
            TFBSModel <- TFBindingModel(project) 
            scales <- as.character(TFBSModel$scales)
            dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
            names(dispModels) <- scales
            
            # Find motif matches across all regions
            if(!is.null(motifs)){
              motifPositions <- getMotifPositions(project, motifs)
            }else{
              motifPositions <- NULL
            }
            
            # To reduce memory usage, we chunk the region list in to smaller chunks
            cat("Chunking data ..\n")
            groupIDs <- mixedsort(groups(project))
            chunkIntervals <- getChunkInterval(regionRanges(project), chunkSize = chunkSize)
            starts <- chunkIntervals[["starts"]]
            ends <- chunkIntervals[["ends"]]
            
            # See whether the user has specified that we only run on selected chunks 
            if(is.null(chunkInds)){chunkInds <- 1:length(starts)}
            
            # Process each chunk
            for(i in chunkInds){
              
              gc()
              
              # Select regions in the current chunk
              print(paste0("Processing region chunk ", i, " out of ", length(starts), " chunks"))
              print(Sys.time())
              
              # Get ATAC insertion data for the current chunk
              chunkRegions <- starts[i]:ends[i]
              if(length(countTensor(project)) == 0){
                chunkTensorDir <- paste0(tmpDir, "chunkedCountTensor/")
                chunkCountTensor <- readRDS(paste(chunkTensorDir, "chunk_",i, ".rds", sep = ""))
              }else{
                chunkCountTensor <- countTensor(project)[chunkRegions]
              }
              names(chunkCountTensor) <- chunkRegions
              
              # Skip current chunk if result already exists
              if(file.exists(paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))){
                next
              }
              
              print(Sys.time(), "\n")
              
              # The outer for loop iterates through each chunk. Within each iteration,
              # we split the data into even smaller chunks so we can process them in parallel
              # Why do we need inner chunks? Main reason: parallel computation speeds 
              # things up when computing time for each worker is significantly larger than the time
              # required for transferring the data to each worker. We want each worker to do more
              # computation in each iteration and fewer iterations in total
              if(length(chunkRegions) <= 1){stop("Must have more than one regions")}
              if(innerChunkSize >= length(chunkRegions)){innerChunkSize <- as.integer(length(chunkRegions) / 2)}
              innerChunkIntervals <- getChunkInterval(chunkRegions, chunkSize = innerChunkSize)
              innerChunkStarts <- chunkRegions[innerChunkIntervals[["starts"]]]
              innerChunkEnds <- chunkRegions[innerChunkIntervals[["ends"]]]
              
              chunkFootprintResults <- Reduce(c,pbmcapply::pbmclapply(
                1:length(innerChunkStarts),
                function(innerChunkInd){
                  lapply(
                    innerChunkStarts[innerChunkInd] : innerChunkEnds[innerChunkInd],
                    function(regionInd){
                      
                      # Get pseudobulk-by-position matrix of ATAC data for the current region
                      regionATAC <- getRegionATAC(chunkCountTensor, as.character(regionInd), groupIDs, width[regionInd])
                      
                      # Get sites for TF binding scoring
                      if(!is.null(motifs)){
                        sites <- subsetByOverlaps(motifPositions, regions[regionInd])
                      }else{
                        sites <- NULL
                      }
                      
                      # Calculate footprint scores for each region of the current batch
                      getRegionTFBS(regionATAC = regionATAC,
                                    Tn5Bias = seqBias[regionInd,],
                                    region = regions[regionInd],
                                    dispModels = dispModels,
                                    TFBSModel = TFBSModel,
                                    sites = sites)
                      
                    }
                  )
                },
                mc.cores = nCores
              ))
              
              # Save results
              saveRDS(chunkFootprintResults,
                      paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))
            }
            
            project
          })

# Get the SummarizedExperiment object of TF binding scores
getTFBindingSE <- function(project, # footprintingProject object
                           selectedTFs = NULL, # Whether to only retrieve data for selected TFs
                           nCores = 16 # Number of cores to use
){
  
  TFBSDir <- paste0(dataDir(project), "chunkedTFBSResults/")
  TFBSChunkFiles <- gtools::mixedsort(list.files(TFBSDir))
  TFBSRanges <- NULL
  TFBSScores <- NULL
  regions <- regionRanges(project)
  if(length(regions) == 0){stop("Must load region ranges first!")}
  
  # Get the genomic ranges of TF binding sites
  TFBSRangesDf <- data.table::rbindlist(
    pbmcapply::pbmclapply(
      TFBSChunkFiles,
      function(file){
        
        # Load a new chunk of data
        TFBSChunkData <- readRDS(paste0(TFBSDir, file))
        
        # Get TF binding site genomic ranges of the current chunk
        data.table::rbindlist(lapply(
          TFBSChunkData, 
          function(x){
            sites <- x$sites
            if(length(sites) > 0){
              df <- data.frame(chr = as.character(seqnames(sites)),
                               start = start(sites),
                               end = end(sites),
                               TF = sites$TF,
                               regionInd = findOverlaps(x$region, regions, type = "equal")@to)
              if(!is.null(selectedTFs)){
                df <- df[df$TF %in% selectedTFs,]
              }
            }else{
              df <- NULL
            }
            df
          }))
        
      },
      mc.cores = 16
    )
  )
  TFBSRanges <- GRanges(seqnames = as.character(TFBSRangesDf$chr),
                        ranges = IRanges(start = TFBSRangesDf$start,
                                         end = TFBSRangesDf$end))
  TFBSRanges$regionInd <- TFBSRangesDf$regionInd
  TFBSRanges$TF <- TFBSRangesDf$TF
  
  # Get TFBS-by-pseudobulk matrix of TF binding scores
  TFBSScores <- data.table::rbindlist(
    pbmcapply::pbmclapply(
      TFBSChunkFiles,
      function(file){
        
        # Load a new chunk of data
        TFBSChunkData <- readRDS(paste0(TFBSDir, file))
        
        # Get TF binding site genomic ranges of the current chunk
        data.table::rbindlist(lapply(
          TFBSChunkData, 
          function(x){
            sites <- x$sites
            if(length(sites) > 0){
              df <- as.data.frame(x$TFBSScores)
              if(!is.null(selectedTFs)){
                df <- df[sites$TF %in% selectedTFs,]
              }
            }else{
              df <- NULL
            }
            df
          }))
        
      },
      mc.cores = 16
    )
  )
  
  TFBSScores <- as.matrix(TFBSScores)
  TFBindingSE <- SummarizedExperiment(assay = list(TFBSScores), rowRanges = TFBSRanges)
  
  if(length(groups(project)) == dim(TFBindingSE)[2]){
    colnames(TFBindingSE) <- groups(project)
  }
  
  TFBindingSE
}

# Get the TF binding score matrix for a specific region
getRegionTFBSMatrix <- function(project, # footprintingProject object
                                regionInd, # Index of the region
                                groupIDs = NULL # Whether to use only selected pseudobulks. If so, provide pseudobulk indices here
){
  
  if(is.null(groupIDs)){
    groupIDs <- groups(project)
  }
  
  ret <- getCountData(project, regionInd)
  countData <- ret[["countData"]]
  adjustedRegionInd <- ret[["regionInd"]]
  
  # Get position-by-pseudobulk ATAC insertion matrix
  regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, regionWidth(project)[regionInd])
  
  TFBSModel <- TFBindingModel(project)
  scales <- TFBSModel$scales # Scales are footprint radii we use for multi-scale footprinting
  dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
  names(dispModels) <- scales
  region <- regionRanges(project)[regionInd]
  
  # Calculate footprint scores for each region of the current batch
  regionTFBS <- getRegionTFBS(regionATAC = regionATAC,
                              Tn5Bias = regionBias(project)[regionInd,],
                              region = region,
                              dispModels = dispModels,
                              TFBSModel = TFBindingModel(project),
                              tileSize = 10)
  
  featureMatrix <- array(0, dim = dim(regionATAC))
  rownames(featureMatrix) <- rep("", dim(featureMatrix)[1])
  for(ind in 1:length(regionTFBS$position)){
    siteStart <- start(regionTFBS$sites[ind]) - start(region) + 1
    siteEnd <- end(regionTFBS$sites[ind]) - start(region) + 1
    for(pos in siteStart:siteEnd){
      featureMatrix[pos, ] <- pmax(featureMatrix[pos, ],regionTFBS$TFBSScores[ind,])
    }
  }
  
  featureMatrix
}

# Prepare training data for TFBS prediction model
getTFBSTrainingData <- function(multikernelFootprints, # Pre-computed multi-scale footprint tracks for all regions. 
                                # See TFBSTrainingData.R for details on how this is generated
                                motifMatches, # Genomic ranges of matched motifs for each TF. Should be a named list of GRanges objects
                                # See TFBSTrainingData.R for details on how this is generated
                                TFChIPRanges, # Genomic ranges of TF ChIP peaks. Should be a named list of GRanges objects
                                regions, # All genomic regions used for footprinting
                                contextRadius = 100, # Local radius of model input (in bp)
                                percentBoundThreshold = 0.2
){
  
  # For each motif matched site, get local multi-kernel footprints and bound/unbound labels
  motifFootprints <- NULL
  metadata <- NULL
  for(TF in names(motifMatches)){
    
    print(paste("Getting TFBS for ", TF))
    
    # Find bound and unbound motif positions
    motifRanges <- motifMatches[[TF]]
    motifRanges <- mergeRegions(motifRanges)
    motifRegionOv <- findOverlaps(motifRanges, regions)
    
    # Skip TFs with a low percentage of motifs overlapping with ChIP
    percentBound <- length(subsetByOverlaps(motifRanges, TFChIPRanges[[TF]])) / length(motifRanges)
    if(percentBound < percentBoundThreshold){next}
    
    # For each motif site, extract multiscale footprints and bound/unbound labels
    width <- unique(width(regions))
    footprintContext <- pbmcapply::pbmclapply(
      1:length(motifRegionOv),
      function(i){
        
        # Get the current region-motif pair
        ovMotifInd <- motifRegionOv@from[i] 
        ovRegionInd <- motifRegionOv@to[i]
        ovMotif <- motifRanges[ovMotifInd]
        ovRegion <- regions[ovRegionInd]
        
        # Get motif match score
        matchScore <- ovMotif$score
        
        # Find out position of the motif relative to start site of the region
        relativePos <- start(resize(ovMotif, 1, fix = "center")) - start(ovRegion) + 1
        contextInds <- (relativePos - contextRadius) : (relativePos + contextRadius)
        
        if((relativePos > contextRadius) & (relativePos <= width - contextRadius)){
          
          # Retrieve multi-scale footprints for the current motif site
          footprints <- lapply(
            names(multikernelFootprints),
            function(kernelSize){
              kernelFootprints <- multikernelFootprints[[as.character(kernelSize)]][ovRegionInd, contextInds]
              # Flip the pattern if the motif is on the negative strand
              if(as.character(strand(ovMotif)) == "-"){kernelFootprints <- rev(kernelFootprints)}
              kernelFootprints
            }
          )
          footprints <- Reduce(c, footprints)
          
          # Binarize the motif into bound/unbound based on overlap with ChIP
          bound <- as.integer(length(findOverlaps(ovMotif, TFChIPRanges[[TF]])) > 0)
          
          list(footprints, matchScore, bound, as.character(ovMotif))
        }else{
          NULL
        }
      },
      mc.cores = 16
    )
    footprintContext <- footprintContext[!sapply(footprintContext, is.null)]
    TFMatchScores <- sapply(footprintContext, function(x){x[[2]]})
    motifFootprints <- rbind(motifFootprints, t(sapply(footprintContext, function(x){x[[1]]})))
    metadata <- rbind(metadata, data.frame(sapply(footprintContext, function(x){x[[3]]}), 
                                           TFMatchScores, TF, 
                                           sapply(footprintContext, function(x){x[[4]]})))
  }
  
  # Cap infinite values
  motifFootprints[!is.finite(motifFootprints)] <- 10
  
  colnames(metadata) <- c("bound", "motifMatchScore", "TF", "range")
  
  list("motifFootprints" = motifFootprints, 
       "metadata" = metadata)
}


# Load TF binding / habitation prediction model
loadTFBSModel <- function(
  h5Path = "../../data/TFBSPrediction/TFBS_model.h5"
){
  
  TFBSModel <- keras::load_model_hdf5(h5Path)
  TFBSModel <- keras::get_weights(TFBSModel)
  h5file <- H5File$new(h5Path, mode="r")
  footprintMean <- h5file[["footprint_mean"]]
  footprintMean <- footprintMean[1:footprintMean$dims]
  footprintSd <- h5file[["footprint_sd"]]
  footprintSd <- footprintSd[1:footprintSd$dims]
  h5file$close_all()
  TFBSModel$footprintMean <- footprintMean
  TFBSModel$footprintSd <- footprintSd
  TFBSModel$scales <- c(10,20,30,50,80,100)
  TFBSModel
  
}

