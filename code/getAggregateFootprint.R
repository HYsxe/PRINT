library(ggplot2)

# Get region-by-position matrix of Tn5 insertion counts (collapse the pseudobulk dimension)
setGeneric("getATACTracks",
           function(project, # footprintingProject object
                    nCores = 16, # Number of cores to use
                    ...) standardGeneric("getATACTracks"))

setMethod("getATACTracks", 
          signature = c(project = "footprintingProject"),
          function(project, 
                   nCores = 16){
            
            chunkSize <- regionChunkSize(project)
            width <- regionWidth(project)
            
            # Get ATAC track for each region
            # To reduce memory usage, we chunk the region list in to smaller chunks
            cat("Chunking data ..\n")
            chunkIntervals <- getChunkInterval(regionRanges(project), chunkSize = chunkSize)
            starts <- chunkIntervals[["starts"]]
            ends <- chunkIntervals[["ends"]]
            
            # Process each chunk
            regionATACTracks <- NULL
            for(i in 1:length(starts)){
              
              # Select regions in the current chunk
              print(paste0("Processing region chunk ", i, " out of ", length(starts), " chunks"))
              print(Sys.time())
              
              # Get ATAC insertion data for the current chunk
              chunkRegions <- starts[i]:ends[i]
              if(length(countTensor(project)) == 0){
                chunkTensorDir <- paste0(dataDir(project), "chunkedCountTensor/")
                filePath <- paste(chunkTensorDir, "chunk_",i, ".rds", sep = "")
                if(!file.exists(filePath)){
                  stop("Need to first run getCountTensor and generate count tensor!")
                }
                chunkCountTensor <- readRDS(paste(chunkTensorDir, "chunk_",i, ".rds", sep = ""))
              }else{
                chunkCountTensor <- projectCountTensor[chunkRegions]
              }
              names(chunkCountTensor) <- chunkRegions
              
              # Generate region-by-position matrix of Tn5 insertion counts
              chunkTracks <- t(pbmcapply::pbmcmapply(
                function(regionInd){
                  regionCounts <- chunkCountTensor[[as.character(regionInd)]]
                  regionCounts <- regionCounts %>% group_by(position) %>% summarise(insertion = sum(count))
                  ATACTrack <- rep(0, width[regionInd])
                  ATACTrack[regionCounts$position] <- regionCounts$insertion
                  ATACTrack
                },
                chunkRegions,
                mc.cores = nCores
              ))
              
              regionATACTracks <- rbind(regionATACTracks, chunkTracks)
            }
            
            ATACTracks(project) <- regionATACTracks
            
            project
            
          })

# Get aggregate footprint around a list of genomic sites
setGeneric("getAggregateFootprint",
           function(project, # footprintingProject object
                    sites, # GRanges object specifying the sites to aggregate
                    nTopSites = NULL, # Whether to only aggregate n top accessible sites
                    nCores = 16, # Number of cores to use
                    ...) standardGeneric("getAggregateFootprint"))

setMethod("getAggregateFootprint", 
          signature = c(project = "footprintingProject"),
          function(project, 
                   sites,
                   nTopSites = NULL,
                   nCores = 16,
                   chunkSize = 2000) {
            
            # Get region ranges and width 
            regions <- regionRanges(project)
            width <- regionWidth(project)
            maxWidth <- max(width)
            genome <- refGenome(project)
            
            # Get the region-by-position ATAC matrix
            if(dim(ATACTracks(project))[1] == 0){
              project <- getATACTracks(project)
            }
            regionATACTracks <- ATACTracks(project)
            
            # Get width of the sites
            siteWidth <- width(sites)
            
            # Pair sites and regions together
            sitePeakOv <- findOverlaps(sites, regions, type = "within")
            sitePeakOv <- data.frame(from = sitePeakOv@from,
                                      to = sitePeakOv@to)
            
            # If specified, select top accessible sites
            # Otherwise, select all sites
            if(!is.null(nTopSites)){
              regionAccessibility <- rowSums(regionATACTracks)
              selectedOvInds <- which(rank(-regionAccessibility[sitePeakOv$to]) <= nTopSites)
            }else{
              selectedOvInds <- 1:length(sitePeakOv$from)
            }
            
            # Get predicted Tn5 bias for each region
            predBias <- regionBias(project)
            
            # For each selected site, shift ATAC signal so it's centered around the site
            centeredATAC <- pbmcapply::pbmclapply(
              selectedOvInds,
              function(ovInd){
                
                ov <- sitePeakOv[ovInd,]
                siteInd <- ov$from
                regionInd <- ov$to
                siteStrand <- as.character(strand(sites)[siteInd])
                
                # Get region and site centers
                siteCenter <- resize(sites[ov$from], 1, fix = "center")
                regionCenter <- resize(regions[regionInd], 1, fix = "center")
                
                # If the site is on the negative strand, adjust the center position
                # For example, for a region of 1000bp, the center for the forward direction is 500
                # while the center for the opposite direction should be 501 (500 looking from the other side)
                # No need to adjust if width of region or site is an odd number
                if(siteStrand == "-"){
                  if(width[regionInd] %% 2 == 0){
                    regionCenter <- GenomicRanges::shift(regionCenter, 1)
                  }
                  if(siteWidth[siteInd] %% 2 == 0){
                    siteCenter <- GenomicRanges::shift(siteCenter, 1)
                  }
                }
                
                # Calculate site offset relative to region center
                if(siteStrand == "-"){
                  centerOffset <- start(regionCenter) - start(siteCenter)
                }else{
                  centerOffset <- start(siteCenter) - start(regionCenter)
                }
                
                # Get Tn5 sequence bias track for the region
                Tn5Bias <- predBias[regionInd, ]
                
                # Get ATAC track for the current region
                ATACTrack <- regionATACTracks[regionInd,]
                
                # Reverse the relevant track if the site is on the negative strand
                if(siteStrand == "-"){
                  Tn5Bias <- rev(Tn5Bias)
                  ATACTrack <- rev(ATACTrack)
                }
                
                # Shift the ATAC data according to the offset
                foregroundTrack <- rep(0, 2 * maxWidth)
                backgroundTrack <- rep(0, 2 * maxWidth)
                correctedTrack <- rep(0, 2 * maxWidth)
                biasTrack <- rep(0, 2 * maxWidth)
                shiftedStart <- maxWidth - round(width[regionInd] / 2) + 1 - centerOffset
                shiftedEnd <- shiftedStart + width[regionInd] - 1
                foregroundTrack[shiftedStart:shiftedEnd] <- ATACTrack
                backgroundTrack[shiftedStart:shiftedEnd] <- 1
                biasTrack[shiftedStart:shiftedEnd] <- Tn5Bias
                correctedTrack[shiftedStart:shiftedEnd] <- ATACTrack / pmax(Tn5Bias, 0.01) # Prevent division by zero
                
                list(foregroundTrack, 
                     correctedTrack,
                     backgroundTrack, 
                     biasTrack)
              },
              mc.cores = 8
            )
            
            # Aggregate the shifted tracks for all sites
            foregroundTrack <- Reduce("+", lapply(centeredATAC, function(x){x[[1]]}))
            correctedTrack <- Reduce("+", lapply(centeredATAC, function(x){x[[2]]}))
            backgroundTrack <- Reduce("+", lapply(centeredATAC, function(x){x[[3]]}))
            backgroundTrack <- pmax(backgroundTrack, 1) # Prevent division by zero
            
            # Normalize the values of biasTrack to our usual range. This is because the background dispersion model would fail 
            # if the bias values are too far outside the usual range
            # We don't normalize foregroundTrack because after aggregation we have in fact many reads aggregated. If we normalize
            # we get the average cutting which is a very low value. Then our dispersion model would predict a very large dispersion
            # which is incorrect.
            biasTrack <- Reduce("+", lapply(centeredATAC, function(x){x[[4]]})) / backgroundTrack
            
            # Return results
            list("correctedATAC" = correctedTrack,
                 "uncorrectedATAC" = foregroundTrack,
                 "background" = backgroundTrack,
                 "Tn5Bias" = biasTrack)
          })
