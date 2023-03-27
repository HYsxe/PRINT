# Method to get the region-by-pseudobulk ATAC count matrix
setGeneric("getGroupATAC",
           function(project, # footprintingProject object
                    nCores = 16 # Number of cores to use
                    ) standardGeneric("getGroupATAC"))

setMethod("getGroupATAC", 
          signature = c(project = "footprintingProject"),
          function(project,
                   nCores = 16
          ) {
            
            chunkSize <- regionChunkSize(project)
            projectCountTensor <- countTensor(project)
            
            # To reduce memory usage, we chunk the region list in to smaller chunks
            cat("Chunking data ..\n")
            groupIDs <- mixedsort(groups(project))
            chunkIntervals <- getChunkInterval(regionRanges(project), chunkSize = chunkSize)
            starts <- chunkIntervals[["starts"]]
            ends <- chunkIntervals[["ends"]]
            
            # Directory for storing intermediate results
            tmpDir <- dataDir(project)
            
            chunkRegionCounts <- list()
            for(i in 1:length(starts)){
              
              # Select regions in the current chunk
              print(paste0("Processing region chunk ", i, " out of ", length(starts), " chunks"))
              print(Sys.time())
              
              # Get ATAC insertion data for the current chunk
              chunkRegions <- starts[i]:ends[i]
              if(length(projectCountTensor) == 0){
                chunkTensorDir <- paste0(tmpDir, "chunkedCountTensor/")
                chunkCountTensor <- readRDS(paste(chunkTensorDir, "chunk_",i, ".rds", sep = ""))
              }else{
                chunkCountTensor <- projectCountTensor[chunkRegions]
              }
              names(chunkCountTensor) <- chunkRegions
              
              print(Sys.time(), "\n")
              
              # For the current chunk of regions, get region-by-pseudobulk ATAC count matrix
              chunkRegionCounts[[i]] <- as.data.frame(t(pbmcapply::pbmcmapply(
                function(regionCountTensor){
                  regionGroupCounts <- regionCountTensor %>% group_by(group) %>% summarise(count = sum(count))
                  regionCounts <- array(0, length(groupIDs))
                  if(class(groupIDs) == "character"){
                    names(regionCounts) <- groupIDs
                  }
                  regionCounts[regionGroupCounts$group] <- regionGroupCounts$count
                  regionCounts
                },
                chunkCountTensor,
                mc.cores = nCores
              )))
              
            }
            
            pseudobulkATAC <- data.table::rbindlist(chunkRegionCounts) 
            
            as.matrix(pseudobulkATAC)

            })

# Get gene-by-pseudobulk matrix
getGroupRNA <- function(RNAMatrix, # Single cell gene-by-cell matrix of RNA levels. Raw counts are preferred
                        barcodeGroups, # Data.frame specifying membership of barcodes in pseudobulks. Two columns should be barcodes and group membership 
                        nCores = 8 # Number of cores to use
                        ){
  
  # Rename columns in barcodeGrouping
  colnames(barcodeGroups) <- c("barcode", "group")
  
  if(!all(barcodeGroups$barcode %in% colnames(RNAMatrix))){
    stop("Not all pseudobulk barcodes are in the RNA matrix!")
  }
  
  groupIDs <- gtools::mixedsort(unique(barcodeGroups$group))
  
  # Get RNA data for each pseudo-bulk
  groupedRNA <- pbmcapply::pbmcmapply(
    
    function(groupID){
      groupBarcodes <- barcodeGroups$barcode[barcodeGroups$group %in% groupID]
      rowSums(RNAMatrix[,groupBarcodes])
    },
    groupIDs,
    mc.cores = nCores
  )
  
  groupedRNA
  
}


