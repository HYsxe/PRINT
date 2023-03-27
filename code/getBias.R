library(hdf5r)

# Use model to predict Tn5 bias for each region
setGeneric("getRegionBias",
           function(project, # footprintingProject object
                    contextRadius = 50, # Radius of sequence context as the input to Tn5 bias model
                    nCores = 16, # Number of cores to use
                    chunkSize = 2000, # Chunk size for parallel processing of regions
                    ...) standardGeneric("getRegionBias"))

setMethod("getRegionBias", 
          signature = c(project = "footprintingProject"),
          function(project,
                   contextRadius = 50,
                   nCores = 16,
                   chunkSize = 2000) {
            
            projectDataDir <- dataDir(project)
            
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
            
            print("Extracting DNA sequence context from regions")
            regions <- regionRanges(project)
            contextLen <- 2 * contextRadius + 1
            
            # Check if the region width is correct
            if(!all(width(regions) == regionWidth(project))){
              stop("Width of provided region ranges differ from the specified width!")
            }
            
            # Extract local DNA sequence context for all positions in regions
            if(!file.exists(paste(projectDataDir, "regionSeqs.txt", sep = ""))){
              regionSeqs <- pbmcapply::pbmclapply(
                1:length(regions),
                function(regionInd){
                  range <- regions[regionInd]
                  # Get genomic sequence of the genomic region
                  # We extract an extra flanking region of length contextLen on both sides (so we can predict bias for edge positions)
                  regionSeq <- Biostrings::getSeq(genome, as.character(seqnames(range)), 
                                                  start = start(range) - contextRadius, 
                                                  width = width(range) + contextLen - 1,
                                                  as.character = T)
                  regionSeq
                },
                mc.cores = nCores
              )
              
              # Save sequence context to a file
              regionSeqs <- as.character(regionSeqs)
              write.table(regionSeqs, paste(projectDataDir, "regionSeqs.txt", sep = ""), quote = F,
                          col.names = F, row.names = F)
            }
            
            # Use pre-trained model to predict Tn5 bias for each position in each region
            # Use something like "py_install("matplotlib")" to install packages if you don't have all of them
            # Installation using reticulate can make sure we install into the same directory used by py_run_file
            projectMainDir <- mainDir(project)
            write.table(c(projectMainDir, projectDataDir, contextRadius, nCores, chunkSize), 
                        "args.txt", quote = F, col.names = F, row.names = F)
            py_run_file(paste0(mainDir(project), "code/predictBias.py"))
            
            # Load Tn5 bias predicted by our neural network
            biasPath <- paste(projectDataDir, "predBias.h5", sep = "")
            h5File <- H5File$new(biasPath, mode="r")
            predBias <- h5File[["predBias"]]
            predBias <- t(predBias[1:predBias$dims[1],])
            stopifnot(dim(predBias)[1] == length(regionRanges(project)))
            h5File$close_all()
            
            # Set the bias slot
            regionBias(project) <- predBias
            
            project
          })

# Get precomputed Tn5 bias for each region
setGeneric("getPrecomputedBias",
           function(project, # footprintingProject object
                    nCores = 4, # Number of cores to use
                    chunkSize = 100, # Chunk size for parallel processing of regions
                    ...) standardGeneric("getPrecomputedBias"))

setMethod("getPrecomputedBias", 
          signature = c(project = "footprintingProject"),
          function(project,
                   nCores = 4,
                   chunkSize = 100) {
            
            projectDataDir <- dataDir(project)
            
            # Get reference genome
            referenceGenome <- refGenome(project)
            availableGenomes <- c("ce11","danRer11","dm6","hg19","hg38","mm10","panTro6","sacCer3")
            if(referenceGenome %in% availableGenomes){
              h5_path <- paste0("../../data/shared/precomputedTn5Bias/", referenceGenome, "Tn5Bias.h5")
            }else{
              stop("Specified reference genome is not available!")
            }
            
            # Retrieve region ranges
            regions <- regionRanges(project)
            
            # First load genome-wide pre-computed Tn5 bias for each chromosome
            h5file <- H5File$new(h5_path, mode="r")
            predBiasList <- lapply(
              names(h5file),
              function(chr){
                chrData <- h5file[[chr]]
                chrData[1:chrData$dims]
              }
            )
            names(predBiasList) <- names(h5file)
            
            # Retrieve bias for each genomic region
            chunkIntervals <- getChunkInterval(regions, chunkSize = chunkSize)
            refBias <- pbmcapply::pbmclapply(
              1:length(chunkIntervals$starts),
              function(chunkInd){
                chunkStart <- chunkIntervals$starts[chunkInd]
                chunkEnd <- chunkIntervals$ends[chunkInd]
                as.data.frame(t(sapply(
                  chunkStart:chunkEnd,
                  function(regionInd){
                    chr <- as.character(seqnames(regions[regionInd]))
                    predBiasList[[chr]][start(regions[regionInd]):end(regions[regionInd])]
                  }
                )))
              },
              mc.cores = nCores
            )
            refBias <- as.matrix(data.table::rbindlist(refBias))
            
            # Set the bias slot
            regionBias(project) <- refBias
            
            project
          })
