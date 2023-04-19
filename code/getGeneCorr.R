splitAndFetch <- function (vec, delim, part) 
{
  if (length(part) == 1) {
    sapply(strsplit(as.character(vec), delim, fixed = TRUE), 
           "[[", part)
  }
  else {
    sapply(strsplit(as.character(vec), delim, fixed = TRUE), 
           function(x) paste(x[part], collapse = delim))
  }
}

geneCorr <- function(ATAC, # ATAC feature (e.g., Footprint or CRE) RangedSummarizedExperiment
                     RNA, # Gene-by-pseudobulk (or -by-single cell) matrix of RNA. Should be log-normalized
                     genome, # Must be one of "hg19", "mm10", or "hg38"
                     geneList = NULL, # 2 or more valid gene symbols
                     windowPadSize = 50000, # Window size for identifying ATAC feature-RNA pairs
                     nCores = 12, # Number of cores to use
                     nBg = 100, # Number of background sampling iterations for each foreground
                     p.cut = NULL, # Optional, if specified, will only return sig hits
                     multimapping = T, # Whether to allow each feature to be mapped to multiple genes,
                     metric = "pearson", # Metric for correlation
                     numBgPairs = 100000, # Total number of precomputed background CRE/footprint-RNA pairs
                     pairsPerChunk = 500 # Number of CRE/footprint-RNA pairs to correlate per chunk
) {
  
  stopifnot(inherits(ATAC,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNA,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC),ncol(RNA)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  # Function needs rownames for both matrices or gives error
  rownames(ATAC) <- paste0("Feature",1:nrow(ATAC))
  
  if(is.null(rownames(RNA)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for ATAC features with 0 signal
  if(any(Matrix::rowSums(assay(ATAC)) == 0)){
    message("Features with 0 signal across cells exist ..")
    message("Removing these features prior to running correlations ..")
    featuresToKeep <- Matrix::rowSums(assay(ATAC)) != 0
    ATAC <- ATAC[featuresToKeep,] # Subset ranges
    message("Important: feature indices in returned gene-feature maps are relative to original input SE")
  }
  
  ATACmat <- assay(ATAC) # Rownames preserved
  
  # Normalize feature signal
  ATACmat <- centerCounts(ATACmat)
  
  featureRanges <- granges(ATAC) # Feature ranges
  
  if(any(Matrix::rowSums(RNA)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNA) != 0
    RNA <- RNA[genesToKeep,]
  }
  
  cat("Number of features in ATAC data:",nrow(ATAC),"\n")
  cat("Number of genes in RNA data:",nrow(RNA),"\n")
  
  if (!genome %in% c("hg19", "hg38", "mm10")) 
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- FigR::hg19TSSRanges
  }, hg38 = {
    TSSg <- FigR::hg38TSSRanges
  }, mm10 = {
    TSSg <- FigR::mm10TSSRanges
  })
  
  # Keep genes that have annotation 
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNA))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNA <- RNA[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg, 
                                   width = windowPadSize,
                                   both = TRUE)
  
  # Get feature summit
  cat("\nTaking feature summits from feature windows ..\n")
  featureSummits <- resize(featureRanges,width = 1,fix = "center")
  
  # Find overlap of all features to all genes given window
  # Subject is Features, query is Gene
  cat("Finding overlapping feature-gene pairs ..\n")
  geneFeatureOv <- findOverlaps(query = TSSflank,subject = featureSummits)
  numPairs <- length(geneFeatureOv)
  
  cat("Found ",numPairs,"total gene-feature pairs for given TSS window ..\n")
  
  cat("Number of feature summits that overlap any gene TSS window: ",length(unique(subjectHits(geneFeatureOv))),"\n")
  cat("Number of gene TSS windows that overlap any feature summit: ",length(unique(queryHits(geneFeatureOv))),"\n\n")
  
  # For each gene, determine observed correlation of each overlapping feature to its associated gene (gene expression)
  
  # For each of those genes, also determine correlation based on background features (run in parallel) and save
  # Do this in parallel, and get data frame of gene-feature-pearson values
  # Fetch background features for each feature tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background features ..\n")
  
  if(is.null(rowData(ATAC)$bias)){
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    
    ATAC <- chromVAR::addGCBias(ATAC, genome=myGenome) }
  
  # Impute NA GC bias values with the mean of the rest of the values
  bias <- ATAC@rowRanges$bias
  ATAC@rowRanges$bias[is.na(bias)] <- mean(bias[!is.na(bias)])
  
  cat("Using ",nBg," iterations ..\n\n")
  
  cat("Computing gene-feature correlations ..\n")
  
  # Randomly select background feature-gene pairs
  bgPairGenes <- sample(1:length(genesToKeep), numBgPairs, replace = T)
  bgPairFeatures <- sample(1:length(featureSummits), numBgPairs, replace = T)
  
  # Get GC content, signal of background features. Also get expression level of background genes.
  bgPairProperties <- data.frame(GC = ATAC@rowRanges$bias[bgPairFeatures],
                                 signal = rowMeans(ATACmat)[bgPairFeatures],
                                 expression = rowMeans(RNA)[bgPairGenes])
  
  # Also get GC content, signal, expression levels for observed feature-gene pairs
  obPairGenes <- geneFeatureOv@from
  obPairFeatures <- geneFeatureOv@to
  obPairProperties <- data.frame(GC = ATAC@rowRanges$bias[obPairFeatures],
                                 signal = rowMeans(ATACmat)[obPairFeatures],
                                 expression = rowMeans(RNA)[obPairGenes])
  
  # Rescale these features so GC/signal/expression roughly have the same weight
  allPairProperties <- scale(rbind(as.matrix(bgPairProperties), as.matrix(obPairProperties)))
  bgPairProperties <- allPairProperties[1:dim(bgPairProperties)[1],]
  obPairProperties <- allPairProperties[(dim(bgPairProperties)[1] + 1):
                                          (dim(bgPairProperties)[1] + dim(obPairProperties)[1]),]
  
  # Find all background pairs of the observed pairs by searching for nearest neighbor pairs in GC-signal-expression space
  bgPairsInds <- FNN::get.knnx(data = bgPairProperties, query = obPairProperties, k=nBg)$nn.index
  
  # Compute the correlation for all background pairs
  pairCorrs <- list()
  rm(ATAC)
  gc()
  for(pairs in c("bg", "ob")){
    
    # Divide the list of feature-gene pairs in to chunks
    if(pairs == "bg") {
      numPairs <- length(bgPairGenes)
      cat("Computing background feature-gene pair correlations")
    } else if(pairs == "ob") {
      numPairs <- length(obPairGenes)
      cat("Computing observed feature-gene pair correlations")
    }
    starts <- seq(1, numPairs, pairsPerChunk)
    ends <- starts  + pairsPerChunk -1
    ends[length(ends)] <- numPairs
    chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)
    
    # Parallelized calculation of feature-gene correlation for each chunk
    if(pairs == "bg") {
      corPairs <- data.frame(Gene = bgPairGenes, Feature = bgPairFeatures)
    } else if(pairs == "ob") {
      corPairs <- data.frame(Gene = obPairGenes, Feature = obPairFeatures)
    }
    
    library(doParallel)
    if(nCores > 1) message("Using ",nCores," cores ..\n")
    opts <- list()
    pb <- txtProgressBar(min = 0, max = length(chunkList), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
    cl <- parallel::makeCluster(nCores)
    clusterEvalQ(cl, .libPaths())
    doSNOW::registerDoSNOW(cl)
    
    corList <- foreach(x=1:length(chunkList),
                       .options.snow = opts, 
                       .export = c(".chunkCore"),
                       .packages = c("dplyr","Matrix")) %dopar%  {
                         corrs <- .chunkCore(chunk=chunkList[[x]],A=ATACmat,R=RNA,O=corPairs,met=metric)
                         return(corrs)
                       }
    
    stopCluster(cl)
    
    if(any(unlist(sapply(corList,is.null)))){
      message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
      message("Please check to see you have enough cores/memory allocated")
      message("Also make sure you have filtered down to non-zero features/genes")
    }
    
    pairCorrs[[pairs]] <- unlist(corList)
  }
  
  # Get correlation pvalues by comparing observes correlation to background distribution
  cat("Calculating pvalues based on background distribution")
  pvals <- pbmcapply::pbmcmapply(
    function(pair_ind){
      obcor <- pairCorrs[["ob"]][pair_ind]
      bgCorrs <- pairCorrs[["bg"]][bgPairsInds[pair_ind,]]
      bgMean <- mean(bgCorrs)
      bgSd <- sd(bgCorrs)
      if(obcor > bgMean){
        pval <- 1-stats::pnorm(q = obcor, mean = bgMean, sd = bgSd)
      }else{
        pval <- stats::pnorm(q = obcor, mean = bgMean, sd = bgSd)
      }
      pval <- pval * 2 # Two-tailed test
    },
    1:length(obPairGenes),
    mc.cores = 8
  )
  
  corrResults <- data.frame(Gene = obPairGenes,
                            Feature = obPairFeatures,
                            rObs = pairCorrs[["ob"]],
                            pvalZ = pvals)
  
  if(!multimapping){
    # Remove multi-mapping features (force 1-1 mapping)
    cat("Keeping max correlation for multi-mapping features ..\n")
    corrResults <- corrResults %>% group_by(Feature) %>% filter(rObs==max(abs(rObs)))
  }
  
  # Swap gene number for gene symbol from TSS annotation lookup
  corrResults$Gene <- as.character(TSSg$gene_name)[corrResults$Gene]
  
  # Swap feature numbers to match reference input feature numbers
  # This only changes if some features had zero signal and were filtered out
  # Use rownames from reference
  corrResults$Feature <- as.numeric(splitAndFetch(rownames(ATACmat)[corrResults$Feature],"Feature",2))
  
  cat("\nFinished!\n")
  
  # If there is a significance cutoff, only keep pairs that are significant
  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    corrResults <- corrResults[corrResults$pvalZ <= p.cut,] # Subset to significant correlations only
  }
  
  corrResults
}

.chunkCore <- function(chunk, # Chunk start and end indices in O
                       A, # ATAC matrix
                       R, # RNA matrix
                       O, # Gene-Feature overlap pairing data.frame (2 columns)
                       met # Correlation method ("spearman" or "pearson")
){
  
  # Get indices of genes and features from overlap object for chunk
  # Assumes query hits are genes and subject hits are features in the overlap object
  geneIndices <- O$Gene[chunk[1]:chunk[2]]
  featureIndices <- O$Feature[chunk[1]:chunk[2]]
  
  pairNames <- cbind(rownames(A)[featureIndices],rownames(R)[geneIndices])
  
  uniqueGenes <- unique(geneIndices)
  uniqueFeatures <- unique(featureIndices)
  
  M1 <- as.matrix(t(A[uniqueFeatures,,drop=FALSE])) # In case only 1 match, keep matrix structure
  M2 <- as.matrix(t(R[uniqueGenes,,drop=FALSE])) # In case only 1 match, keep matrix structure
  
  # Feature x Gene correlation matrix, subset by feature-gene pair names to get corresponding correlation vector
  # NOTE: This correlation call fails if you have maps with just 1 gene / feature. This is unlikely for large chunk sizes
  cor(x = M1,y = M2,method = met)[pairNames]
  
}

