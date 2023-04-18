library(GenomicRanges)
library(dplyr)
library(reticulate)
library(SummarizedExperiment)
library(ggplot2)

setClassUnion("groupsClass", c("integer", "character"))
setClass("footprintingProject", 
         slots = c(
           projectName = "character", # Name of the project
           dataDir = "character", # Where we store input and output data
           mainDir = "character", # The main folder containing the code/analyses/data subfolders
           fragFile = "character", # Path to the ATAC fragments file
           refGenome = "character", # Reference genome. Example "hg38", "mm10"
           countTensor = "list", # Sample-by-region-by-position 3D Tensor of ATAC insertions
           regionRanges = "GRanges", # Regions within which we wish to perform footprinting
           barcodeGrouping = "data.frame", # Two-column data.frame (barcode and group)
           groups = "groupsClass", # Names of all samples/pseudobulks. Can be integer or character
           regionChunkSize = "integer", # We sometimes chunk the region list into batches. This is batch size
           regionWidth = "integer", # Width of each region. Should be an integer vector
           footprints = "list", # Footprinting results
           regionBias = "matrix", # Region-by-position predicted Tn5 bias matrix. All regions should be the same size
           dispModel = "list", # List of dispersion models for each window size
           groupATAC = "matrix", # Region-by-sample matrix of ATAC data
           groupRNA = "matrix", # Gene-by-sample matrix of RNA data
           groupCellType = "character", # Cell type of each sample
           groupUMAP = "matrix", # UMAP coordinate of each sample
           groupPseudoTime = "numeric", # Pseudotime assigned to each sample
           ATACTracks = "matrix", # Region-by-position matrix of ATAC data
           TFBindingModel = "list" # Model that uses multi-scale footprints to predict TF binding
         ))

# Constructor fucntion for our main class
footprintingProject <- function(refGenome, # Reference genome. Example "hg38", "mm10"
                                projectName, # Character. Name of the project
                                mainDir = "../",
                                dataDir = "../") {
  
  new("footprintingProject", 
      projectName = projectName,
      refGenome = refGenome,
      mainDir = mainDir,
      dataDir = paste0(dataDir, projectName, "/"),
      regionChunkSize = 2000L)
  
}

# Method to get region width from the footprintingProject object
setGeneric("regionWidth", function(x) standardGeneric("regionWidth"))

setMethod("regionWidth", "footprintingProject", function(x) x@regionWidth)

# Method to set region width in the footprintingProject object
setGeneric("regionWidth<-", function(x, value) standardGeneric("regionWidth<-"))

setMethod("regionWidth<-", "footprintingProject", function(x, value) {
  x@regionWidth <- value
  x
})

# Method to get region ranges from the footprintingProject object
setGeneric("regionRanges", function(x) standardGeneric("regionRanges"))

setMethod("regionRanges", "footprintingProject", function(x) x@regionRanges)

# Method to set region ranges in the footprintingProject object
setGeneric("regionRanges<-", function(x, value) standardGeneric("regionRanges<-"))

setMethod("regionRanges<-", "footprintingProject", function(x, value) {
  x@regionRanges <- value
  x@regionWidth <- width(x@regionRanges)
  x
})

# Method to get barcode grouping from the footprintingProject object
setGeneric("barcodeGrouping", function(x) standardGeneric("barcodeGrouping"))

setMethod("barcodeGrouping", "footprintingProject", function(x) x@barcodeGrouping)

# Method to set barcode grouping in the footprintingProject object
setGeneric("barcodeGrouping<-", function(x, value) standardGeneric("barcodeGrouping<-"))

setMethod("barcodeGrouping<-", "footprintingProject", function(x, value) {
  x@barcodeGrouping <- value
  x
})

# Method to get region chunk size from the footprintingProject object
setGeneric("regionChunkSize", function(x) standardGeneric("regionChunkSize"))

setMethod("regionChunkSize", "footprintingProject", function(x) x@regionChunkSize)

# Method to set region chunk size in the footprintingProject object
setGeneric("regionChunkSize<-", function(x, value) standardGeneric("regionChunkSize<-"))

setMethod("regionChunkSize<-", "footprintingProject", function(x, value) {
  x@regionChunkSize <- value
  x
})

# Method to get group names from the footprintingProject object
setGeneric("groups", function(x) standardGeneric("groups"))

setMethod("groups", "footprintingProject", function(x) x@groups)

# Method to set group names in the footprintingProject object
setGeneric("groups<-", function(x, value) standardGeneric("groups<-"))

setMethod("groups<-", "footprintingProject", function(x, value) {
  x@groups <- value
  x
})

# Method to get the path to the fragment file from the footprintingProject object
setGeneric("fragFile", function(x) standardGeneric("fragFile"))

setMethod("fragFile", "footprintingProject", function(x) x@fragFile)

# Method to set the path to the fragment file in the footprintingProject object
setGeneric("fragFile<-", function(x, value) standardGeneric("fragFile<-"))

setMethod("fragFile<-", "footprintingProject", function(x, value) {
  x@fragFile <- value
  x
})

# Method to get the directory for storing data from the footprintingProject object
setGeneric("dataDir", function(x) standardGeneric("dataDir"))

setMethod("dataDir", "footprintingProject", function(x) x@dataDir)

# Method to set the directory for storing data for the footprintingProject object
setGeneric("dataDir<-", function(x, value) standardGeneric("dataDir<-"))

setMethod("dataDir<-", "footprintingProject", function(x, value) {
  x@dataDir <- value
  if(!dir.exists(value)) system(paste("mkdir", value))
  x
})

# Method to get the main directory
setGeneric("mainDir", function(x) standardGeneric("mainDir"))

setMethod("mainDir", "footprintingProject", function(x) x@mainDir)

# Method to set the main directory
setGeneric("mainDir<-", function(x, value) standardGeneric("mainDir<-"))

setMethod("mainDir<-", "footprintingProject", function(x, value) {
  x@mainDir <- value
  x
})

# Method to get the reference genome name from the footprintingProject object
setGeneric("refGenome", function(x) standardGeneric("refGenome"))

setMethod("refGenome", "footprintingProject", function(x) x@refGenome)

# Method to set the reference genome name in the footprintingProject object
setGeneric("refGenome<-", function(x, value) standardGeneric("refGenome<-"))

setMethod("refGenome<-", "footprintingProject", function(x, value) {
  x@refGenome <- value
  x
})

# Method to get the Tn5 bias stored in the footprintingProject object
setGeneric("regionBias", function(x) standardGeneric("regionBias"))

setMethod("regionBias", "footprintingProject", function(x) x@regionBias)

# Method to set the Tn5 bias in the footprintingProject object
setGeneric("regionBias<-", function(x, value) standardGeneric("regionBias<-"))

setMethod("regionBias<-", "footprintingProject", function(x, value) {
  if(length(regionRanges(project)) == 0){
    stop("Need to set the regionRanges slot first!")
  }else if(dim(value)[1] != length(regionRanges(project))){
    stop("Number of regions in region ranges and region bias don't match")
  }
  x@regionBias <- value
  x
})

# Method to get the ATAC data for each cell group
setGeneric("groupATAC", function(x) standardGeneric("groupATAC"))

setMethod("groupATAC", "footprintingProject", function(x) x@groupATAC)

# Method to set the the ATAC data for each cell group
setGeneric("groupATAC<-", function(x, value) standardGeneric("groupATAC<-"))

setMethod("groupATAC<-", "footprintingProject", function(x, value) {
  if(length(regionRanges(project)) == 0){
    stop("Need to set the regionRanges slot first!")
  }else if(dim(value)[1] != length(regionRanges(project))){
    stop("Number of regions in region ranges and  bias don't match")
  }
  x@groupATAC <- value
  x
})

# Method to get the RNA data for each cell group
setGeneric("groupRNA", function(x) standardGeneric("groupRNA"))

setMethod("groupRNA", "footprintingProject", function(x) x@groupRNA)

# Method to set the the RNA data for each cell group
setGeneric("groupRNA<-", function(x, value) standardGeneric("groupRNA<-"))

setMethod("groupRNA<-", "footprintingProject", function(x, value) {
  x@groupRNA <- value
  x
})

# Method to get the cell type label for each cell group
setGeneric("groupCellType", function(x) standardGeneric("groupCellType"))

setMethod("groupCellType", "footprintingProject", function(x) x@groupCellType)

# Method to set the the cell type label for each cell group
setGeneric("groupCellType<-", function(x, value) standardGeneric("groupCellType<-"))

setMethod("groupCellType<-", "footprintingProject", function(x, value) {
  x@groupCellType <- value
  x
})

# Method to get the UMAP coordinates for each cell group
setGeneric("groupUMAP", function(x) standardGeneric("groupUMAP"))

setMethod("groupUMAP", "footprintingProject", function(x) x@groupUMAP)

# Method to set the the UMAP coordinates for each cell group
setGeneric("groupUMAP<-", function(x, value) standardGeneric("groupUMAP<-"))

setMethod("groupUMAP<-", "footprintingProject", function(x, value) {
  x@groupUMAP <- value
  x
})

# Method to get the pseudotime for each cell group
setGeneric("groupPseudoTime", function(x) standardGeneric("groupPseudoTime"))

setMethod("groupPseudoTime", "footprintingProject", function(x) x@groupPseudoTime)

# Method to set the pseudotime for each cell group
setGeneric("groupPseudoTime<-", function(x, value) standardGeneric("groupPseudoTime<-"))

setMethod("groupPseudoTime<-", "footprintingProject", function(x, value) {
  x@groupPseudoTime <- value
  x
})

# Method to get the ATAC Tracks for each region
setGeneric("ATACTracks", function(x) standardGeneric("ATACTracks"))

setMethod("ATACTracks", "footprintingProject", function(x) x@ATACTracks)

# Method to set the ATAC Tracks for each 
setGeneric("ATACTracks<-", function(x, value) standardGeneric("ATACTracks<-"))

setMethod("ATACTracks<-", "footprintingProject", function(x, value) {
  x@ATACTracks <- value
  x
})

# Prepares cluster for parallel computing using foreach
prep_cluster <- function(len, # Number of elemnts in the iterable list
                         n_cores = 12 # Number of cores to use
){
  library(doParallel)
  opts <- list()
  pb <- txtProgressBar(min = 0, max = len, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(n_cores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  list("opts" = opts, "cl" = cl)
}

# Chunk a vector/list x into chunks. Return starts and ends of chunks
getChunkInterval <- function(x, # Vector or list
                             chunkSize = 2000 # Size of a single chunk
){
  
  chunkSize <- min(length(x), chunkSize)
  nData <- length(x)
  starts <- seq(1, nData, chunkSize)
  ends <- starts + chunkSize - 1
  ends[length(ends)] <- nData
  
  list("starts" = starts,
       "ends" = ends)
}

# Calculate cosine distance among observations
distCosine <- function(data # observation-by-feature matrix
){
  Matrix <- as.matrix(data)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  dist <- as.dist(1 - sim)
  dist[is.na(dist)] <- 1
  dist
}

# Plotting function used for debugging purposes
plt <- function(data){
  ggplot(data = data.frame(x = 1:length(data), y = data)) +
    geom_line(aes(x = x, y = y))
}

# Performs column-normalization on matrix-like data
# After normalization each column will sum up to the same number
centerCounts <- function(obj, # Input data
                         doInChunks = TRUE, # Whether chunk data when processing 
                         chunkSize = 1000 # Number of cells per chunk
){
  if (!class(obj) %in% c("SummarizedExperiment", "RangedSummarizedExperiment", 
                         "dgCMatrix", "dgeMatrix", "Matrix", "matrix")) 
    stop("Supplied object must be either of class SummarizedExperiment, matrix or sparse Matrix ..\n")
  if (ncol(obj) > 10000) 
    doInChunks <- TRUE
  if (doInChunks) {
    cat("Centering counts for cells sequentially in groups of size ", 
        chunkSize, " ..\n\n")
    starts <- seq(1, ncol(obj), chunkSize)
  }
  else {
    starts <- 1
  }
  counts.l <- list()
  for (i in 1:length(starts)) {
    beginning <- starts[i]
    if (i == length(starts)) {
      ending <- ncol(obj)
    }
    else {
      ending <- starts[i] + chunkSize - 1
    }
    cat("Computing centered counts for cells: ", beginning, 
        " to ", ending, "..\n")
    if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
        "SummarizedExperiment") {
      m <- SummarizedExperiment::assay(obj[, beginning:ending])
    }
    else {
      m <- obj[, beginning:ending]
    }
    obsMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    cCounts <- Matrix::t(Matrix::t(m)/obsMeans)
    counts.l[[i]] <- cCounts
    gc()
  }
  cat("Merging results..\n")
  centered.counts <- do.call("cbind", counts.l)
  cat("Done!\n")
  if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
      "SummarizedExperiment") {
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  }
  else {
    return(centered.counts)
  }
}

# Merge overlapping regions in a GRanges object
# The input should have $score attribute so we can merge based on this score
mergeRegions <- function(regions # The GRanges object to be merged. It must have a $score attribute
){
  
  "%ni%" <- Negate("%in%")
  
  # Sort
  regions <- sortSeqlevels(regions)
  regions <- sort(regions)
  
  # Filter regions based on summit score/quantile normalized summit score
  keptRegions <- 1:length(regions)
  while (!(isDisjoint(regions[keptRegions]))) {
    
    # Fast variable access
    chrNames <- as.character(seqnames(regions[keptRegions]))
    starts <- start(regions[keptRegions])
    ends <- end(regions[keptRegions])
    scores <- regions$score
    
    # See if consecutive regions are overlapping
    overlapNext <- intersect(
      which(chrNames[1:(length(keptRegions) - 1)] == chrNames[2:(length(keptRegions))]),
      which(ends[1:(length(keptRegions) - 1)] >= starts[2:(length(keptRegions))] )
    )
    
    # Compare consectuive regions
    overlapPrevious <- overlapNext + 1
    overlapComparison <- scores[keptRegions[overlapPrevious]] > scores[keptRegions[overlapNext]]
    discard <- keptRegions[c(overlapPrevious[!overlapComparison], overlapNext[overlapComparison])]
    keptRegions <- keptRegions[keptRegions %ni% discard]
  }
  regions <- sortSeqlevels(regions); regions <- sort(regions)
  
  regions[keptRegions]
  
}

# Cluster TF motifs based on Jaccard distance of their genome-wide motif matching patterns
clusterMotifs <- function(motifs, # PWMatrixList object. 
                          regions, # GRanges object. Specifies the regions for motif matching
                          genome, # Ref genome. Example: "hg38"
                          nClusters = 100 # Number of clusters
){
  
  # Find matches of the remaining motifs in regions of interest
  print("Finding matches of the remaining motifs in footprint regions")
  motifMatches <- t(pbmcapply::pbmcmapply(
    function(TF){
      as.matrix(assay(motifmatchr::matchMotifs(subject = regions, 
                                               pwms = motifs[TF], 
                                               genome = genome)))
    },
    names(motifs),
    mc.cores = 16
  ))
  motifMatches <- t(as.matrix(motifMatches * 1)) # Convert to a matrix with 1s and 0s
  
  # Calculating simialrity among motifs
  # Similarity = Jaccard index between motif matches of two TFs
  print("Calculating simialrity among motifs")
  nTFs <- dim(motifMatches)[2]
  intersection <- t(motifMatches) %*% motifMatches
  unmatched <- 1 - motifMatches
  union <- array(dim(unmatched)[1], dim = c(nTFs, nTFs)) - t(unmatched) %*% unmatched
  seqSimilarity <- intersection / union
  
  # K-means clustering
  set.seed(42)
  print("K-means clustering of motifs")
  kmClustering <- kmeans(seqSimilarity, nClusters, iter.max = 100, nstart = 5)
  clusterLabels <- kmClustering$cluster
  
  # Use hierarchical clustering to reorder the K-mean clusters
  clusterCenters <- t(sapply(
    sort(unique(clusterLabels)),
    function(cluster){
      colMeans(seqSimilarity[clusterLabels %in% cluster, , drop = F])
    }
  ))
  hclustTree <- hclust(dist(clusterCenters))
  clusterOrder <- hclustTree$order
  reOrderMap <- 1:nClusters
  names(reOrderMap) <- clusterOrder
  clusterLabels <- reOrderMap[as.character(clusterLabels)]
  motifOrder <- order(clusterLabels)
  
  # Summarize motif clustering results into a data.frame
  motifClustering <- as.data.frame(cbind(names(motifs)[motifOrder], clusterLabels[motifOrder]))
  colnames(motifClustering) <- c("TF", "cluster")
  rownames(motifClustering) <- motifClustering$TF
  motifClustering
  
}

# Calculate pathway enrichment using Fisher's test
pathwayEnrichment <- function(fgGenes, # Foreground gene symbols
                              bgGenes, # Background gene symbols, must contain all genes in the foreground
                              geneSets, # Gene sets. A named list of genes in each pathway
                              pvalThrshold = 0.01, # Threshold to filter results
                              nGeneThreshold = 10 # Remove pathways with foreground genes fewer than this number
){
  
  # Remove repetitive entries from gene lists
  fgGenes <- unique(fgGenes)
  bgGenes <- unique(bgGenes)
  
  # Calculate pathway enrichment using Fisher's exact test
  enrichment <- data.table::rbindlist(pbmcapply::pbmclapply(
    names(geneSets),
    function(pathway){
      isPathway <- bgGenes %in% geneSets[[pathway]]
      isForeground <- bgGenes %in% fgGenes
      contingencyTable <- c(sum(isForeground & isPathway),
                            sum(isForeground & (!isPathway)),
                            sum((!isForeground) & isPathway),
                            sum((!isForeground) & (!isPathway)))
      contingencyTable  <- array(contingencyTable , dim = c(2,2))
      oddsRatio <- (contingencyTable[1,1] / contingencyTable[1,2] ) * (contingencyTable[2,2] / contingencyTable[2,1])
      data.frame(pathway = pathway, 
                 pval = fisher.test(contingencyTable)$p.value, 
                 logOR = log2(oddsRatio),
                 nPathway = sum(isPathway),
                 nDiff = sum(isPathway & isForeground),
                 expected = sum(isPathway) * sum(isForeground) / length(bgGenes),
                 genes = paste(sort(bgGenes[isForeground & isPathway]), collapse = ","))
    },
    mc.cores = 16
  ))
  enrichment <- enrichment[(enrichment$logOR > 0) & (enrichment$nDiff > nGeneThreshold), ]
  enrichment$fdr <- p.adjust(enrichment$pval, method = "fdr")
  enrichment <- enrichment[(enrichment$pval < pvalThrshold), ]
  enrichment <- enrichment[order(enrichment$pval),]
  enrichment
}

# Perform t-test for each feature in the two feature-by-observation matrices x and y
# Here we are using unequal variance Welch test
# For instructions see https://www.statsdirect.co.uk/help/parametric_methods/utt.htm
twoSampTTest <- function(x, # Feature-by-observation matrices x 
                         y # Feature-by-observation matrices y
){
  if(is.null(dim(x)) & is.null(dim(y))){
    n1 <- length(x)
    n2 <- length(y)
    s1 <- sum((x - mean(x)) ^ 2) / (n1 - 1)
    s2 <- sum((y - mean(y)) ^ 2) / (n2 - 1)
    statistic <- (mean(x) - mean(y)) /  sqrt(s1 / n1 + s2 / n2)
  }else{
    n1 <- dim(x)[2]
    n2 <- dim(y)[2]
    s1 <- rowSums((x - rowMeans(x)) ^ 2) / (n1 - 1)
    s2 <- rowSums((y - rowMeans(y)) ^ 2) / (n2 - 1)
    statistic <- (rowMeans(x) - rowMeans(y)) /  sqrt(s1 / n1 + s2 / n2)
  }
  df <- (s1 / n1 + s2 / n2) ^ 2 / ((s1 / n1) ^ 2 / (n1 - 1) + (s2 / n2) ^ 2 / (n2 - 1))
  statistic[is.na(statistic)] <- 0
  
  pvals <- pt(statistic, df = df)
  pvals[is.na(pvals)] <- 0.5
  pvals[pvals > 0.5] <- 1 - pvals[pvals > 0.5]
  pvals <- pvals * 2 # Two tailed test
  
  pvals
}

# Helper function for making density plot
get_density <- function(x, y, n = 100) {
  library(ggpointdensity)
  dens <- MASS::kde2d(x = x, y = y, n = n, 
                      h = c((max(x) - min(x)) / 50,
                            (max(y) - min(y)) / 50))
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

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

getCountsFromFrags <- function(fragFile, peaks, barcodeList = NULL, maxFragLength = NULL, 
                               addColData = TRUE) 
{
  start_time <- Sys.time()
  if (class(fragFile) == "character") {
    GA <- fragsToRanges(fragFile, barcodeList = barcodeList, 
                        startsAre0based = TRUE)
  }
  else if (class(fragFile) == "GRanges") {
    GA <- fragFile
    if (!is.null(barcodeList)) {
      cat("Retaining only select barcodes specified within list ..\\n")
      GA <- GA[mcols(GA)[, 1] %in% barcodeList]
    }
  }
  if (ncol(mcols(GA)) > 1) {
    if (inherits(mcols(GA)[2][, 1], "character")) {
      colnames(mcols(GA)) <- c("barcodeID", "readID")
    }
    else if (inherits(mcols(GA)[2][, 1], "integer")) {
      colnames(mcols(GA)) <- c("barcodeID", "pcrDup")
    }
  }
  else {
    colnames(mcols(GA)) <- "barcodeID"
  }
  if (!is.null(maxFragLength)) {
    cat("Removing frags with length > ", maxFragLength, " bp ..\n")
    GA <- GA[width(GA) <= maxFragLength]
    if (length(GA) == 0) 
      stop("Fragment filtering resulting in 0 aligned fragments. Please check / change the provided filter size ..\n")
  }
  barcodes <- as.character(GA$barcodeID)
  denom <- table(barcodes)
  uniqueBarcodes <- names(denom)
  id <- factor(barcodes, levels = uniqueBarcodes)
  cat("Finding overlap between peaks and fragments in data ..\n")
  ovPEAKStarts <- findOverlaps(query = peaks, subject = resize(GA, 
                                                               width = 1, fix = "start"))
  ovPEAKEnds <- findOverlaps(query = peaks, subject = resize(GA, 
                                                             width = 1, fix = "end"))
  cat("Filtering for valid fragment-peak overlaps based on cut site start/end coordinates ..\n")
  validHits <- unique.data.frame(rbind(as.data.frame(ovPEAKStarts), 
                                       as.data.frame(ovPEAKEnds)))
  require(dplyr)
  cat("Generating matrix of counts ..\n")
  countdf <- data.frame(peaks = validHits$queryHits, sample = as.numeric(id)[validHits$subjectHits]) %>% 
    dplyr::group_by(peaks, sample) %>% dplyr::summarise(count = n()) %>% 
    data.matrix()
  m <- Matrix::sparseMatrix(i = c(countdf[, 1], length(peaks)), 
                            j = c(countdf[, 2], length(uniqueBarcodes)), x = c(countdf[, 
                                                                                       3], 0))
  colnames(m) <- uniqueBarcodes
  if (addColData) {
    cat("Computing sample read depth and FRIP ..\n")
    colData <- data.frame(sample = uniqueBarcodes, depth = as.numeric(denom), 
                          FRIP = Matrix::colSums(m)/as.numeric(denom), stringsAsFactors = FALSE)
    stopifnot(all.equal(colData$sample, colnames(m)))
    if (any(colData$FRIP > 1)) 
      warning("One or more barcodes ended up with FRIP score > 1 .. check your fragment file as it may contain some abnormally large fragments that should be removed ..\n")
    cat("Generating SummarizedExperiment object ..\n")
    SE <- SummarizedExperiment(rowRanges = peaks, assays = list(counts = m), 
                               colData = colData)
  }
  else {
    cat("Generating SummarizedExperiment object ..\n")
    SE <- SummarizedExperiment(rowRanges = peaks, assays = list(counts = m))
  }
  cat("Done!\n")
  end_time <- Sys.time()
  cat("Time elapsed: ", end_time - start_time, units(end_time - 
                                                       start_time), " \n\n")
  return(SE)
}

getGeneScoresFromPeaks <- function (SE, geneList = NULL, genome = c("hg19", "hg38", "mm10"), 
                                    TSSwindow = 10000, getWeightsOnly = FALSE) 
{
  if (length(genome) > 1) 
    stop("Must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  if (!genome %in% c("hg19", "hg38", "mm10")) 
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- BuenRTools::hg19TSSRanges
  }, hg38 = {
    TSSg <- BuenRTools::hg38TSSRanges
  }, mm10 = {
    TSSg <- BuenRTools::mm10TSSRanges
  })
  TSSg$gene_name <- as.character(TSSg$gene_name)
  if (!is.null(geneList)) {
    if (!(all(geneList %in% as.character(TSSg$gene_name)))) 
      stop("One or more of the gene names supplied is not present in the annotation provided..\n")
    cat("Running gene-peak mapping for genes:", geneList, 
        sep = "\n")
    cat("........\n")
    TSSg <- TSSg[TSSg$gene_name %in% geneList, ]
  }
  else {
    cat("Running gene-peak mapping for all genes in annotation! (n = ", 
        length(TSSg), ") This is bound to take more time than querying specific markers ..\n", 
        sep = "")
  }
  cat("Using window of: ", TSSwindow, " bp (total) around TSS per gene ..\n")
  TSSflank <- GenomicRanges::flank(TSSg, width = TSSwindow/2, 
                                   both = TRUE)
  if (!all(start(TSSflank) > 0)) {
    cat("WARNING: ", sum(start(TSSflank) < 0), " flanked TSS window(s) found to extend beyond the chromosomal boundary .. Resetting start coordinates will result in uneven bins for these TSSs..\n")
    start(TSSflank)[start(TSSflank) < 0] <- 0
  }
  peakg <- GenomicRanges::granges(SE)
  peakSummits <- GenomicRanges::start(peakg) + GenomicRanges::width(peakg)/2
  GenomicRanges::start(peakg) <- GenomicRanges::end(peakg) <- peakSummits
  TSSPeakOverlap <- GenomicRanges::findOverlaps(TSSflank, peakg)
  cat("Determining peak weights based on exponential inverse distance to TSS ..\n")
  distToTSS <- abs(start(peakg)[subjectHits(TSSPeakOverlap)] - 
                     start(TSSg)[queryHits(TSSPeakOverlap)])
  weights <- exp(-1 * (distToTSS/1000))
  cat("Assembling Peaks x Genes weights matrix ..\n")
  m <- Matrix::sparseMatrix(i = c(subjectHits(TSSPeakOverlap), 
                                  length(peakg)), j = c(queryHits(TSSPeakOverlap), length(TSSflank)), 
                            x = c(weights, 0))
  colnames(m) <- TSSg$gene_name
  if (getWeightsOnly) 
    return(m)
  weights.m <- m[, Matrix::colSums(m) != 0]
  cat("Assembling Gene x Cells scores matrix ..\n")
  geneScoresMat <- Matrix::t(weights.m) %*% SummarizedExperiment::assay(SE)
  cat("Done!\n\n")
  return(geneScoresMat)
}