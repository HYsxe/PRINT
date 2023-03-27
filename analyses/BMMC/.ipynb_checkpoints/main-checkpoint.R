# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getGroupData.R")
source("../../code/footprintTracking.R")
source("../../code/getTFBS.R")

# Initialize a footprintingProject object
projectName <- "BMMC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regions <- readRDS(paste0(projectDataDir, "regionRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getRegionBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each pseudo-bulk
pathToFragGrouping <- paste0(projectDataDir, "barcodeGrouping.txt")
barcodeGroups <- read.table(pathToFragGrouping, header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectMainDir, "data/", projectName, "/merged.fragments.tsv.gz")
getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = F)

# Load dispersion model for footprinting
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- 
    readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Get region-by-pseudobulk ATAC matrix 
pseudobulkATACPath <- paste0(projectDataDir, "pseudobulkATAC.rds")
if(!file.exists(pseudobulkATACPath)){
  groupATAC(project) <- getGroupATAC(project)
  groupATAC(project) <- centerCounts(groupATAC(project))
  saveRDS(groupATAC(project), pseudobulkATACPath)
}else{
  groupATAC(project) <- readRDS(pseudobulkATACPath)
}

# Get gene-by-pseudobulk RNA matrix, followed by LogNormalization
pseudobulkRNAPath <- paste0(projectDataDir, "pseudobulkRNA.rds")
if(!file.exists(pseudobulkRNAPath)){
  scRNA <- readRDS(paste0(projectDataDir, "scRNA.rds"))
  RNAMatrix <- scRNA@assays$RNA@data
  pseudobulkRNA <- getGroupRNA(RNAMatrix, barcodeGroups)
  pseudobulkRNA <- t(t(pseudobulkRNA) * 1e6 / colSums(pseudobulkRNA))
  pseudobulkRNA <- log2(pseudobulkRNA + 1)
  groupRNA(project) <- pseudobulkRNA
  saveRDS(pseudobulkRNA, pseudobulkRNAPath)
}else{
  groupRNA(project) <- readRDS(pseudobulkRNAPath)
}

# Load pseudo-bulk metadata
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
groupUMAP(project) <- as.matrix(cbind(groupInfo$UMAP1, groupInfo$UMAP2))
groupCellType(project) <- groupInfo$cellType
groupPseudoTime(project) <- groupInfo$Pseudotime

# Get color mapping for cell types
cellTypes <- unique(groupInfo$cellType)
cellTypeColors <- groupInfo$color[match(cellTypes, groupInfo$cellType)]
names(cellTypeColors) <- cellTypes

# Save the footprintingProject object with loaded data
projectPath <- paste0(projectDataDir, "footprintingProject.rds")
if(!file.exists(projectPath)){
  saveRDS(project, projectPath)
}

# Load TFBS prediction model
h5Path <- "../../data/TFBSPrediction/TFBS_model.h5"
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
TFBindingModel(project) <- TFBSModel

# Specify which programs to run
args <- commandArgs(trailingOnly=TRUE)
visualize <- F
runFootprinting <- F
runTFBS <- T
runNucleosomeTracking <- F

# Visualize footprint
if(visualize){
  
  regionInd <- 107440
  
  cellTypeOrder <- c("NK", "CD4", "CD8",
                     "CLP", "LMPP", "HSC/MPP", "CMP", "MEP", "early-Ery", "late-Ery",
                     "GMP", "CD16mono", "CD14mono", "pDC", "DC", "Baso",
                     "NaiveB", "plasmaB", "MemoryB", "pro/pre-B")
  rowOrder <- order(match(groupCellType(project), cellTypeOrder))
  
  # Visualize predicted TF binding landscape
  plotFeatureHeatmap(project,
                     regionInd = regionInd,
                     feature = "TFBS",
                     cellTypeColors = cellTypeColors,
                     rowOrder = rowOrder
  )
  
  # Visualize nucleosome positions
  plotFeatureHeatmap(project, 
                     regionInd = regionInd, 
                     feature = "footprint",
                     footprintRadius = 50,
                     cellTypeColors = cellTypeColors,
                     rowOrder = rowOrder)
}

# Run footprinting
if(runFootprinting){
  
  # Set footprint radius
  args <- commandArgs(trailingOnly=TRUE)
  if(is.na(args[2])){
    footprintRadius <- 20
  }else{
    footprintRadius <- as.integer(args[2])
  }
  
  project <- getFootprints(
    project, 
    mode = as.character(footprintRadius),
    nCores = 16, 
    footprintRadius = footprintRadius,
    flankRadius = footprintRadius)
  
  system(paste0("mkdir ../../data/", projectName, "/footprints/"))
  saveRDS(footprints(project, as.character(footprintRadius)), 
          paste0(projectDataDir, "footprints/", footprintRadius, "bp_footprints.rds"))
}

# Run TFBS detection
if(runTFBS){
  
  # Get argument from command line
  args <- commandArgs(trailingOnly=TRUE)
  if(is.na(args[2])){
    nChunks <- length(getChunkInterval(regions, regionChunkSize(project))$ends)
    chunkInds <- 1:nChunks
  }else{
    chunkInds <- 10 * (as.integer(args[2]) - 1) + 1:10
  }
  print(paste("Processing chunks", chunkInds[1], "to", chunkInds[length(chunkInds)]))
  
  getTFBS(project,
          chunkInds = chunkInds,
          nCores = 12)
  
}

if(runNucleosomeTracking){
  
  # Track nucleosomes across pseudotime
  LineageInd <- 2 # Select a specific lineage
  lineageProbs <- groupInfo[, c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")]
  lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                              function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})
  lineageGroups <- which(lineageIdentities %in% LineageInd)
  nucleosomeTracks <- footprintTracking(project, 
                                        lineageGroups,
                                        threshold = 2,
                                        pseudoTimeWindowSize = 10,
                                        footprintRadius = 50,
                                        nCores = 12)
  
  saveRDS(nucleosomeTracks, paste0(projectDataDir, "nucleosomeTracks.rds"))
}
