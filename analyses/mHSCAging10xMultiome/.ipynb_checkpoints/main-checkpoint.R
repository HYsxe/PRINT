# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getGroupData.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getTFBS.R")

# Initialize a footprintingProject object
projectName <- "mHSCAging10xMultiome/"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "mm10")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName)
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regions <- readRDS(paste0(projectMainDir, "data/mHSCAging10xMultiome/regionRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
biasPath <- paste0(projectMainDir, "data/mHSCAging10xMultiome/predBias.rds")
if(file.exists(biasPath)){
  regionBias(project) <- readRDS(biasPath)
}else{
  project <- getRegionBias(project, nCores = 24)
  saveRDS(regionBias(project), biasPath)
}

# Load barcodes for each pseudo-bulk
pathToFragGrouping <- paste0(projectDataDir, "barcodeGrouping.txt")
barcodeGroups <- read.table(pathToFragGrouping, header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectMainDir, "data/mHSCAging10xMultiome/all.frags.filt.tsv.gz")
getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = F)

# Get gene-by-pseudobulk RNA matrix, followed by LogNormalization
pseudobulkRNAPath <- "../../data/mHSCAging10xMultiome/pseudobulkRNA.rds"
if(!file.exists(pseudobulkRNAPath)){
  scRNA <- readRDS("../../data/mHSCAging10xMultiome/scRNA.rds")
  RNAMatrix <- scRNA@assays$RNA@data
  pseudobulkRNA <- getGroupRNA(RNAMatrix, barcodeGroups)
  pseudobulkRNA <- preprocessCore::normalize.quantiles(pseudobulkRNA)
  rownames(pseudobulkRNA) <- rownames(RNAMatrix)
  groupRNA(project) <- pseudobulkRNA
  saveRDS(pseudobulkRNA, pseudobulkRNAPath)
}else{
  groupRNA(project) <- readRDS(pseudobulkRNAPath)
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

# Load Tn5 background dispersion model
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- 
    readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
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

# Save the project object with pre-loaded slots
saveRDS(project, "../../data/mHSCAging10xMultiome/project.rds")

# Specify which programs to run
args <- commandArgs(trailingOnly=TRUE)
options <- rep(F, 2)
options[as.integer(args[1])] <- T
runTFBS <- options[1]
runFootprinting <- options[2]

# Run footprinting
if(runFootprinting){
  
  # Set footprint radius
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
  if(is.na(args[2])){
    nChunks <- length(getChunkInterval(regions, regionChunkSize(project))$ends)
    chunkInds <- 1:nChunks
  }else{
    chunkInds <- 10 * (as.integer(args[2]) - 1) + 1:10
  }
  print(paste("Processing chunks", chunkInds[1], "to", chunkInds[length(chunkInds)]))
  
  project <- getTFBS(project,
                     chunkInds = chunkInds,
                     nCores = 12)
  
}


